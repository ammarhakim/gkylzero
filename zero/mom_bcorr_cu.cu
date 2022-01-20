/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>    
#include <gkyl_array_ops.h>
#include <gkyl_mom_bcorr.h>
#include <gkyl_mom_bcorr_priv.h>
#include <gkyl_mom_type.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_mom_bcorr_advance_cu_ker(const gkyl_mom_bcorr* bcorr,
  const struct gkyl_range conf_rng, struct gkyl_range vel_rng,
  enum gkyl_vel_edge edge, const struct gkyl_array* fin, struct gkyl_array* out)
{
  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < vel_rng.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&vel_rng, tid, pidx);
    gkyl_rect_grid_cell_center(&bcorr->grid, pidx, xc);

    long lincP = gkyl_range_idx(&vel_rng, pidx);
    const double* fptr = (const double*) gkyl_array_cfetch(fin, lincP);
    double momLocal[96]; // hard-coded to max confBasis.num_basis (3x p=3 Ser) for now.
    for (unsigned int k=0; k<96; ++k)
      momLocal[k] = 0.0;

    // reduce local f to local mom
    bcorr->momt->kernel(bcorr->momt, xc, bcorr->grid.dx, pidx, fptr, &momLocal[0], &edge);

    // get conf-space linear index.
    for (unsigned int k = 0; k < conf_rng.ndim; k++)
      cidx[k] = pidx[k];
    long lincC = gkyl_range_idx(&conf_rng, cidx);

    double* mptr = (double*) gkyl_array_fetch(out, lincC);
    for (unsigned int k = 0; k < out->ncomp; ++k) {
       if (tid < vel_rng.volume)
         atomicAdd(&mptr[k], momLocal[k]);
    }
  }
}

void
gkyl_mom_bcorr_advance_cu(gkyl_mom_bcorr *bcorr,
  struct gkyl_range phase_rng, struct gkyl_range conf_rng,
  const struct gkyl_array *fin, struct gkyl_array *out)
{
  struct gkyl_range vel_rng;
  int nblocks, nthreads;
  int viter_idx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  enum gkyl_vel_edge edge;
  
  gkyl_array_clear_range(out, 0.0, conf_rng);
  
  for(int d=0; d<phase_rng.ndim - conf_rng.ndim; ++d) {
    rem_dir[conf_rng.ndim + d] = 1;

    edge = gkyl_vel_edge(d + GKYL_MAX_CDIM);
    viter_idx[conf_rng.ndim + d] = phase_rng.upper[conf_rng.ndim + d];
    gkyl_range_deflate(&vel_rng, &phase_rng, rem_dir, viter_idx);
    nblocks = vel_rng.nblocks;
    nthreads = vel_rng.nthreads;

    gkyl_mom_bcorr_advance_cu_ker<<<nblocks, nthreads>>>(bcorr->on_dev,
      conf_rng, vel_rng, edge, fin->on_dev, out->on_dev);

    edge = gkyl_vel_edge(d);
    viter_idx[conf_rng.ndim + d] = phase_rng.lower[conf_rng.ndim + d];
    gkyl_range_deflate(&vel_rng, &phase_rng, rem_dir, viter_idx);
    nblocks = vel_rng.nblocks;
    nthreads = vel_rng.nthreads;

    gkyl_mom_bcorr_advance_cu_ker<<<nblocks, nthreads>>>(bcorr->on_dev,
      conf_rng, vel_rng, edge, fin->on_dev, out->on_dev);
    
    rem_dir[conf_rng.ndim + d] = 0;
    viter_idx[conf_rng.ndim + d] = 0;
  }
}

gkyl_mom_bcorr*
gkyl_mom_bcorr_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt)
{
  gkyl_mom_bcorr *up = (gkyl_mom_bcorr*) gkyl_malloc(sizeof(gkyl_mom_bcorr));
  up->grid = *grid;
  
  struct gkyl_mom_type *mt = gkyl_mom_type_acquire(momt);
  up->momt = mt->on_dev;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  gkyl_mom_bcorr *up_cu = (gkyl_mom_bcorr*) gkyl_cu_malloc(sizeof(gkyl_mom_bcorr));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_mom_bcorr), GKYL_CU_MEMCPY_H2D);

  up->momt = mt;
  up->on_dev = up_cu;
  
  return up;
}
