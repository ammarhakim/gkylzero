/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>    
#include <gkyl_array_ops.h>
#include <gkyl_ghost_surf_calc.h>
#include <gkyl_ghost_surf_calc_priv.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_ghost_surf_calc_advance_cu_ker(const gkyl_ghost_surf_calc* gcalc,
  int dir, int edge, struct gkyl_range edge_rng,
  const struct gkyl_array* fIn, struct gkyl_array* rhs)
{
  double xcc[GKYL_MAX_DIM], xcl[GKYL_MAX_DIM], xcr[GKYL_MAX_DIM];
  int idxc[GKYL_MAX_DIM], idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  double* fBlank;
  fBlank = (double*) gkyl_array_cfetch(fIn, 0);
  for (int i=0; i<fIn->ncomp; ++i) {
    fBlank[i] = 0.0;
  }
  
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < edge_rng.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&edge_rng, tid, idxc);
    gkyl_copy_int_arr(edge_rng.ndim, idxc, idxl);
    gkyl_copy_int_arr(edge_rng.ndim, idxc, idxr);
    idxl[dir] = idxl[dir] - 1; idxr[dir] = idxr[dir] + 1;

    gkyl_rect_grid_cell_center(&gcalc->grid, idxc, xcc);
    gkyl_rect_grid_cell_center(&gcalc->grid, idxl, xcl);
    gkyl_rect_grid_cell_center(&gcalc->grid, idxr, xcr);

    long lincP = gkyl_range_idx(&edge_rng, idxc);
    long linlP = gkyl_range_idx(&edge_rng, idxl);
    long linrP = gkyl_range_idx(&edge_rng, idxr);

    const double* fcptr = (const double*) gkyl_array_cfetch(fIn, lincP);
    
    const double* flptr;
    const double* frptr;
    // const double* fcptr;

    if (edge == 1) {
      flptr = (const double*) gkyl_array_cfetch(fIn, linlP);
      frptr = (const double*) fBlank;
    } else {
      flptr = (const double*) fBlank;
      frptr = (const double*) gkyl_array_cfetch(fIn, linrP);
    }

    // reduce local f to local mom
    gcalc->equation->surf_term(gcalc->equation,
      dir, xcl, xcc, xcr, gcalc->grid.dx, gcalc->grid.dx, gcalc->grid.dx,
      idxl, idxc, idxr, flptr, (const double*) fBlank, frptr, (double*) gkyl_array_fetch(rhs, lincP)
    );
  }
}

void
gkyl_ghost_surf_calc_advance_cu(gkyl_ghost_surf_calc *gcalc,
  const struct gkyl_range *phase_rng,
  const struct gkyl_array *fIn, struct gkyl_array *rhs)
{
  struct gkyl_range edge_rng;
  int nblocks, nthreads;
  int edge;
  int clower_idx[GKYL_MAX_DIM], cupper_idx[GKYL_MAX_DIM] = { 0 };
  for (int dim=0; dim<phase_rng->ndim; ++dim) {
    clower_idx[dim] = phase_rng->lower[dim];
    cupper_idx[dim] = phase_rng->upper[dim];
  }
  
  for(int dir=0; dir<gcalc->cdim; ++dir) {
    edge = 1;
    clower_idx[dir] = phase_rng->upper[dir];
    cupper_idx[dir] = phase_rng->upper[dir];
    gkyl_sub_range_init(&edge_rng, phase_rng, clower_idx, cupper_idx);
    nblocks = edge_rng.nblocks;
    nthreads = edge_rng.nthreads;

    gkyl_ghost_surf_calc_advance_cu_ker<<<nblocks, nthreads>>>(gcalc->on_dev, dir,
      edge, edge_rng, fIn->on_dev, rhs->on_dev);

    edge = 0;
    clower_idx[dir] = phase_rng->lower[dir];
    cupper_idx[dir] = phase_rng->lower[dir];
    gkyl_sub_range_init(&edge_rng, phase_rng, clower_idx, cupper_idx);
    nblocks = edge_rng.nblocks;
    nthreads = edge_rng.nthreads;

    gkyl_ghost_surf_calc_advance_cu_ker<<<nblocks, nthreads>>>(gcalc->on_dev, dir,
      edge, edge_rng, fIn->on_dev, rhs->on_dev);

    // Reset indices for loop over each velocity dimension
    clower_idx[dir] = phase_rng->lower[dir];
    cupper_idx[dir] = phase_rng->upper[dir];
  }
}

gkyl_ghost_surf_calc*
gkyl_ghost_surf_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_dg_eqn *equation, int cdim)
{
  gkyl_ghost_surf_calc *up = (gkyl_ghost_surf_calc*) gkyl_malloc(sizeof(gkyl_ghost_surf_calc));
  up->grid = *grid;
  up->cdim = cdim;
  
  struct gkyl_dg_eqn *eqn = gkyl_dg_eqn_acquire(equation);
  up->equation = eqn->on_dev;
  
  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  gkyl_ghost_surf_calc *up_cu = (gkyl_ghost_surf_calc*) gkyl_cu_malloc(sizeof(gkyl_ghost_surf_calc));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_ghost_surf_calc), GKYL_CU_MEMCPY_H2D);
  up->on_dev = up_cu;
  
  up->equation = eqn;
  
  return up;
}
