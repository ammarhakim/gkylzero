/* -*- c++ -*- */

extern "C" {
#include <gkyl_bc_basic.h>
#include <gkyl_bc_basic_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
}

__global__ static void
bc_sheath_gyrokinetic_advance_cu_ker(int cdim, int dir, struct gkyl_range skin_r, struct gkyl_range ghost_r,
  const struct gkyl_basis *basis, struct gkyl_rect_grid grid,
  double q2Dm, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_array *distf)
{
  int fidx[GKYL_MAX_DIM]; // Flipped index.
  int pidx[GKYL_MAX_DIM], cidx[3];
  double xc[GKYL_MAX_DIM];

  int vpar_dir = cdim;
  double dvpar = grid.dx[vpar_dir];
  int uplo = skin_r.upper[vpar_dir]+skin_r.lower[vpar_dir];

  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < skin_r.volume; linc += blockDim.x*gridDim.x) {

    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&skin_r, linc, pidx);

    gkyl_copy_int_arr(skin_r.ndim, pidx, fidx);
    fidx[vpar_dir] = uplo - pidx[vpar_dir];
    // Turn this skin fidx into a ghost fidx.
    fidx[up->dir] = ghost_r.lower[up->dir];

    gkyl_rect_grid_cell_center(grid, pidx, xc);
    double vpar_c = xc[vpar_dir];
    double vparAbsSq_lo = vpar_c > 0.? pow(vpar_c-0.5*dvpar,2) : pow(vpar_c+0.5*dvpar,2);
    double vparAbsSq_up = vpar_c > 0.? pow(vpar_c+0.5*dvpar,2) : pow(vpar_c-0.5*dvpar,2);

    long skin_loc = gkyl_range_idx(&skin_r, pidx);
    long ghost_loc = gkyl_range_idx(&ghost_r, fidx);

    const double *inp = (const double*) gkyl_array_cfetch(distf, skin_loc);
    double *out = (double*) gkyl_array_fetch(distf, ghost_loc);

    for (int d=0; d<cdim; d++) cidx[d] = pidx[d];
    long conf_loc = gkyl_range_idx(&conf_r, cidx);
    const double *phi_p = (const double*) gkyl_array_cfetch(phi, conf_loc);
    const double *phi_wall_p = (const double*) gkyl_array_cfetch(phi_wall, conf_loc);

    // Calculate reflected distribution function fhat.
    // note: reflected distribution can be
    // 1) fhat=0 (no reflection, i.e. absorb),
    // 2) fhat=f (full reflection)
    // 3) fhat=c*f (partial reflection)
    up->ker_reflectedf(vpar_c, dvpar, vparAbsSq_lo, vparAbsSq_up, q2Dm, phi_p, phi_wall_p, inp, out);

    // Reflect fhat into skin cells.
    bc_gksheath_reflect(dir, basis, cdim, out, out);
  }
}

void
bc_sheath_gyrokinetic_advance_cu(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_array *distf)
{

  int nblocks = up->skin_r.nblocks, nthreads = up->skin_r.nthreads;

  bc_sheath_gyrokinetic_advance_cu_ker<<<nblocks, nthreads>>>(up->cdim, up->dir, up->skin_r, up->ghost_r,
    up->basis, up->gridi, up->q2Dm, phi, phi_wall, distf);
}
