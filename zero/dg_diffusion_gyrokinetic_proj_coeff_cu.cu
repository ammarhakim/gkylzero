/* -*- c++ -*- */

extern "C" {
#include <gkyl_dg_diffusion_gyrokinetic_proj_coeff.h>
#include <gkyl_dg_diffusion_gyrokinetic_proj_coeff_priv.h>
#include <gkyl_array_ops.h>
#include <float.h>
}

// CUDA kernel to set device pointers to kernels.
__global__ static void
gkyl_dg_diff_gk_projC_set_cu_ker_ptrs(struct gkyl_dg_diffusion_gyrokinetic_proj_coeff_kernels *kernels,
  int cdim, struct gkyl_basis pbasis, int dirs_linidx)
{
  int pdim = pbasis.ndim;
  enum gkyl_basis_type basis_type = pbasis.b_type;
  int poly_order = pbasis.poly_order;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->projC = ser_diff_proj_coeff_kernels[pdim-2].list[poly_order-1].kernels[dirs_linidx];
      break;
    default:
      assert(false);
  }
};

void
dg_diff_gk_projC_choose_kernel_cu(struct gkyl_dg_diffusion_gyrokinetic_proj_coeff_kernels *kernels,
  int cdim, struct gkyl_basis pbasis, const bool *diff_in_dir)
{
  int dirs_linidx = dg_diff_gk_projC_diffdirs_linidx(diff_in_dir, cdim);
  gkyl_dg_diff_gk_projC_set_cu_ker_ptrs<<<1,1>>>(kernels, cdim, pbasis, dirs_linidx);
}

__global__ static void
gkyl_dg_diffusion_gyrokinetic_proj_coeff_advance_cu_ker(int cdim, double mass, double vtsq_min, struct gkyl_rect_grid grid,
  const double *nu, const double *xi, struct gkyl_dg_diffusion_gyrokinetic_proj_coeff_kernels *kernels,
  const struct gkyl_range conf_rng, const struct gkyl_range vel_rng, const struct gkyl_range phase_rng,
  const struct gkyl_array *GKYL_RESTRICT gijJ, struct gkyl_array *GKYL_RESTRICT vmap,
  struct gkyl_array *GKYL_RESTRICT vmapSq, struct gkyl_array *GKYL_RESTRICT bmag, struct gkyl_array *GKYL_RESTRICT vtsq,
  struct gkyl_array *GKYL_RESTRICT out)
{
  int idx[GKYL_MAX_DIM];
  int idx_vel[2];
  double xc[GKYL_MAX_DIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_rng.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_rng, tid, idx);

    gkyl_rect_grid_cell_center(&grid, idx, xc);

    for (int d=cdim; d<phase_rng.ndim; d++) idx_vel[d-cdim] = idx[d];

    long loc_conf = gkyl_range_idx(&conf_rng, idx);
    long loc_vel = gkyl_range_idx(&vel_rng, idx_vel);
    long loc_phase = gkyl_range_idx(&phase_rng, idx);

    const double *bmag_c = (const double*) gkyl_array_cfetch(bmag, loc_conf);
    const double *gijJ_c = (const double*) gkyl_array_cfetch(gijJ, loc_conf);
    const double *vtsq_c = (const double*) gkyl_array_cfetch(vtsq, loc_conf);
    const double *vmap_c = (const double*) gkyl_array_cfetch(vmap, loc_vel);
    const double *vmapSq_c = (const double*) gkyl_array_cfetch(vmapSq, loc_vel);
    double *out_c = (double*) gkyl_array_fetch(out, loc_phase);

    kernels->projC(xc, grid.dx, nu, xi, mass, vtsq_min,
      gijJ_c, bmag_c, vtsq_c, vmap_c, vmapSq_c, out_c);

  }
}

void
gkyl_dg_diffusion_gyrokinetic_proj_coeff_advance_cu(gkyl_dg_diffusion_gyrokinetic_proj_coeff* up,
  const struct gkyl_range *conf_rng, const struct gkyl_range *vel_rng, const struct gkyl_range *phase_rng,
  const struct gkyl_array *GKYL_RESTRICT gijJ, struct gkyl_array *GKYL_RESTRICT vmap,
  struct gkyl_array *GKYL_RESTRICT vmapSq, struct gkyl_array *GKYL_RESTRICT bmag, struct gkyl_array *GKYL_RESTRICT vtsq,
  struct gkyl_array *GKYL_RESTRICT out)
{
  int nblocks = phase_rng->nblocks, nthreads = phase_rng->nthreads;

  gkyl_dg_diffusion_gyrokinetic_proj_coeff_advance_cu_ker<<<nblocks, nthreads>>>
    (up->cdim, up->mass, up->vtsq_min, *up->grid, up->nu, up->xi, up->kernels,
     *conf_rng, *vel_rng, *phase_rng, gijJ->on_dev, vmap->on_dev, vmapSq->on_dev,
     bmag->on_dev, vtsq->on_dev, out->on_dev);
}
