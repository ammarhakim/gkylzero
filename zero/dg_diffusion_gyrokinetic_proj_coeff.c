#include <gkyl_dg_diffusion_gyrokinetic_proj_coeff.h>
#include <gkyl_dg_diffusion_gyrokinetic_proj_coeff_priv.h>
#include <gkyl_alloc.h>
#include <string.h>

struct gkyl_dg_diffusion_gyrokinetic_proj_coeff*
gkyl_dg_diffusion_gyrokinetic_proj_coeff_new(int cdim, struct gkyl_basis pbasis, struct gkyl_rect_grid *grid,
  const double *D, const double *chi, double mass, double vtsq_min, const bool *diff_in_dir, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_dg_diffusion_gyrokinetic_proj_coeff *up = gkyl_malloc(sizeof(*up));

  up->use_gpu = use_gpu;
  up->cdim = cdim;
  up->pdim = pbasis.ndim;
  up->grid = grid;
  up->mass = mass;
  up->vtsq_min = vtsq_min;

  // Set nu and chi to give the requested D and chi in the d/dx(D*dn/dx)
  // and d/dx(chi*n*dT/dx) terms.
  double nu[3], xi[3];
  for (int d=0; d<cdim; d++) {
    nu[d] = (2.0/7.0)*(5.0*D[d] - chi[d]);
    xi[d] = (2.0/7.0)*(2.0*chi[d] - 3.0*D[d])/(6.0*nu[d]);
  }

  if (!use_gpu) {
    up->nu = gkyl_malloc(cdim*sizeof(double)); 
    up->xi = gkyl_malloc(cdim*sizeof(double)); 
    memcpy(up->nu, nu, cdim*sizeof(double)); 
    memcpy(up->xi, xi, cdim*sizeof(double)); 

    up->kernels = gkyl_malloc(sizeof(struct gkyl_dg_diffusion_gyrokinetic_proj_coeff_kernels));
  }
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->nu = gkyl_cu_malloc(cdim*sizeof(double)); 
    up->xi = gkyl_cu_malloc(cdim*sizeof(double)); 
    gkyl_cu_memcpy(up->nu, nu, cdim*sizeof(double), GKYL_CU_MEMCPY_H2D); 
    gkyl_cu_memcpy(up->xi, xi, cdim*sizeof(double), GKYL_CU_MEMCPY_H2D); 

    up->kernels = gkyl_cu_malloc(sizeof(struct gkyl_dg_diffusion_gyrokinetic_proj_coeff_kernels));
  }
#endif

  // Choose kernels that translates the DG coefficients.
  dg_diff_gk_projC_choose_kernel(up->kernels, cdim, pbasis, diff_in_dir, use_gpu);

  return up;
}

void
gkyl_dg_diffusion_gyrokinetic_proj_coeff_advance(gkyl_dg_diffusion_gyrokinetic_proj_coeff* up,
  const struct gkyl_range *conf_rng, const struct gkyl_range *vel_rng, const struct gkyl_range *phase_rng,
  const struct gkyl_array *GKYL_RESTRICT gijJ, struct gkyl_array *GKYL_RESTRICT vmap,
  struct gkyl_array *GKYL_RESTRICT vmapSq, struct gkyl_array *GKYL_RESTRICT bmag, struct gkyl_array *GKYL_RESTRICT vtsq,
  struct gkyl_array *GKYL_RESTRICT out)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_dg_diffusion_gyrokinetic_proj_coeff_advance_cu(up, conf_rng, vel_rng, phase_rng, gijJ, vmap, vmapSq, bmag, vtsq, out);
    return;
  }
#endif

  int cdim = up->cdim;
  int pdim = up->pdim;

  int idx_vel[2];
  double xc[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_rng);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(up->grid, iter.idx, xc);

    for (int d=cdim; d<pdim; d++) idx_vel[d-cdim] = iter.idx[d];

    long loc_conf = gkyl_range_idx(conf_rng, iter.idx);
    long loc_vel = gkyl_range_idx(vel_rng, idx_vel);
    long loc_phase = gkyl_range_idx(phase_rng, iter.idx);

    const double *bmag_c = gkyl_array_cfetch(bmag, loc_conf);
    const double *gijJ_c = gkyl_array_cfetch(gijJ, loc_conf);
    const double *vtsq_c = gkyl_array_cfetch(vtsq, loc_conf);
    const double *vmap_c = gkyl_array_cfetch(vmap, loc_vel);
    const double *vmapSq_c = gkyl_array_cfetch(vmapSq, loc_vel);
    double *out_c = gkyl_array_fetch(out, loc_phase);

    up->kernels->projC(xc, up->grid->dx, up->nu, up->xi, up->mass, up->vtsq_min,
      gijJ_c, bmag_c, vtsq_c, vmap_c, vmapSq_c, out_c);
  }
}

void
gkyl_dg_diffusion_gyrokinetic_proj_coeff_release(gkyl_dg_diffusion_gyrokinetic_proj_coeff* up)
{
  // Release memory associated with this updater.
  if (!up->use_gpu) {
    gkyl_free(up->nu); 
    gkyl_free(up->xi); 
    gkyl_free(up->kernels);
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->nu); 
    gkyl_cu_free(up->xi); 
    gkyl_cu_free(up->kernels);
  }
#endif
  gkyl_free(up);
}

