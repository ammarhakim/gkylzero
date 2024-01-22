#include <gkyl_bc_neutral_recycling.h>
#include <gkyl_bc_neutral_recycling_priv.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_ghost_surf_calc.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_const.h>
#include <gkyl_alloc.h>
#include <assert.h>

struct gkyl_bc_neutral_recycling*
gkyl_bc_neutral_recycling_new(int dir, enum gkyl_edge_loc edge, const struct gkyl_basis *pbasis,
  const struct gkyl_basis *cbasis, const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  const struct gkyl_range *ghost_conf_r, const struct gkyl_rect_grid *grid,
  const struct gkyl_dg_eqn *vlasov_eqn, double recyc_temp, double recyc_frac,
  double ion_mass, bool use_gpu)
{

  // Allocate space for new updater.
  struct gkyl_bc_neutral_recycling *up = gkyl_malloc(sizeof(struct gkyl_bc_neutral_recycling));

  up->dir = dir;
  up->cdim = cdim;
  up->edge = edge;
  up->use_gpu = use_gpu;
  up->pbasis = pbasis;
  up->grid = grid;
  up->skin_r = skin_r;
  up->ghost_r = ghost_r;
  up->ghost_conf_r = ghost_conf_r;
  up->recyc_frac = recyc_frac;

  int cdim = up->cbasis->ndim;
  int pdim = up->pbasis->ndim;
  int vdim = pdim - cdim; 
  
  // Assume dir is same direction as B

  // Define arrays
  struct gkyl_array *recyc_fmax_flux = gkyl_array_new(GKYL_DOUBLE, pbasis->num_basis, up->ghost_r->volume); // different range??
  struct gkyl_array *moms_recyc = gkyl_array_new(GKYL_DOUBLE, (2+vdim)*cbasis->num_basis, up->ghost_conf_r->volume);
  struct gkyl_array *m0_recyc = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->ghost_conf_r->volume);
  struct gkyl_array *m2_recyc = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->ghost_conf_r->volume);
  // Will be passed at advance method
  /* up->recyc_fmax = gkyl_array_new(GKYL_DOUBLE, pbasis->num_basis, up->ghost_r->volume); */
  /* up->recyc_flux_m0 = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->ghost_conf_r->volume); */
  up->scale_fmax_m0 = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, up->ghost_conf_r->volume); 
  
  struct gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_new(up->grid, up->cbasis, up->pbasis, poly_order+1, use_gpu);
  struct gkyl_ghost_surf_calc *flux_slvr = gkyl_ghost_surf_calc_new(grid, vlasov_eqn, cdim, use_gpu);
  struct gkyl_mom_type *m0_flux_t = gkyl_mom_vlasov_new(up->cbasis, up->pbasis, "M0", use_gpu);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(grid, m0_flux_t, false);

  // Define weak division function
  up->mem = gkyl_dg_bin_op_mem_cu_dev_new(up->recyc_flux_m0->size, cbasis.num_basis);
  
  // Create unit-density Maxwellian with recycle temp
  double dg_fac = 1/pow(1/sqrt(2),cdim);
  double m2_val_coef = vdim*recyc_frac*GKYL_ELEMENTARY_CHARGE/ion_mass*dg_fac;
  gkyl_array_shiftc(m0_recyc, dg_fac, 0);
  gkyl_array_shiftc(m2_recyc, m2_val_coef, 0);
  gkyl_array_set_offset_range(moms_recyc, 1., m0_recyc, 0, up->conf_rng);
  gkyl_array_set_offset_range(moms_recyc, 1., m2_recyc, (1+vdim)*up->cbasis->num_basis, up->conf_rng);

  gkyl_proj_maxwellian_on_basis_lab_mom(proj_max, ghost_r, ghost_conf_r, moms_recyc, up->recyc_fmax);

  // Calculate flux and zeroth moment of flux.
  // How to get correct direction? 
  gkyl_ghost_surf_calc_advance(flux_slvr, ghost_r, recyc_fmax, recyc_fmax_flux);
  gkyl_mom_calc_advance(m0calc, ghost_r, ghost_conf_r, recyc_fmax_flux, up->recyc_flux_m0);
  
  return up;
}

// assume advance is only for one boundary at a time
void
gkyl_bc_neutral_recycling_advance(const struct gkyl_bc_neutral_recycling *up, const struct gkyl_array *ion_flux_m0,
  struct gkyl_array *distf, const struct gkyl_range *conf_r)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_bc_neutral_recycling_advance_cu(up, ion_flux_m0, distf, conf_r);
    return;
  }
#endif

  // scale_fmax_m0 = ion_flux_m0 / recyc_flux_m0
  gkyl_dg_div_op_range(up->mem, up->cbasis, 0, up->scale_fmax_m0, 0, ion_flux_m0, 0, up->recyc_flux_m0, up->ghost_conf_r);
  // fhat = scale_fmax_m0 * recyc_fmax
  gkyl_dg_mul_conf_phase_op_range(up->cbasis, up->pbasis, distf, up->scale_fmax_m0, up->recyc_fmax, up->ghost_conf_r, up->ghost_r);

  // distf for ghost is created. now what? flip and copy to buffer? 
  
}

void gkyl_bc_neutral_recycling_release(struct gkyl_bc_neutral_recycling *up)
{
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
  
  }
#endif
  
  gkyl_free(up);
}
