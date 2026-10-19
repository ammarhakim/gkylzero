#include <gkyl_dg_calc_flash_coupling.h>
#include <gkyl_dg_calc_flash_coupling_priv.h>

#include <gkyl_array_rio.h>

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

// Initialize updater object
struct gkyl_dg_calc_flash_coupling*
gkyl_dg_calc_flash_coupling_new(struct gkyl_rect_grid *phase_grid,
   struct gkyl_rect_grid *conf_grid_kinetic, const struct gkyl_basis *conf_basis)
{
  struct gkyl_dg_calc_flash_coupling *up = gkyl_malloc(sizeof(struct gkyl_dg_calc_flash_coupling));

  up->num_conf_basis = conf_basis->num_basis;

  up->phase_grid = phase_grid;
  up->conf_grid_kinetic = conf_grid_kinetic;

  up->pdim = phase_grid->ndim;
  up->cdim = up->conf_grid_kinetic->ndim;

  for (int d=0; d<up->cdim; ++d) {
    up->ker.mom_flux_ker[d] = choose_flash_coupling_mom_kernel(d, up->cdim);
    up->ker.ener_flux_ker[d] = choose_flash_coupling_ener_kernel(d, up->cdim);
    up->ker.div_q_flux_ker[d] = choose_flash_coupling_div_q_kernel(d, up->cdim);
  }

  return up;
}

// Compute the flux contributions from the kinetic equation,
// 1) Momentum equation pressure force: div(Pi)
// 2) Energy equation viscous heating:   div(Pi.u) = Pi : grad(u) + u . div(Pi) = pi_gradu + u_divpi.
// 3) Energy equation heat flux: div(q)
//
// p_gradu is computed from quantities used in the kinetic update
// u_divp uses the Euler bulk velocity and div(P) computed directly from the kinetic flux  
void
gkyl_dg_calc_flash_coupling_advance(struct gkyl_dg_calc_flash_coupling *up,
  struct gkyl_range *phase_range, struct gkyl_range *conf_range, double mass,
  const struct gkyl_array *fin,const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *pkpm_prim, const struct gkyl_array *pkpm_prim_surf, 
  const struct gkyl_array *pkpm_div_ppar, const struct gkyl_array *bvar, 
  const struct gkyl_array *bvar_surf, const struct gkyl_array *max_b,
  struct gkyl_array *flash_coupling_fluxes) 
{
  gkyl_array_clear(flash_coupling_fluxes, 0.0);
  int idxc_l[GKYL_MAX_DIM], idxc_r[GKYL_MAX_DIM];
  int idxp_l[GKYL_MAX_DIM], idxp_r[GKYL_MAX_DIM];

  double xmu[GKYL_MAX_DIM];
  int idxp[GKYL_MAX_DIM];

  // Configuration space dimensions to remove to create velocity space range
  int rem_dir[GKYL_MAX_DIM] = {0};
  for (int d=0; d<up->cdim; ++d) rem_dir[d] = 1;

  struct gkyl_range vel_range;

  // Iterate over kinetic grid.
  // Outer loop over configuration space to compute div(Pi) and div(Pi.u)
  // Inner loop over velocity space to compute div(q)
  struct gkyl_range_iter conf_iter, vel_iter;
  gkyl_range_iter_init(&conf_iter, conf_range);
  while (gkyl_range_iter_next(&conf_iter)) {
    gkyl_copy_int_arr(up->cdim, conf_iter.idx, idxc_l);
    gkyl_copy_int_arr(up->cdim, conf_iter.idx, idxc_r);

    for (int d=0; d<up->cdim; ++d) {
      idxc_l[d] -= 1;
      idxc_r[d] += 1;

      long linc_l = gkyl_range_idx(conf_range, idxc_l);
      long linc = gkyl_range_idx(conf_range, conf_iter.idx);
      long linc_r = gkyl_range_idx(conf_range, idxc_r);

      const double *pkpm_prim_c = (const double *)gkyl_array_cfetch(pkpm_prim, linc);
      const double *pkpm_div_ppar_c = (const double *)gkyl_array_cfetch(pkpm_div_ppar, linc);

      // Need vlasov_pkpm_moms and pkpm_prim_surf in 3-cell stencil
      const double *vlasov_pkpm_moms_l = (const double *)gkyl_array_cfetch(vlasov_pkpm_moms, linc_l);
      const double *vlasov_pkpm_moms_c = (const double *)gkyl_array_cfetch(vlasov_pkpm_moms, linc);
      const double *vlasov_pkpm_moms_r = (const double *)gkyl_array_cfetch(vlasov_pkpm_moms, linc_r);
      const double *pkpm_prim_surf_l = (const double *)gkyl_array_cfetch(pkpm_prim_surf, linc_l);
      const double *pkpm_prim_surf_c = (const double *)gkyl_array_cfetch(pkpm_prim_surf, linc);
      const double *pkpm_prim_surf_r = (const double *)gkyl_array_cfetch(pkpm_prim_surf, linc_r);

      double *coupling_fluxes_c = (double *)gkyl_array_fetch(flash_coupling_fluxes, linc);
      double *flux_div_pi_c = &coupling_fluxes_c[0];
      double *flux_pi_gradu_c = &coupling_fluxes_c[up->num_conf_basis];

      // Compute div(Pi)
      up->ker.mom_flux_ker[d](up->phase_grid->dx, 
        vlasov_pkpm_moms_l, vlasov_pkpm_moms_c, vlasov_pkpm_moms_r,
        pkpm_div_ppar_c, flux_div_pi_c);
      
      // Compute Pi:grad(u)
      up->ker.ener_flux_ker[d](up->phase_grid->dx,
      vlasov_pkpm_moms_c, pkpm_prim_c, 
      pkpm_prim_surf_l, pkpm_prim_surf_c, pkpm_prim_surf_r, flux_pi_gradu_c);

      // Iterate in velocity space over velocity directions to compute total div(q)
      gkyl_range_deflate(&vel_range, phase_range, rem_dir, conf_iter.idx);
      gkyl_range_iter_init(&vel_iter, &vel_range);
      while (gkyl_range_iter_next(&vel_iter)) {
        copy_idx_arrays(up->cdim, up->pdim, conf_iter.idx, vel_iter.idx, idxp);
        copy_idx_arrays(up->cdim, up->pdim, idxc_l, vel_iter.idx, idxp_l);
        copy_idx_arrays(up->cdim, up->pdim, idxc_r, vel_iter.idx, idxp_r);

        long linp_l = gkyl_range_idx(phase_range, idxp_l);
        long linp = gkyl_range_idx(phase_range, idxp);
        long linp_r= gkyl_range_idx(phase_range, idxp_r);

        // Get cell-centered physical coordinates in central cell
        gkyl_rect_grid_cell_center(up->phase_grid, idxp, xmu);

        const double *f_l = (const double *)gkyl_array_cfetch(fin, linp_l);
        const double *f_c = (const double *)gkyl_array_cfetch(fin, linp);
        const double *f_r = (const double *)gkyl_array_cfetch(fin, linp_r);
        const double *bvar_c = (const double *)gkyl_array_cfetch(bvar, linc);
        const double *bvar_surf_l = (const double *)gkyl_array_cfetch(bvar_surf, linc_l);
        const double *bvar_surf_c = (const double *)gkyl_array_cfetch(bvar_surf, linc);
        const double *bvar_surf_r = (const double *)gkyl_array_cfetch(bvar_surf, linc_r);
        const double *max_b_c = (const double *)gkyl_array_cfetch(max_b, linc);
        double *flux_div_q_c = &coupling_fluxes_c[2*up->num_conf_basis];

        up->ker.div_q_flux_ker[d](xmu, up->phase_grid->dx, mass, f_l, f_c, f_r, bvar_c, max_b_c, 
          bvar_surf_l, bvar_surf_c, bvar_surf_r, flux_div_q_c);
      }
    }
  }
}

void
gkyl_dg_calc_flash_coupling_release(struct gkyl_dg_calc_flash_coupling *up)
{
  gkyl_free(up);
}