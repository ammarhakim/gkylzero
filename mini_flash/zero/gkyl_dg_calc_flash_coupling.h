#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

typedef struct gkyl_dg_calc_flash_coupling gkyl_dg_calc_flash_coupling;

/**
 * Initialize updater to compute flux contributions to Euler equations to 
 * couple PKPM kinetic equations:
 * 
 * Momentum equation flux: div(P)
 * Energy equation flux:   div(P.u) = P : grad(u) + u.div(P)
 * 
 * Free with gkyl_dg_calc_flash_coupling_release
 * 
 * @param conf_grid_kinetic 
 * @param conf_basis
 * @return struct gkyl_dg_calc_flash_coupling* 
 */
struct gkyl_dg_calc_flash_coupling*
gkyl_dg_calc_flash_coupling_new(struct gkyl_rect_grid *phase_grid, 
  struct gkyl_rect_grid *conf_grid_kinetic, const struct gkyl_basis *conf_basis);

/**
 * Compute flux contributions to Euler equations to 
 * couple PKPM kinetic equations:
 * 
 * Momentum equation flux: div(P)
 * Energy equation flux:   div(P.u) = P : grad(u) + u.div(P)
 * 
 * @param up Updater object
 * @param range_kinetic Range over which to compute flux terms. MUST include ghost cells.
 * @param vlasov_pkpm_moms Moments of the distribution function, [rho, p_parallel, p_perp].
 * @param pkpm_prim Primitive moments of the distribution function, 
 *                  [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp, 3*T_xx/m, 3*T_yy/m, 3*T_zz/m].
 * @param pkpm_prim_surf Surface expansions of primitive variables,
 *                       [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 3.0*Txx_xl/m, 3.0*Txx_xr/m, 
 *                        ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 3.0*Tyy_yl/m, 3.0*Tyy_yr/m, 
 *                        ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr, 3.0*Tzz_zl/m, 3.0*Tzz_zr/m] 
 * @param pkpm_div_ppar div(p_parallel b_hat) computed from kinetic flux term for consistent pressure force.
 * @param flash_coupling_fluxes Flux contribution terms, [div(P)-grad(p), (P - pI) : grad(u)].
 */
void
gkyl_dg_calc_flash_coupling_advance(struct gkyl_dg_calc_flash_coupling *up,
  struct gkyl_range *phase_range, struct gkyl_range *conf_range, double mass,
  const struct gkyl_array *fin, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *pkpm_prim, const struct gkyl_array *pkpm_prim_surf, 
  const struct gkyl_array *pkpm_div_ppar, const struct gkyl_array *bvar, 
  const struct gkyl_array *bvar_surf, const struct gkyl_array *max_b,
  struct gkyl_array *flash_coupling_fluxes);

/**
 * Release memory allocated to updater object.
 * 
 * @param up Updater to free
 */
void
gkyl_dg_calc_flash_coupling_release(struct gkyl_dg_calc_flash_coupling *up);