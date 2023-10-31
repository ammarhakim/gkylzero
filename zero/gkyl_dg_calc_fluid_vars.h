#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

// Object type
typedef struct gkyl_dg_calc_fluid_vars gkyl_dg_calc_fluid_vars;

/**
 * Create new updater to compute fluid variables needed in 
 * update of DG fluid equations (e.g., isothermal Euler, Euler, 10 moment). Methods compute:
 * p : Isothermal Euler -> p = vth*rho; Euler -> p = (gas_gamma - 1)*(E - 1/2 rho u^2)
 * p_surf : [p_xl, p_xr, p_yl, p_yr, p_zl, p_zr]
 * u : [ux, uy, uz]
 * u_surf : [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
 *           ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 
 *           ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr] 
 * 
 * Updater also stores the kernels to compute fluid source terms and fluid integrated moments.
 * 
 * @param conf_grid Configuration space grid (for getting cell spacing and cell center)
 * @param cbasis Configuration space basis functions
 * @param mem_range Configuration space range that sets the size of the bin_op memory
 *                  for computing primitive moments. Note range is stored so 
 *                  updater loops over consistent range for primitive moments
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_fluid_vars* 
gkyl_dg_calc_fluid_vars_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_range *mem_range, 
  bool use_gpu);

/**
 * Create new updater to compute fluid variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_fluid_vars* 
gkyl_dg_calc_fluid_vars_cu_dev_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_range *mem_range);

/**
 * Compute flow velocity from mass density and momentum density.
 *
 * @param up Updater for computing fluid variables 
 * @param fluid Input array of fluid variables [rho, rho ux, rho uy, rho uz, ...]
 * @param cell_avg_prim Array for storing boolean value of whether rho is negative at corners 
 *                      Note: Only used for diagnostic purposes (not for adjusting solution)
 * @param u      Output array of volume expansion of flow velocity [ux, uy, uz]
 * @param u_surf Output array of surface expansion of flow velocity
 *                  [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
 *                   ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 
 *                   ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr] 
 */
void gkyl_dg_calc_fluid_vars_advance(struct gkyl_dg_calc_fluid_vars *up, const struct gkyl_array* fluid, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* u, struct gkyl_array* u_surf);

/**
 * Compute pressure from fluid variables in the volume and at needed surfaces
 *
 * @param up Updater for computing fluid variables 
 * @param param Input parameter needed for computing pressure (vth for isothermal Euler, gas_gamma for Euler)
 * @param conf_range Configuration space range
 * @param fluid  Input array of fluid variables [rho, rho ux, rho uy, rho uz, ...]
 * @param u      Input array of volume expansion of flow velocity [ux, uy, uz]
 * @param p      Output array of volume expansion of pressure 
 * @param p_surf Output array of surface expansion of pressure [p_xl, p_xr, p_yl, p_yr, p_zl, p_zr] 
 */
void gkyl_dg_calc_fluid_vars_pressure(struct gkyl_dg_calc_fluid_vars *up, 
  double param, const struct gkyl_range *conf_range, 
  const struct gkyl_array* fluid, const struct gkyl_array* u, 
  struct gkyl_array* p, struct gkyl_array* p_surf);

/**
 * Limit slopes for fluid variables
 *
 * @param up Updater for computing fluid variables 
 * @param param Input parameter needed for computing pressure (vth for isothermal Euler, gas_gamma for Euler)
 * @param conf_range Configuration space range
 * @param p      Input array of volume expansion of pressure 
 * @param fluid  Input (and Output after limiting) array of fluid variables [rho, rho ux, rho uy, rho uz, ...]

 */
void gkyl_dg_calc_fluid_vars_limiter(struct gkyl_dg_calc_fluid_vars *up, 
  double param, const struct gkyl_range *conf_range, 
  struct gkyl_array* p, struct gkyl_array* fluid);

/**
 * Delete pointer to updater to compute fluid variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_fluid_vars_release(struct gkyl_dg_calc_fluid_vars *up);

/**
 * Host-side wrappers for fluid vars operations on device
 */

void gkyl_dg_calc_fluid_vars_advance_cu(struct gkyl_dg_calc_fluid_vars *up, const struct gkyl_array* fluid, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* u, struct gkyl_array* u_surf);

void gkyl_dg_calc_fluid_vars_pressure_cu(struct gkyl_dg_calc_fluid_vars *up, 
  double param, const struct gkyl_range *conf_range, 
  const struct gkyl_array* fluid, const struct gkyl_array* u, 
  struct gkyl_array* p, struct gkyl_array* p_surf);

void gkyl_dg_calc_fluid_vars_limiter_cu(struct gkyl_dg_calc_fluid_vars *up, 
  double param, const struct gkyl_range *conf_range, 
  struct gkyl_array* p, struct gkyl_array* fluid);
