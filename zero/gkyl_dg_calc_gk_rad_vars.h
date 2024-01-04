#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_calc_gk_rad_vars gkyl_dg_calc_gk_rad_vars;

/**
 * Create new updater to compute the drag coefficients needed for 
 * radiation in gyrokinetic equations. Method computes:
 * 1. vnu = 2/pi*|v|*nu(v)
 * 2. vsqnu = 1/2*(m/B)^(3/2)*sqrt(mu)*|v|^2*nu(v)
 *    where |v| = sqrt(v_par^2 + 2 mu B/m)
 *    Note that through the spatial variation of B, both these drag coefficients depend on the full phase space
 * 
 * @param phase_grid Phase space grid (for getting cell spacing and cell center)
 * @param conf_basis Configuration space basis functions
 * @param phase_basis Phase space basis functions
 * @param charge Species charge
 * @param mass Species mass
 * @param gk_geom Geometry struct
 * @param a, alpha, beta, gamma, v0 Input fitting parameters for a given collision type
 * @return New updater pointer.
 */
struct gkyl_dg_calc_gk_rad_vars* 
gkyl_dg_calc_gk_rad_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  double charge, double mass, const struct gk_geometry *gk_geom, 
  double a, double alpha, double beta, double gamma, double v0);

/**
 * Compute drag coefficients needed for radiation in gyrokinetic equations
 * 
 * @param up Updater for computing gyrokinetic radiation variables 
 * @param conf_range Configuration space range (should only be local range because geometry only defined on local range)
 * @param phase_range Phase space range 
 * @param vnu Output volume expansion of vpar drag coefficient, vnu = 2/pi*|v|*nu(v)
 * @param vsqnu Output volume expansion of mu drag coefficient, vsqnu = 1/2*(m/B)^(3/2)*sqrt(mu)*|v|^2*nu(v)
 */
void gkyl_dg_calc_gk_rad_vars_advance(const struct gkyl_dg_calc_gk_rad_vars *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  struct gkyl_array* vnu, struct gkyl_array* vsqnu);

/**
 * Delete pointer to updater to compute gyrokinetic variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_gk_rad_vars_release(struct gkyl_dg_calc_gk_rad_vars *up);
