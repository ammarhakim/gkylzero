#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Object type
typedef struct gkyl_gk_neut_fluid_prim_vars gkyl_gk_neut_fluid_prim_vars;

/**
 * Create new updater to compute primitive variables for the GK neutral fluid
 * model. Methods compute:
 *  - udrift (ux, uy, uz).
 *  - pressure p = (gas_gamma - 1)*(E - 1/2 rho u^2).
 *  - udrift and pressure.
 *
 * @param gas_gammaAdiabatic index.
 * @param cbasis Configuration space basis functions
 * @param mem_range Configuration space range that sets the size of the bin_op memory
 *                  for computing primitive moments. Note range is stored so
 *                  updater loops over consistent range for primitive moments
 * @param use_gpu Whether to run on the GPU.
 * @return New updater pointer.
 */
struct gkyl_gk_neut_fluid_prim_vars*
gkyl_gk_neut_fluid_prim_vars_new(double gas_gamma, const struct gkyl_basis* cbasis,
  const struct gkyl_range *mem_range, bool use_gpu);

/**
 * Compute the drift velocity vector (ux, uy, uz).
 *
 * @param up Updater to run.
 * @param moms Fluid moments (rho, rho*ux, rho*uy, rho*uz, totalE).
 * @param udrift Ouput drift velocity.
 */
void gkyl_gk_neut_fluid_prim_vars_udrift_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *udrift);

/**
 * Compute the pressure p = (gas_gamma - 1)*(E - 1/2 rho u^2).
 *
 * @param up Updater to run.
 * @param moms Fluid moments (rho, rho*ux, rho*uy, rho*uz, totalE).
 * @param udrift Ouput drift velocity.
 */
void gkyl_gk_neut_fluid_prim_vars_pressure_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *pressure);

/**
 * Compute the drift velocity vector (ux, uy, uz)
 * and the pressure p = (gas_gamma - 1)*(E - 1/2 rho u^2).
 *
 * @param up Updater to run.
 * @param moms Fluid moments (rho, rho*ux, rho*uy, rho*uz, totalE).
 * @param udrift Ouput drift velocity.
 */
void gkyl_gk_neut_fluid_prim_vars_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *prim_vars);

/**
 * Free memory associated with this updater.
 * @param up Updater to free.
 */
void gkyl_gk_neut_fluid_prim_vars_release(gkyl_gk_neut_fluid_prim_vars *up);
