#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Type of prim_vars calculation to perform.
enum gkyl_gk_neut_fluid_prim_vars_type {
  GKYL_GK_NEUT_FLUID_PRIM_VARS_UDRIFT = 0,
  GKYL_GK_NEUT_FLUID_PRIM_VARS_PRESSURE,
  GKYL_GK_NEUT_FLUID_PRIM_VARS_TEMP,
  GKYL_GK_NEUT_FLUID_PRIM_VARS_UDRIFT_PRESSURE, // (ux, uy, uz, p).
  GKYL_GK_NEUT_FLUID_PRIM_VARS_UDRIFT_TEMP, // (ux, uy, uz, T).
  GKYL_GK_NEUT_FLUID_PRIM_VARS_LTE, // // (n, ux, uy, uz, T/m).
};

// Object type
typedef struct gkyl_gk_neut_fluid_prim_vars gkyl_gk_neut_fluid_prim_vars;

/**
 * Create new updater to compute primitive variables for the GK neutral fluid
 * model. Methods compute:
 *  - udrift (ux, uy, uz) = (rho u) / rho.
 *  - pressure p = (gas_gamma - 1)*(E - 1/2 rho u^2).
 *  - temperature T = (gas_gamma - 1)*(mass * E - 1/2 (rho u)^2) / rho.
 *  - udrift and pressure.
 *  - udrift and temperature.
 *
 * @param gas_gamma Adiabatic index.
 * @param mass Species mass.
 * @param cbasis Configuration space basis functions
 * @param mem_range Configuration space range that sets the size of the bin_op memory
 *                  for computing primitive moments. Note range is stored so
 *                  updater loops over consistent range for primitive moments
 * @param use_gpu Whether to run on the GPU.
 * @return New updater pointer.
 */
struct gkyl_gk_neut_fluid_prim_vars*
gkyl_gk_neut_fluid_prim_vars_new(double gas_gamma, double mass, const struct gkyl_basis* cbasis,
  const struct gkyl_range *mem_range, enum gkyl_gk_neut_fluid_prim_vars_type prim_vars_type,
  bool use_gpu);

/**
 * Compute the drift velocity vector (ux, uy, uz).
 *
 * @param up Updater to run.
 * @param moms Fluid moments (rho, rho*ux, rho*uy, rho*uz, totalE).
 * @param out Output drift velocity.
 * @param out_coff Offset in out where to place drift velocity.
 */
void gkyl_gk_neut_fluid_prim_vars_udrift_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff);

/**
 * Compute the pressure p = (gas_gamma - 1)*(E - 1/2 rho u^2).
 *
 * @param up Updater to run.
 * @param out Output pressure.
 * @param out_coff Offset in out where to place pressure.
 */
void gkyl_gk_neut_fluid_prim_vars_pressure_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff);

/**
 * Compute the temperature T = p/n = (gas_gamma - 1)*(mass * E - 1/2 (rho u)^2)/rho.
 *
 * @param up Updater to run.
 * @param moms Fluid moments (rho, rho*ux, rho*uy, rho*uz, totalE).
 * @param out Output temperature.
 * @param out_coff Offset in out where to place temperature.
 */
void gkyl_gk_neut_fluid_prim_vars_temp_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff);

/**
 * Compute the drift velocity vector (ux, uy, uz)
 * and the pressure p = (gas_gamma - 1)*(E - 1/2 rho u^2).
 *
 * @param up Updater to run.
 * @param moms Fluid moments (rho, rho*ux, rho*uy, rho*uz, totalE).
 * @param out Output primitive moments.
 * @param out_coff Offset in out where to place primitive moments.
 */
void gkyl_gk_neut_fluid_prim_vars_udrift_pressure_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff);

/**
 * Compute the drift velocity vector (ux, uy, uz)
 * and the temperature T = p/n = (gas_gamma - 1)*(mass * E - 1/2 (rho u)^2)/rho.
 *
 * @param up Updater to run.
 * @param moms Fluid moments (rho, rho*ux, rho*uy, rho*uz, totalE).
 * @param out Output primitive moments.
 * @param out_coff Offset in out where to place primitive moments.
 */
void gkyl_gk_neut_fluid_prim_vars_udrift_temp_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff);

/**
 * Compute the LTE moments: density, drift velocity vector (ux, uy, uz)
 * and the temperature T = p/n = (gas_gamma - 1)*(mass * E - 1/2 (rho u)^2)/rho.
 *
 * @param up Updater to run.
 * @param moms Fluid moments (rho, rho*ux, rho*uy, rho*uz, totalE).
 * @param out Output primitive moments.
 * @param out_coff Offset in out where to place primitive moments.
 */
void gkyl_gk_neut_fluid_prim_vars_lte_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff);

/**
 * Free memory associated with this updater.
 * @param up Updater to free.
 */
void gkyl_gk_neut_fluid_prim_vars_release(gkyl_gk_neut_fluid_prim_vars *up);
