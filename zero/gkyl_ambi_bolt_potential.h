#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_ambi_bolt_potential gkyl_ambi_bolt_potential;

/**
 * Create new updater to compute the electrostatic potential assuming ambipolar
 * sheath particle fluxes and Boltzmann isothermal electrons.
 *
 * @param grid Cartesian grid dynamic field is defined on.
 * @param basis Basis object (configuration space).
 * @param local_range_ext Local extended range.
 * @param num_ghosts Number of ghosts in each dimension.
 * @param mass_e Electron mass.
 * @param charge_e Electron charge.
 * @param temp_e Electron temperature.
 * @param use_gpu Boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_ambi_bolt_potential* gkyl_ambi_bolt_potential_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_range *local_range_ext,
  const int *num_ghosts, double mass_e, double charge_e, double temp_e, bool use_gpu);

/**
 * Compute the ion density and electrostatic potential at the sheath entrance.
 * Below we use jac to refer to the fact that some of these quantities are
 * multiplied by the configuration space Jacobian (e.g. for field aligned
 * coordinates). In reality they are also multiplied by the phase space
 * Jacobian (approx B), but we focus on jac because need to divide by it to
 * get the potential (or multiply by 1/jac).
 *
 * @param Ambipolar, Boltzmann electron sheath potential updater. 
 * @param edge Lower (-1) or upper (1) boundary along the field line.
 * @param jacob_geo_inv Reciprocal of the configuration space Jacobian.
 * @param gamma_jac_i Ion particle flux at the sheath entrance times the
 *                    conf-space Jacobian.
 * @param m0_jac_i Ion number density times the conf-space Jacobian.
 * @param sheath_vals Ion number density and potential at the sheath entrance.
 */
void gkyl_ambi_bolt_potential_sheath_calc(struct gkyl_ambi_bolt_potential *up,
  enum gkyl_edge_loc edge, struct gkyl_array *jacob_geo_inv, struct gkyl_array *gamma_jac_i,
  struct gkyl_array *m0_jac_i, struct gkyl_array *sheath_vals);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_ambi_bolt_potential_release(gkyl_ambi_bolt_potential *up);
