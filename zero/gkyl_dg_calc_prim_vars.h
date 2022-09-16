#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>

/**
 * Compute u from state vector which contains *both* density and momentum.
 * Data always organized such that density is 0th component of statevec
 * & momentum is 1st, 2nd, and 3rd component of statevec.
 * Examples include: 
 * Isothermal Euler state vector [rho, rhoux, rhouy, rhouz],
 * Euler state vector [rho, rhoux, rhouy, rhouz, E],
 * Nonrelativistic Vlasov Five Moments array [M0, M1i, M2] (u_i = M1i/M0),
 * Special Relativistic Vlasov N_i array [Gamma*n, Gamma*n*u] where u is the bulk flow
 *
 * @param mem Pre-allocated space for use in the division 
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param statevec Input state vector which contains *both* density and momentum
 * @param u_i Output array of bulk flow velocity
 */
void gkyl_calc_prim_vars_u_from_statevec(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* statevec, struct gkyl_array* u_i);

/**
 * Compute u from input density and momentum vectors.
 * Examples include: 
 * Nonrelativistic Vlasov M0 & M1i (u_i = M1i/M0),
 * Special Relativistic M0 (Gamma*n) and M1i (Gamma*n*u) where u is the bulk flow
 * Parallel-kinetic-perpendicular-moment model with vlasov_pkpm_moms and euler_pkpm
 *
 * @param mem Pre-allocated space for use in the division 
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param rho Input density 
 * @param rhou Input momentum
 * @param u_i Output array of bulk flow velocity
 */
void gkyl_calc_prim_vars_u_from_rhou(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* rho, const struct gkyl_array* rhou, struct gkyl_array* u_i);

/**
 * Compute pressure from state vector which contains density, momentum, and energy.
 * Data always organized such that density is 0th component of statevec, momentum is 
 * 1st, 2nd, and 3rd component of statevec, and energy is 4th component of statevec.
 * Note: Pressure computation requires computation of u, bulk flow velocity.
 * Examples include: 
 * Euler state vector [rho, rhoux, rhouy, rhouz, E],
 * Nonrelativistic Vlasov Five Moments array [M0, M1i, M2]
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param p_fac Factor for obtaining pressure (gas_gamma for Euler)
 * @param u_i Input array of bulk flow velocity
 * @param statevec Input state vector which contains density, momentum, and energy
 * @param p_ij Output array of pressure (scalar pressure)
 */
void gkyl_calc_prim_vars_p_from_statevec(struct gkyl_basis basis, const struct gkyl_range *range,
  const double p_fac, const struct gkyl_array* u_i, const struct gkyl_array* statevec, 
  struct gkyl_array* p_ij);

/**
 * Compute pressure from parallel-kinetic-perpendicular-moment inputs.
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param u_i Input array of bulk flow velocity
 * @param bvar Input array of magnetic field unit vector and unit tensor
 * @param vlasov_pkpm_moms Input array parallel-kinetic-perpendicular-moment kinetic moments
 * @param euler_pkpm Input array parallel-kinetic-perpendicular-moment fluid variables
 * @param p_ij Output array of pressure tensor
 */
void gkyl_calc_prim_vars_p_pkpm(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* u_i, const struct gkyl_array* bvar, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* p_ij);
