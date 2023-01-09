#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
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
 * Compute primitive variables for parallel-kinetic-perpendicular-moment model.
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param bvar Input array of magnetic field unit vector and unit tensor
 * @param vlasov_pkpm_moms Input array of parallel-kinetic-perpendicular-moment kinetic moments [rho, p_parallel, p_perp]
 * @param euler_pkpm Input array of parallel-kinetic-perpendicular-moment fluid variables [rho ux, rho uy, rho uz]
 * @param u_i Output array of flow velocity 
 * @param p_ij Output array of pressure tensor
 * @param T_ij Output array of temperature tensor
 * @param rho_inv Output array of 1/rho
 * @param T_perp_over_m Output array of p_perp/rho = T_perp/m
 * @param T_perp_over_m_inv Output array of (T_perp/m)^-1
 */
void gkyl_calc_prim_vars_pkpm(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* u_i, struct gkyl_array* p_ij, struct gkyl_array* T_ij, 
  struct gkyl_array* rho_inv, struct gkyl_array* T_perp_over_m, struct gkyl_array* T_perp_over_m_inv);

/**
 * Compute parallel-kinetic-perpendicular-moment model source terms.
 *
 * @param basis Basis functions used in expansions
 * @param range Range to calculate source terms (configuration space range)
 * @param qmem Input array of q/m*EM fields
 * @param vlasov_pkpm_moms Input array of parallel-kinetic-perpendicular-moment kinetic moments [rho, p_parallel, p_perp]
 * @param euler_pkpm Input array of parallel-kinetic-perpendicular-moment fluid variables [rho ux, rho uy, rho uz]
 * @param rhs Output increment to fluid variables
 */
void gkyl_calc_prim_vars_pkpm_source(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* qmem, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs);

/**
 * Compute parallel-kinetic-perpendicular-moment model mirror force distribution function.
 * In the mirror force for the T_perp/m*G = T_perp/m*(F_0 - F_1) kinetic equation,
 * g_dist_source = 2.0*(T_perp/m*G) - T_perp/m*F_0 + T_perp/m*F_2
 * Note that T_perp/m*G is the evolved quantity for the first Laguerre moment. 
 * Also outputs F_1 from T_perp/m*G for the evolution of F_2 if F_2 is present. 
 * To simplify internal Gkeyll logic, assume F_2 is present and output F_1 even if F_2 = 0.0.
 *
 * @param basis Basis functions used in expansions
 * @param conf_range Range to index configuration space fields
 * @param phase_range Range to calculate mirror force distribution function 
 * @param T_perp_over_m Input array of p_perp/rho = T_perp/m
 * @param T_perp_over_m_inv Input array of (T_perp/m)^-1
 * @param fIn Input array of pkpm distribution functions: [F_0, T_perp/m G = T_perp/m (F_0 - F_1)]
 * @param F_k_p_1 Input array of k+1 distribution function. F_2 expansion is the first NP coefficients.
 * @param g_dist_source Output array: 2.0*(T_perp/m*G) + T_perp/m*(F_2 - F_0)
 * @param F_k_m_1 Output array of k-1 distribution function. F_1 expansion is the first NP coefficients.
 */
void gkyl_calc_prim_vars_pkpm_dist_mirror_force(struct gkyl_basis basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* T_perp_over_m, const struct gkyl_array* T_perp_over_m_inv, 
  const struct gkyl_array* fIn, const struct gkyl_array* F_k_p_1,
  struct gkyl_array* g_dist_source, struct gkyl_array* F_k_m_1);

/**
 * Compute needed gradient quantities with recovery for discretization of the 
 * parallel-kinetic-perpendicular-moment (pkpm) model. These include div(p) and the forces and sources in the pkpm kinetic equation
 *
 * @param grid Grid (for getting cell spacing)
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param bvar Input array of magnetic field unit vector and unit tensor
 * @param u_i Input array of flow velocity 
 * @param p_ij Input array of pressure tensor
 * @param vlasov_pkpm_moms Input array of parallel-kinetic-perpendicular-moment kinetic moments [rho, p_parallel, p_perp]
 * @param euler_pkpm Input array parallel-kinetic-perpendicular-moment fluid variables [rho ux, rho uy, rho uz]
 * @param rho_inv Input array of 1/rho
 * @param T_perp_over_m Input array of p_perp/rho = T_perp/m
 * @param T_perp_over_m_inv Input array of (T_perp/m)^-1
 * @param nu Input array of collisionality
 * @param nu_vthsq Input array of nu*vth^2
 * @param div_p Output array of divergence of pressure tensor
 * @param pkpm_accel_vars Output arrary of pkpm acceleration variables ordered as:
 *        0: div_b (divergence of magnetic field unit vector)
          1: bb_grad_u (bb : grad(u))
          2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
          3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - nu + nu rho vth^2/p_perp)
          4: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
 */
void gkyl_calc_prim_vars_pkpm_recovery(const struct gkyl_rect_grid *grid, 
  struct gkyl_basis basis, const struct gkyl_range *range, double nuHyp, 
  const struct gkyl_array* bvar, const struct gkyl_array* u_i, 
  const struct gkyl_array* p_ij, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* rho_inv, const struct gkyl_array* T_perp_over_m, const struct gkyl_array* T_perp_over_m_inv, 
  const struct gkyl_array* nu, const struct gkyl_array* nu_vthsq, 
  struct gkyl_array* div_p, struct gkyl_array* pkpm_accel_vars);

/**
 * Host-side wrappers for prim vars operations on device
 */

void gkyl_calc_prim_vars_pkpm_cu(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* u_i, struct gkyl_array* p_ij, struct gkyl_array* T_ij, 
  struct gkyl_array* rho_inv, struct gkyl_array* T_perp_over_m, struct gkyl_array* T_perp_over_m_inv);

void gkyl_calc_prim_vars_pkpm_source_cu(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* qmem, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs);

void gkyl_calc_prim_vars_pkpm_dist_mirror_force_cu(struct gkyl_basis basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* T_perp_over_m, const struct gkyl_array* T_perp_over_m_inv, 
  const struct gkyl_array* fIn, const struct gkyl_array* F_k_p_1,
  struct gkyl_array* g_dist_source, struct gkyl_array* F_k_m_1);

void gkyl_calc_prim_vars_pkpm_recovery_cu(const struct gkyl_rect_grid *grid, 
  struct gkyl_basis basis, const struct gkyl_range *range, double nuHyp, 
  const struct gkyl_array* bvar, const struct gkyl_array* u_i, 
  const struct gkyl_array* p_ij, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* rho_inv, const struct gkyl_array* T_perp_over_m, const struct gkyl_array* T_perp_over_m_inv, 
  const struct gkyl_array* nu, const struct gkyl_array* nu_vthsq, 
  struct gkyl_array* div_p, struct gkyl_array* pkpm_accel_vars);
