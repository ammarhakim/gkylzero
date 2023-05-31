#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

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
void gkyl_calc_pkpm_vars_prim(struct gkyl_basis basis, const struct gkyl_range *range,
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
void gkyl_calc_pkpm_vars_source(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* qmem, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs);

/**
 * Compute parallel-kinetic-perpendicular-moment model distribution function
 * in the mirror force for the T_perp/m*G = T_perp/m*(F_0 - F_1) kinetic equation,
 * along with the vperp characteristics which are a pure source term in the first Laguerre moment update.
 * g_dist_source = [2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)), 
 *                  (-vpar div(b) + bb:grad(u) - div(u) - 2 nu) T_perp/m G + 2 nu vth^2 F_0 ]
 * First output is mirror force source *distribution*, second output is *total* vperp characteristics source.
 * Note that T_perp/m*G is the evolved quantity for the first Laguerre moment. 
 * Also outputs F_1 from T_perp/m*G for the evolution of F_2 if F_2 is present. 
 * To simplify internal Gkeyll logic, assume F_2 is present and output F_1 even if F_2 = 0.0.
 *
 * @param grid Grid (for getting cell spacing and cell center for vperp characteristics)
 * @param basis Basis functions used in expansions
 * @param conf_range Range to index configuration space fields
 * @param phase_range Range to calculate mirror force distribution function 
 * @param T_perp_over_m Input array of p_perp/rho = T_perp/m
 * @param T_perp_over_m_inv Input array of (T_perp/m)^-1
 * @param nu_vthsq Input array of nu*vth^2
 * @param pkpm_accel_vars Input arrary of pkpm acceleration variables ordered as (this updater needs div_b and p_perp_source):
 *        0: div_b (divergence of magnetic field unit vector)
          1: bb_grad_u (bb : grad(u))
          2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
          3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2 nu)
          4: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
 * @param fIn Input array of pkpm distribution functions: [F_0, T_perp/m G = T_perp/m (F_0 - F_1)]
 * @param F_k_p_1 Input array of k+1 distribution function. F_2 expansion is the first NP coefficients.
 * @param g_dist_source Output array: 2.0*(T_perp/m*G) + T_perp/m*(F_2 - F_0)
 * @param F_k_m_1 Output array of k-1 distribution function. F_1 expansion is the first NP coefficients.
 */
void gkyl_calc_pkpm_vars_dist_mirror_force(const struct gkyl_rect_grid *grid, struct gkyl_basis basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* T_perp_over_m, const struct gkyl_array* T_perp_over_m_inv, 
  const struct gkyl_array* nu_vthsq, const struct gkyl_array* pkpm_accel_vars, 
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
 * @param pkpm_div_ppar Input array of div(p_parallel b_hat) for computing pressure force
 * @param rho_inv Input array of 1/rho
 * @param T_perp_over_m Input array of p_perp/rho = T_perp/m
 * @param T_perp_over_m_inv Input array of (T_perp/m)^-1
 * @param nu Input array of collisionality
 * @param div_p Output array of divergence of pressure tensor
 * @param pkpm_accel_vars Output arrary of pkpm acceleration variables ordered as:
 *        0: div_b (divergence of magnetic field unit vector)
          1: bb_grad_u (bb : grad(u))
          2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
          3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2 nu)
          4: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
 */
void gkyl_calc_pkpm_vars_recovery(const struct gkyl_rect_grid *grid, 
  struct gkyl_basis basis, const struct gkyl_range *range, double nuHyp, 
  const struct gkyl_array* bvar, const struct gkyl_array* u_i, 
  const struct gkyl_array* p_ij, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* pkpm_div_ppar, const struct gkyl_array* rho_inv, const struct gkyl_array* T_perp_over_m, 
  const struct gkyl_array* T_perp_over_m_inv, const struct gkyl_array* nu, 
  struct gkyl_array* div_p, struct gkyl_array* pkpm_accel_vars);

/**
 * Compute divergence of parallel pressure for use in pressure force in kinetic equation.
 * Needed for self-consistency and to avoid spurious development of first moment of kinetic equation.
 *
 * @param grid Grid (for getting cell spacing and cell center)
 * @param basis Basis functions used in expansions
 * @param conf_range Range to index configuration space fields
 * @param phase_range Range to index phase space fields
 * @param bvar Input array of magnetic field unit vector and unit tensor
 * @param fIn Input array of pkpm distribution functions: [F_0, T_perp/m G = T_perp/m (F_0 - F_1)]
 * @param pkpm_div_ppar Output array of divergence of p_parallel b_hat
 */
void gkyl_calc_pkpm_vars_pressure(const struct gkyl_rect_grid *grid, struct gkyl_basis basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* bvar, const struct gkyl_array* f, 
  struct gkyl_array* pkpm_div_ppar);

/**
 * Host-side wrappers for prim vars operations on device
 */

void gkyl_calc_pkpm_vars_prim_cu(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* u_i, struct gkyl_array* p_ij, struct gkyl_array* T_ij, 
  struct gkyl_array* rho_inv, struct gkyl_array* T_perp_over_m, struct gkyl_array* T_perp_over_m_inv);

void gkyl_calc_pkpm_vars_source_cu(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* qmem, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs);

void gkyl_calc_pkpm_vars_dist_mirror_force_cu(const struct gkyl_rect_grid *grid, struct gkyl_basis basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* T_perp_over_m, const struct gkyl_array* T_perp_over_m_inv, 
  const struct gkyl_array* nu_vthsq, const struct gkyl_array* pkpm_accel_vars, 
  const struct gkyl_array* fIn, const struct gkyl_array* F_k_p_1,
  struct gkyl_array* g_dist_source, struct gkyl_array* F_k_m_1);

void gkyl_calc_pkpm_vars_recovery_cu(const struct gkyl_rect_grid *grid, 
  struct gkyl_basis basis, const struct gkyl_range *range, double nuHyp, 
  const struct gkyl_array* bvar, const struct gkyl_array* u_i, 
  const struct gkyl_array* p_ij, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* pkpm_div_ppar, const struct gkyl_array* rho_inv, const struct gkyl_array* T_perp_over_m, 
  const struct gkyl_array* T_perp_over_m_inv, const struct gkyl_array* nu, 
  struct gkyl_array* div_p, struct gkyl_array* pkpm_accel_vars);

void gkyl_calc_pkpm_vars_pressure_cu(const struct gkyl_rect_grid *grid, struct gkyl_basis basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* bvar, const struct gkyl_array* f, 
  struct gkyl_array* pkpm_div_ppar);

