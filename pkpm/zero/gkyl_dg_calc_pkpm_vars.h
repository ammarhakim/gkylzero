#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_calc_pkpm_vars gkyl_dg_calc_pkpm_vars;

/**
 * Create new updater for the configuration-space variables of the energy-constrained PKPM
 * solver (p=1 serendipity): the weak division M/rho for the flow, the nodal primitive
 * variables, the closures with the velocity-space edges, the weak gradients, the explicit
 * source, the initial-condition normalization and the diagnostics.
 *
 * @param conf_grid Configuration space grid (for getting cell spacing and cell center)
 * @param cbasis    Configuration space basis functions
 * @param mem_range Configuration space range of the weak divisions (the cells updated)
 * @param use_gpu   Whether to run on host or device
 * @return New updater pointer.
 */
struct gkyl_dg_calc_pkpm_vars *gkyl_dg_calc_pkpm_vars_new(
  const struct gkyl_rect_grid *conf_grid, const struct gkyl_basis *cbasis,
  const struct gkyl_range *mem_range, bool use_gpu
);

/**
 * Create new updater on NV-GPU.
 */
struct gkyl_dg_calc_pkpm_vars *gkyl_dg_calc_pkpm_vars_cu_dev_new(
  const struct gkyl_rect_grid *conf_grid, const struct gkyl_basis *cbasis,
  const struct gkyl_range *mem_range
);

/**
 * Compute volume expansion of flow velocity u in the PKPM system.
 *
 * @param up Updater for computing pkpm variables 
 * @param vlasov_pkpm_moms Input array of pkpm kinetic moments [rho, p_parallel, p_perp]
 * @param euler_pkpm Input array of pkpm fluid variables [rho ux, rho uy, rho uz]
 * @param cell_avg_prim Array for storing boolean value of whether rho, p_parallel, p_perp are negative at corners 
 *                      Note: Only used for diagnostic purposes (not for adjusting solution)
 * @param pkpm_u        Output array of volume expansion of flow velocity [ux, uy, uz]
 */
void gkyl_dg_calc_pkpm_vars_u(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *euler_pkpm, struct gkyl_array *cell_avg_prim, struct gkyl_array *pkpm_u
);

/**
 * Compute the weak (surface-averaged) gradient of the flow velocity, the p=1 expansions
 * of d_i u_j for every configuration direction i and component j, stored per cell as
 * [direction][component][mode] (cdim * 3 * num_basis entries). The gradient uses the
 * average of the two traces at each face, so it needs u in the neighbors of every cell
 * of conf_range (the ghost cells). It is the gradient the kinetic sources bb:grad(u)
 * and div(u) are formed from. p=1 serendipity only.
 *
 * @param up         Updater for computing pkpm variables
 * @param conf_range Configuration space range (interior cells)
 * @param u          Input flow velocity [ux, uy, uz], defined on the extended range
 * @param grad_u     Output weak gradient of u
 */
void gkyl_dg_calc_pkpm_vars_grad_u(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *u, struct gkyl_array *grad_u
);

/**
 * Closures c_0, alpha_0 = a_P/xi and c of the normalized kinetic update, with the
 * characteristic velocity-space edges: a 3x3 solve per Gauss-Legendre node for the three
 * constraints M0(F) = 1, M1(F) = 0, M2(F) + 2 M0(G) = 3 of the complete right-hand side
 * (configuration pass + closure terms + the upwind edge fluxes against an empty state
 * outside the velocity domain), see pkpm-vars-closure.mac. The edge fluxes are written in
 * the stored face format: the lower one into the lowest velocity cell's slot of
 * vel_flux_surf, the upper one into the slot of the velocity ghost cell above the
 * uppermost cell (the array layout must have that ghost cell), and their time-step rates
 * are added to cflrate in the two edge cells. The upwind switch of
 * the edges uses the closures of the previous stage (the output arrays on input). At a
 * node with M0(F) <= 0 or a non-positive determinant all closures and edge fluxes of the
 * cell are zero; the number of such nodes is returned. p=1 serendipity only.
 *
 * @param up            Updater for computing pkpm variables
 * @param conf_range    Configuration space range
 * @param phase_range   Phase-space range (its velocity extent locates the edge cells)
 * @param phase_grid    Phase-space grid (velocity extent and cell size)
 * @param norm_moms     Raw moments [M0, M2, M0_G, M1, M3, M1_G] of the normalized distributions
 * @param rhs_moms      The same moments of the configuration-pass right-hand side
 * @param bvar          Magnetic field unit vector
 * @param grad_b        Weak gradient of b
 * @param grad_u        Weak gradient of u
 * @param xi            Thermal speed [xi, xi^2, p]
 * @param grad_xi       Weak gradient of [xi, 0, 0]
 * @param theta         T_perp/T
 * @param fin           Normalized distributions [F, G]
 * @param c_0           Output normalization closure (input: previous stage)
 * @param alpha_0       Output pressure force a_P/xi (input: previous stage)
 * @param c_h           Output heating rate (input: previous stage)
 * @param vel_flux_surf Stored velocity face fluxes: the slots of the lowest cells (lower
 *                      edge) and of the velocity ghost cells above the uppermost cells
 *                      (upper edge) are written; phase_range must be laid out with that
 *                      ghost cell
 * @param cflrate       Phase-space time-step rates, incremented in the edge cells
 * @return Number of nodes where the closures are undefined
 */
int gkyl_dg_calc_pkpm_vars_closure(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_array *norm_moms, const struct gkyl_array *rhs_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *grad_u,
  const struct gkyl_array *xi, const struct gkyl_array *grad_xi, const struct gkyl_array *theta,
  const struct gkyl_array *fin, struct gkyl_array *c_0, struct gkyl_array *alpha_0, struct gkyl_array *c_h,
  struct gkyl_array *vel_flux_surf, struct gkyl_array *cflrate
);

/**
 * Compute the nodal primitive variables of the energy-constrained fluid state
 * [rho, rho ux, rho uy, rho uz, E]: the flow u = M/rho and [xi, xi^2, p] with
 * xi^2 = (2 E/rho - |u|^2)/3, p = rho xi^2, by the nodal quotient at the 2^cdim
 * Gauss-Legendre nodes interpolated to p=1. Nodes with rho <= 0 or xi^2 <= 0 get zero
 * and are counted. p=1 serendipity only.
 *
 * @param up         Updater for computing pkpm variables
 * @param conf_range Configuration space range
 * @param fluid      Input fluid state
 * @param u          Output flow velocity (3 expansions per cell)
 * @param xi         Output [xi, xi^2, p] (3 expansions per cell)
 * @return Number of inadmissible nodes in this call
 */
int gkyl_dg_calc_pkpm_vars_prim(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, struct gkyl_array *u, struct gkyl_array *xi
);

/**
 * Cell integrals of the diagnostics of the [rho, M, E] state: [rho, rho ux, rho uy,
 * rho uz, rho ux^2, rho uy^2, rho uz^2, p_par, p_perp] (accumulated). p=1 serendipity only.
 *
 * @param up               Updater for computing pkpm variables
 * @param conf_range       Configuration space range
 * @param fluid            Input fluid state [rho, rho ux, rho uy, rho uz, E]
 * @param vlasov_pkpm_moms Input kinetic coupling moments [p_par, p_perp, q_par, q_perp]
 * @param u                Input flow velocity
 * @param int_pkpm_vars    Output cell integrals (9 per cell)
 */
void gkyl_dg_calc_pkpm_vars_int(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *u, struct gkyl_array *int_pkpm_vars
);

/**
 * The three constraints of the normalized distributions, [M0(F), M1(F), (M2(F) + 2 M0(G))/3],
 * from the raw normalized moments [M0, M2, M0_G, M1, M3, M1_G]: they are [1, 0, 1] in every
 * mode of every cell when the closures hold.
 *
 * @param up          Updater (unused; kept for symmetry)
 * @param conf_range  Configuration space range
 * @param num_basis   Number of configuration basis functions
 * @param norm_moms   Input normalized moments (6 expansions per cell)
 * @param constraints Output constraints (3 expansions per cell)
 */
void gkyl_dg_calc_pkpm_vars_constraints(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, int num_basis,
  const struct gkyl_array *norm_moms, struct gkyl_array *constraints
);

/**
 * Compute grad(ln n) = grad(rho)/rho, the configuration-space factor of the source
 * -v_par b.grad(ln n) of the factored kinetic equations, by the nodal quotient of the
 * weak gradient of rho (computed with gkyl_dg_calc_pkpm_vars_grad_u on the array
 * [rho, 0, 0]) at the 2^cdim Gauss-Legendre nodes. At a node where rho <= 0 the quotient is
 * set to zero; the number of such nodes is returned. p=1 serendipity only.
 *
 * @param up         Updater for computing pkpm variables
 * @param conf_range Configuration space range
 * @param rho        Input mass density
 * @param grad_rho   Input weak gradient of [rho, 0, 0] (cdim x 3 expansions per cell)
 * @param grad_ln_n  Output grad(ln n) (cdim expansions per cell)
 * @return Number of nodes where rho <= 0 in this call
 */
int gkyl_dg_calc_pkpm_vars_grad_ln_n(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *rho, const struct gkyl_array *grad_rho, struct gkyl_array *grad_ln_n
);

/**
 * Divide the distributions [F_0, G] weakly by the number density n to form the
 * normalized state [F_0/n, G/n] of the factored update (initial conditions). The
 * division is exact for a p=1 n positive at the 2^cdim Gauss-Legendre nodes; at a node
 * where n <= 0 the quotient is set to zero and counted. p=1 serendipity only.
 *
 * @param up          Updater for computing pkpm variables
 * @param conf_range  Configuration space range for indexing n
 * @param phase_range Phase space range to operate on
 * @param n           Input number density
 * @param f           Input distributions [F_0, G]
 * @param fbar        Output normalized distributions [F_0/n, G/n]
 * @return Number of nodes where n <= 0 in this call
 */
int gkyl_dg_calc_pkpm_vars_normalize_dist(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *n, const struct gkyl_array *f,
  struct gkyl_array *fbar
);

/**
 * Construct the diagnostic arrays of the refactored solver:
 * fluid_io = [rho, M_i, P_ij + M_i u_j] (10 components) and
 * pkpm_vars_io = [u_x, u_y, u_z, theta, a_P, div b, bb:grad u, theta div b] (8 components),
 * with the products projected onto p=1. p=1 serendipity only.
 *
 * @param up               Updater for computing pkpm variables 
 * @param conf_range       Configuration space range
 * @param vlasov_pkpm_moms Input array of pkpm kinetic moments [rho, p_parallel, p_perp, M1]
 * @param euler_pkpm       Input momentum [rho ux, rho uy, rho uz]
 * @param bvar             Input magnetic field unit vector [bx, by, bz]
 * @param grad_b           Input weak gradient of b (gkyl_dg_calc_pkpm_bvar_grad)
 * @param u                Input flow velocity [ux, uy, uz]
 * @param grad_u           Input weak gradient of u (gkyl_dg_calc_pkpm_vars_grad_u)
 * @param theta            Input T_perp/m
 * @param a_P              Input pressure force
 * @param fluid_io         Output conserved fluid variables (10 components)
 * @param pkpm_vars_io     Output pkpm variables (8 components)
 */
void gkyl_dg_calc_pkpm_vars_diag(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *euler_pkpm,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *theta, const struct gkyl_array *alpha_0,
  const struct gkyl_array *xi, const struct gkyl_array *c_h, struct gkyl_array *fluid_io,
  struct gkyl_array *pkpm_vars_io
);

/**
 * Compute integrated PKPM variables (rho, rhoux, rhouy, rhouz, rhoux^2, rhouy^2, rhouz^2, p_parallel, p_perp).
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param vlasov_pkpm_moms Input array of pkpm kinetic moments [rho, p_parallel, p_perp]
 * @param euler_pkpm Input array of pkpm fluid variables [rho ux, rho uy, rho uz]
 * @param prim Input array of primitive moments [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp]
 * @param int_pkpm_vars Output array of integrated variables (6 components)
 */
void gkyl_dg_calc_pkpm_integrated_vars(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *euler_pkpm,
  const struct gkyl_array *prim, struct gkyl_array *pkpm_int_vars
);

/**
 * Compute pkpm model source terms.
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param qmem Input array of q/m*EM fields
 * @param vlasov_pkpm_moms Input array of pkpm kinetic moments [rho, p_parallel, p_perp]
 * @param euler_pkpm Input array of pkpm fluid variables [rho ux, rho uy, rho uz]
 * @param rhs Output increment to fluid variables
 */
void gkyl_dg_calc_pkpm_vars_source(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *qmem, const struct gkyl_array *fluid, struct gkyl_array *rhs
);

/**
 * Delete pointer to updater to compute pkpm variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_pkpm_vars_release(struct gkyl_dg_calc_pkpm_vars *up);

/**
 * Host-side wrappers for pkpm vars operations on device
 */

void gkyl_dg_calc_pkpm_vars_u_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *euler_pkpm, struct gkyl_array *cell_avg_prim, struct gkyl_array *pkpm_u
);

int gkyl_dg_calc_pkpm_vars_closure_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_array *norm_moms, const struct gkyl_array *rhs_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *grad_u,
  const struct gkyl_array *xi, const struct gkyl_array *grad_xi, const struct gkyl_array *theta,
  const struct gkyl_array *fin, struct gkyl_array *c_0, struct gkyl_array *alpha_0, struct gkyl_array *c_h,
  struct gkyl_array *vel_flux_surf, struct gkyl_array *cflrate
);
int gkyl_dg_calc_pkpm_vars_prim_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, struct gkyl_array *u, struct gkyl_array *xi
);
void gkyl_dg_calc_pkpm_vars_int_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *u, struct gkyl_array *int_pkpm_vars
);
int gkyl_dg_calc_pkpm_vars_grad_ln_n_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *rho, const struct gkyl_array *grad_rho, struct gkyl_array *grad_ln_n
);
int gkyl_dg_calc_pkpm_vars_normalize_dist_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *n, const struct gkyl_array *f,
  struct gkyl_array *fbar
);
void gkyl_dg_calc_pkpm_vars_diag_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *euler_pkpm,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *theta, const struct gkyl_array *alpha_0,
  const struct gkyl_array *xi, const struct gkyl_array *c_h, struct gkyl_array *fluid_io,
  struct gkyl_array *pkpm_vars_io
);

void gkyl_dg_calc_pkpm_vars_grad_u_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *u, struct gkyl_array *grad_u
);

void gkyl_dg_calc_pkpm_integrated_vars_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *euler_pkpm,
  const struct gkyl_array *prim, struct gkyl_array *pkpm_int_vars
);

void gkyl_dg_calc_pkpm_vars_source_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *qmem, const struct gkyl_array *fluid, struct gkyl_array *rhs
);

