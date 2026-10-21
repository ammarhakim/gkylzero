#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_pkpm_harm_exp gkyl_dg_pkpm_harm_exp;

// Input packaged as a struct
struct gkyl_dg_pkpm_harm_exp_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid (the velocity cell width).
  const struct gkyl_basis *conf_basis; // Configuration-space basis (p=1 serendipity).
  double qbym; // Charge to mass ratio of the species.
  double mix; // Gyrophase-mixing rate in units of the local |Omega| (0: no mixing).
  bool use_gpu; // Whether to run on the GPU.
};

/**
 * Create a new updater for the exponential Runge-Kutta update of the first-harmonic blocks
 * M0_perp and N_perp of the normalized PKPM distributions [F, G, M0_perp, N_perp] (plan 17.6,
 * 17.7). The linear part
 *
 *   L_n U = (q/m) U x B_n
 *
 * uses the magnetic field B_n at the start of the step (the gyration at Omega_n about b_n). It is
 * applied node by node at the Gauss-Legendre nodes of the phase cell (2 per configuration direction
 * times 3 in w), where the exponential and phi functions of the rotation are closed form. With z = h L_n
 * and the stage right-hand sides N_i = rhs_i - L_n U_i, where rhs_i contains the gyration in the
 * field of the stage, so N_i keeps its change across the step:
 *
 *   U2 = e^z M + h phi_1(z) N_1,
 *   U3 = e^(z/2) M + (h/4) phi_1(z/2) (N_1 + N_2),
 *   M^(n+1) = e^z M + h [(phi_1 - 3 phi_2 + 4 phi_3) N_1 + (-phi_2 + 4 phi_3) N_2
 *                        + (4 phi_2 - 8 phi_3) N_3].
 *
 * The component along b_n is not rotated (its weights are those of SSP-RK3), and with B_n = 0 the
 * update is SSP-RK3. F and G are not touched. p=1 serendipity only.
 *
 * Gyrophase mixing (plan 17.24): agyrotropy beyond what one harmonic can represent relaxes to that
 * bound, d(M0_perp, N_perp)/dt = -nu ((M0_perp, N_perp) - P(M0_perp, N_perp)) with nu = mix |Omega|
 * and P the projection onto r <= 1 at each phase node, r = max(sqrt2 |M0_perp|/sqrt(F G),
 * |N_perp| sqrt(F)/G^(3/2)). As a rate on the first harmonic it is nu (1 - 1/r)_+, frozen at the
 * start of the step (capped at ln(r)/h, which takes the unforced state to the bound in one step) and
 * part of L on the plane perpendicular to b_n: z = h(-nu (1 - 1/r)_+ + i |Omega|) at that node.
 * It changes int M0_perp dw, which gkyl_dg_pkpm_harm_exp_restore restores.
 *
 * @param inp Input parameters.
 * @return New updater pointer.
 */
struct gkyl_dg_pkpm_harm_exp *gkyl_dg_pkpm_harm_exp_new(const struct gkyl_dg_pkpm_harm_exp_inp *inp);

/**
 * Create the updater on NV-GPU. See the new method for documentation.
 */
struct gkyl_dg_pkpm_harm_exp *gkyl_dg_pkpm_harm_exp_cu_dev_new(const struct gkyl_dg_pkpm_harm_exp_inp *inp);

/**
 * Number of components of the configuration-space coefficient array (the frozen linear part at
 * the configuration nodes).
 *
 * @param up Updater.
 * @return Number of components.
 */
int gkyl_dg_pkpm_harm_exp_num_coef(const struct gkyl_dg_pkpm_harm_exp *up);

/**
 * Number of components of the phase-space stage arrays S and A (nodal, six vector components per
 * phase node).
 *
 * @param up Updater.
 * @return Number of components.
 */
int gkyl_dg_pkpm_harm_exp_num_stage(const struct gkyl_dg_pkpm_harm_exp *up);

/**
 * Freeze the linear part at the start of the step: (q/m) B_n, its direction and the exponential and
 * phi functions of z = h L_n at the configuration nodes.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param h Time step.
 * @param em Total physical electromagnetic field at the start of the step.
 * @param coef Coefficients at the configuration nodes (output).
 */
void gkyl_dg_pkpm_harm_exp_coef(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range, double h,
  const struct gkyl_array *em, struct gkyl_array *coef
);

/**
 * Accumulate the right-hand side of stage s (1, 2 or 3) of the first-harmonic blocks:
 * N_s = rhs_s - L_n U_s, then S = N_1 (s = 1) or S += N_2 (s = 2), and A = h w_1(z) N_1 (s = 1)
 * or A += h w_s(z) N_s (s = 2, 3), with w_s the weights of the final update. At stage 1 (U_1 the
 * state at the start of the step) the mixing rate of every phase node is set first, and act is set to
 * 1 in the configuration cells with a damped node (act is not cleared here).
 *
 * @param up Updater.
 * @param stage Stage (1, 2 or 3).
 * @param h Time step (the one the coefficients were frozen with).
 * @param conf_range Configuration-space range.
 * @param phase_range Phase-space range.
 * @param coef Frozen coefficients.
 * @param fin Stage state U_s (eight blocks).
 * @param rhs Stage right-hand side (eight blocks; the first-harmonic blocks are read).
 * @param nu Mixing rate at the phase nodes (written at stage 1; unused without mixing).
 * @param act Configuration-space flag of damped columns (set at stage 1; unused without mixing).
 * @param S Nodal N_1 (+ N_2), updated.
 * @param A Nodal weighted sum of the final update, updated.
 */
void gkyl_dg_pkpm_harm_exp_accum(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *fin,
  const struct gkyl_array *rhs, struct gkyl_array *nu, struct gkyl_array *act, struct gkyl_array *S,
  struct gkyl_array *A
);

/**
 * Write the first-harmonic blocks of the output of stage s in the Shu-Osher form of the SSP-RK3
 * driver: the value X_s whose SSP-RK3 combination with the state M at the start of the step gives
 * the exponential stage, X_1 = U2, X_2 = 4 U3 - 3 M (3/4 M + 1/4 X_2 = U3), X_3 = (3 M^(n+1) - M)/2
 * (1/3 M + 2/3 X_3 = M^(n+1)). F and G of out are not touched.
 *
 * @param up Updater.
 * @param stage Stage (1, 2 or 3).
 * @param h Time step (the one the coefficients were frozen with).
 * @param conf_range Configuration-space range.
 * @param phase_range Phase-space range.
 * @param coef Frozen coefficients.
 * @param nu Mixing rate at the phase nodes (unused without mixing).
 * @param f0 State M at the start of the step (eight blocks).
 * @param S Nodal N_1 + N_2 (stage 2) or N_1 (stage 1).
 * @param A Nodal weighted sum of the final update (stage 3).
 * @param out Stage output (eight blocks; blocks 2-7 written).
 */
void gkyl_dg_pkpm_harm_exp_state(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *nu,
  const struct gkyl_array *f0, const struct gkyl_array *S, const struct gkyl_array *A,
  struct gkyl_array *out
);

/**
 * Number of components of the configuration-space array of the w integrals of the restoration.
 *
 * @param up Updater.
 * @return Number of components.
 */
int gkyl_dg_pkpm_harm_exp_num_col(const struct gkyl_dg_pkpm_harm_exp *up);

/**
 * Restore int M0_perp dw = 0 in the configuration cells flagged by act (those with a damped node):
 * with m = int M0_perp dw / int F dw at the configuration nodes, M0_perp -= F m and
 * N_perp -= 2 G m at every w. (F, 2G) dv is the first-harmonic change of a perpendicular frame shift
 * dv, the direction of the perpendicular closure force, so the momentum the mixing removes is
 * returned as a frame shift; F and G (and with them the other constraints and the energy) are not
 * changed. The products are exact at the phase nodes, so int M0_perp dw = 0 holds to round-off. Applied
 * to a stage output it also holds for the SSP-RK3 combinations of the driver (they are linear).
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param phase_range Phase-space range (its velocity extent is the column).
 * @param act Configuration-space flag of damped columns.
 * @param col Configuration-space work array (gkyl_dg_pkpm_harm_exp_num_col components).
 * @param out Stage output (eight blocks), restored in place.
 */
void gkyl_dg_pkpm_harm_exp_restore(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *act, struct gkyl_array *col,
  struct gkyl_array *out
);

/**
 * Host-side wrappers of the device kernels (see the host methods for documentation).
 */
void gkyl_dg_pkpm_harm_exp_coef_cu(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range, double h,
  const struct gkyl_array *em, struct gkyl_array *coef
);
void gkyl_dg_pkpm_harm_exp_accum_cu(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *fin,
  const struct gkyl_array *rhs, struct gkyl_array *nu, struct gkyl_array *act, struct gkyl_array *S,
  struct gkyl_array *A
);
void gkyl_dg_pkpm_harm_exp_state_cu(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *nu,
  const struct gkyl_array *f0, const struct gkyl_array *S, const struct gkyl_array *A,
  struct gkyl_array *out
);
void gkyl_dg_pkpm_harm_exp_restore_cu(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *act, struct gkyl_array *col,
  struct gkyl_array *out
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_pkpm_harm_exp_release(struct gkyl_dg_pkpm_harm_exp *up);
