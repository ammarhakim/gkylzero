#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Object type
typedef struct gkyl_dg_calc_pkpm_vars gkyl_dg_calc_pkpm_vars;

/**
 * Create new updater to compute pkpm variables needed in 
 * updates and used for diagnostics. Methods compute:
 * pkpm_u : flow velocity, order p, [ux, uy, uz] 
 * pkpm_u_surf : surface expansion of flow velocity, order p
 *               [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
 *                ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 
 *                ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr] 
 * 
 * p_ij : pressure tensor, order 2*p, (p_par - p_perp) b_i b_j + p_perp g_ij
 * prim : primitive variables, order 2*p, [1/rho*div(p_par b), T_perp/m, m/T_perp, 3*Txx/m, 3*Tyy/m, 3*Tzz/m]
 * 
 * pkpm_accel_vars : 0 : p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
                     1 : bb_grad_u (bb : grad(u))
                     2 : p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
                     3 : p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2 nu)
                     all order 2*p

 * pkpm_int_vars : integrated PKPM variables, order 2*p, (rho, rhoux, rhouy, rhouz, rhoux^2, rhouy^2, rhouz^2, p_parallel, p_perp)
 * pkpm fluid source : Explicit source terms in momentum solve, order p, q/m rho (E_i + epsilon_ijk u_j B_k)
 * pkpm limiter : Limit slopes of fluid variables in PKPM system, order p (since momentum equation is order p)
 * 
 * Note : Quantities derived from, or coupled to, the kinetic equation are order 2*p
 *        Quantities derived from the momentum equation are order p
 * 
 * @param conf_grid   Configuration space grid (for getting cell spacing and cell center)
 * @param cbasis      Configuration space basis functions for momentum (order p basis)
 * @param cbasis_2p   Configuration space basis functions for kinetic variables (order 2*p basis)
 * @param mem_range   Configuration space range that sets the size of the bin_op memory
 *                    for computing primitive moments. Note range is stored so 
 *                    updater loops over consistent range for primitive moments
 * @param wv_eqn      Wave equation (stores function pointers for computing waves and limiting solution)
 * @param geom        Wave geometry object for computing waves in local coordinate system
 * @param limiter_fac Optional parameter for changing diffusion in sloper limiter 
 *                    by changing relationship between slopes and cell average differences.
 *                    By default, this factor is 1/sqrt(3) because cell_avg(f) = f0/sqrt(2^cdim)
 *                    and a cell slope estimate from two adjacent cells is (for the x variation): 
 *                    integral(psi_1 [cell_avg(f_{i+1}) - cell_avg(f_{i})]*x) = sqrt(2^cdim)/sqrt(3)*[cell_avg(f_{i+1}) - cell_avg(f_{i})]
 *                    where psi_1 is the x cell slope basis in our orthonormal expansion psi_1 = sqrt(3)/sqrt(2^cdim)*x
 *                    This factor can be made smaller (larger) to increase (decrease) the diffusion from the slope limiter
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_pkpm_vars* 
gkyl_dg_calc_pkpm_vars_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_basis* cbasis_2p, 
  const struct gkyl_range *mem_range, 
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *geom, 
  double limiter_fac, bool use_gpu);

/**
 * Create new updater to compute pkpm variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_pkpm_vars* 
gkyl_dg_calc_pkpm_vars_cu_dev_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_basis* cbasis_2p, 
  const struct gkyl_range *mem_range, 
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *geom, 
  double limiter_fac);

/**
 * Compute all of the pkpm primitive moments. All quantities order 2*p
 *
 * @param up Updater for computing pkpm variables 
 * @param vlasov_pkpm_moms Input array of pkpm kinetic moments [rho, p_parallel, p_perp]
 * @param p_ij Input pressure tensor p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij
 * @param pkpm_div_ppar Input array of div(p_parallel b_hat) for computing first component of pressure force
 * @param div_b Input array of div(b) for computing second component of pressure force (p_perp div(b)/rho)
 * @param cell_avg_prim Array for storing boolean value of whether [rho, p = p_parallel + 2*p_perp] 
 *                      are negative at control points. 
 *                      Note: Only used for diagnostic purposes (not for adjusting solution)
 * @param prim Output array of volume expansion of primitive variables
 *             [1/rho*div(p_par b), T_perp/m, m/T_perp, 3*Txx/m, 3*Tyy/m, 3*Tzz/m, 1/rho*p_perp div(b)]
 */
void gkyl_dg_calc_pkpm_vars_advance(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* p_ij, 
  const struct gkyl_array* pkpm_div_ppar, const struct gkyl_array* div_b, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* prim, struct gkyl_array* pkpm_accel);

/**
 * Compute the needed surface expansions of the pkpm primitive moments. All quantities order 2*p.
 * 2*cdim components:  (2 components) in each dimension (cdim components)
 * 3*Txx/m at the left and right x surfaces 
 * 3*Tyy/m at the left and right y surfaces 
 * 3*Tzz/m at the left and right z surfaces
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param vlasov_pkpm_moms Input array of pkpm kinetic moments [rho, p_parallel, p_perp]
 * @param p_ij Input pressure tensor p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij
 * @param prim_surf Output array of surface expansion of primitive moments 
 *                  [3.0*Txx_xl/m, 3.0*Txx_xr/m, 3.0*Tyy_yl/m, 3.0*Tyy_yr/m, 3.0*Tzz_zl/m, 3.0*Tzz_zr/m] 
 */
void gkyl_dg_calc_pkpm_vars_surf_advance(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* p_ij, 
  struct gkyl_array* prim_surf);

/**
 * Compute volume expansion of flow velocity u in the PKPM system. Flow velocity is order p.
 *
 * @param up Updater for computing pkpm variables 
 * @param vlasov_pkpm_moms Input array of pkpm kinetic moments [rho, p_parallel, p_perp]
 * @param euler_pkpm Input array of pkpm fluid variables [rho ux, rho uy, rho uz]
 * @param cell_avg_prim Array for storing boolean value of whether [rho, p = p_parallel + 2*p_perp] 
 *                      are negative at control points. 
 *                      Note: Only used for diagnostic purposes (not for adjusting solution)
 * @param pkpm_u        Output array of volume expansion of flow velocity [ux, uy, uz]
 */
void gkyl_dg_calc_pkpm_vars_u(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* pkpm_u);

/**
 * Compute surface expansion of flow velocity u in the PKPM system. Flow velocity is order p.
 *
 * @param up Updater for computing pkpm variables 
 * @param pkpm_u Input array of volume expansion of flow velocity [ux, uy, uz]
 * @param pkpm_u_surf Output array of surface expansion of flow velocity 
 *                    [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
 *                     ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 
 *                     ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr] 
 */
void gkyl_dg_calc_pkpm_vars_u_surf(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* pkpm_u, struct gkyl_array* pkpm_u_surf);

/**
 * Compute pkpm pressure p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij in the volume, order 2*p
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param bb Input array of volume expansion of magnetic field unit tensor
 * @param vlasov_pkpm_moms Input array of pkpm kinetic moments [rho, p_parallel, p_perp]
 * @param p_ij Output array of volume expansion of pressure tensor p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij
 */
void gkyl_dg_calc_pkpm_vars_pressure(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bb, const struct gkyl_array* vlasov_pkpm_moms, struct gkyl_array* p_ij);

/**
 * Compute pkpm acceleration variables. All quantities order 2*p
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param pkpm_u_surf Input array of surface expansions of flow velocity
 * @param pkpm_u Input array of volume expansion of flow velocity
 * @param prim Input array of volume expansion of primitive moments 
 *             [1/rho*div(p_par b), T_perp/m, m/T_perp, 3*Txx/m, 3*Tyy/m, 3*Tzz/m]
 * @param bb Input array of volume expansion of magnetic field unit tensor
 * @param div_b Input array of div(b)
 * @param nu Input array of collisionality
 * @param pkpm_lax Output array of surface expansion of Lax penalization lambda_i = |u_i| + sqrt(3.0*T_ii/m)
 * @param pkpm_accel Output arrary of pkpm acceleration variables ordered as:
 *        0: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
          1: bb_grad_u (bb : grad(u))
          2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
          3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2 nu)
 */
void gkyl_dg_calc_pkpm_vars_accel(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* pkpm_u_surf, const struct gkyl_array* pkpm_u, 
  const struct gkyl_array* prim, const struct gkyl_array* bb, 
  const struct gkyl_array* div_b, const struct gkyl_array* nu, 
  struct gkyl_array* pkpm_lax, struct gkyl_array* pkpm_accel);

/**
 * Compute integrated PKPM variables (rho, rhoux, rhouy, rhouz, rhoux^2, rhouy^2, rhouz^2, p_parallel, p_perp).
 * For diagnostic purposes, all quantities are expanded in the 2*p basis 
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param vlasov_pkpm_moms Input array of parallel-kinetic-perpendicular-moment kinetic moments [rho, p_parallel, p_perp]
 * @param pkpm_u Input array of volume expansion of flow velocity
 * @param int_pkpm_vars Output array of integrated variables (9 components)
 */
void gkyl_dg_calc_pkpm_integrated_vars(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, 
  struct gkyl_array* pkpm_int_vars);

/**
 * Compute pkpm model explicit source terms for order p momentum equation. 
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param qmem Input array of q/m*EM fields
 * @param vlasov_pkpm_moms Input array of parallel-kinetic-perpendicular-moment kinetic moments [rho, p_parallel, p_perp]
 * @param euler_pkpm Input array of parallel-kinetic-perpendicular-moment fluid variables [rho ux, rho uy, rho uz]
 * @param rhs Output increment to fluid variables
 */
void gkyl_dg_calc_pkpm_vars_explicit_source(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* qmem, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs);

/**
 * Construct PKPM variables for I/O. Computes the conserved fluid variables 
 * [rho, rho ux, rho uy, rho uz, 
 *  Pxx + rho ux^2, Pxy + rho ux uy, Pxz + rho ux uz, Pyy + rho uy^2, Pyz + rho uy uz, Pzz + rho uz^2]
 * And copies the pkpm primitive and acceleration variables into an array for output
 * [T_perp/m, m/T_perp, 1/rho div(p_par b), 1/rho p_perp div(b), bb : grad(u)]
 * For diagnostic purposes, all quantities are expanded in the 2*p basis 
 *
 * @param up Updater for computing pkpm variables 
 * @param conf_range Configuration space range
 * @param vlasov_pkpm_moms Input array of parallel-kinetic-perpendicular-moment kinetic moments [rho, p_parallel, p_perp]
 * @param pkpm_u Input array of volume expansion of flow velocity
 * @param p_ij Input pressure tensor p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij
 * @param prim Input array of primitive moments [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp]
 * @param pkpm_accel Input arrary of pkpm acceleration variables ordered as:
 *        0: p_perp_div_b (1/rho p_perp*div(b) = T_perp/m*div(b))
          1: bb_grad_u (bb : grad(u))
          2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
          3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2 nu)
 * @param fluid_io Output array of conserved fluid variables (10 components)
 * @param pkpm_vars_io Output array of pkpm variables, primitive and acceleration (5 components)
 */
void gkyl_dg_calc_pkpm_vars_io(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, const struct gkyl_array* p_ij, 
  const struct gkyl_array* prim, const struct gkyl_array* pkpm_accel, 
  struct gkyl_array* fluid_io, struct gkyl_array* pkpm_vars_io);

/**
 * Limit slopes in the PKPM system for order p momentum equation. 
 *
 * @param up               Updater for computing pkpm variables 
 * @param conf_range       Configuration space range
 * @param vlasov_pkpm_moms Input array of parallel-kinetic-perpendicular-moment kinetic moments [rho, p_parallel, p_perp]
 * @param pkpm_u           Input array of volume expansion of flow velocity
 * @param p_ij             Input pressure tensor p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij
 * @param euler_pkpm       Input (and Output after limiting) array of fluid variables [rho ux, rho uy, rho uz]
 */
void gkyl_dg_calc_pkpm_vars_limiter(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, const struct gkyl_array* p_ij, 
  struct gkyl_array* euler_pkpm);

/**
 * Delete pointer to updater to compute pkpm variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_pkpm_vars_release(struct gkyl_dg_calc_pkpm_vars *up);

/**
 * Host-side wrappers for pkpm vars operations on device
 */

void gkyl_dg_calc_pkpm_vars_advance_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* p_ij, 
  const struct gkyl_array* pkpm_div_ppar, const struct gkyl_array* div_b, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* prim);

void gkyl_dg_calc_pkpm_vars_surf_advance_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* p_ij, 
  struct gkyl_array* prim_surf);

void gkyl_dg_calc_pkpm_vars_u_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* pkpm_u);

void gkyl_dg_calc_pkpm_vars_u_surf_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* pkpm_u, struct gkyl_array* pkpm_u_surf);

void gkyl_dg_calc_pkpm_vars_pressure_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bb, const struct gkyl_array* vlasov_pkpm_moms, struct gkyl_array* p_ij);

void gkyl_dg_calc_pkpm_vars_accel_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* pkpm_u_surf, const struct gkyl_array* pkpm_u, 
  const struct gkyl_array* prim, const struct gkyl_array* bb, 
  const struct gkyl_array* div_b, const struct gkyl_array* nu, 
  struct gkyl_array* pkpm_lax, struct gkyl_array* pkpm_accel);

void gkyl_dg_calc_pkpm_integrated_vars_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, 
  struct gkyl_array* pkpm_int_vars);

void gkyl_dg_calc_pkpm_vars_explicit_source_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* qmem, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs);

void gkyl_dg_calc_pkpm_vars_io_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, const struct gkyl_array* p_ij, 
  const struct gkyl_array* prim, const struct gkyl_array* pkpm_accel, 
  struct gkyl_array* fluid_io, struct gkyl_array* pkpm_vars_io);

void gkyl_dg_calc_pkpm_vars_limiter_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, const struct gkyl_array* p_ij, 
  struct gkyl_array* euler_pkpm);

