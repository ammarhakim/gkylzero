#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_moment gkyl_dg_updater_moment;

// return type for drag and diffusion timers
struct gkyl_dg_updater_moment_tm {
  double moment_tm; // time for moment updates
};

/**
 * Create new updater to compute moments of distribution function.
 * Supports Vlasov-Maxwell, special relativistic Vlasov-Maxwell,
 * and parallel-kinetic-perpendicular-moment (pkpm) Vlasov
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis function
 * @param conf_range Config space range
 * @param vel_range Velocity space range
 * @param model_id Enum identifier for model type (e.g., SR, PKPM, see gkyl_eqn_type.h)
 * @param is_integrated Boolean for if the moment is an integrated moment
 * @param mom Name of moment
 * @param mass Mass of species 
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 * 
 * @return New moment updater object
 */
struct gkyl_dg_updater_moment*
gkyl_dg_updater_moment_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range,
  enum gkyl_model_id model_id, const char *mom, 
  bool is_integrated, double mass, bool use_gpu);

/**
 * Acquire moment object
 *
 * @param moment moment updater object
 * 
 * @return moment type object
 */
struct gkyl_mom_type* 
gkyl_dg_updater_moment_acquire_type(const struct gkyl_dg_updater_moment* moment);

/**
 * Acquire number of moments
 *
 * @param moment moment updater object
 * 
 * @return number of moments
 */
int 
gkyl_dg_updater_moment_num_mom(const struct gkyl_dg_updater_moment* moment);

/**
 * Compute moment. The update_phase_rng and update_conf_rng MUST be a sub-range of the
 * be a sub-range of the range on which the array is defined. 
 * That is, it must be either the same range as the array range, or one created using the
 * or one created using the gkyl_sub_range_init method.
 *
 * @param moemnt moment updater object
 * @param update_phase_rng Phase space range on which to compute.
 * @param update_conf_rng Configuration space range on which to compute.
 * Auxiliary variables used by special relativistic vlasov solver
 * @param p_over_gamma p/gamma (velocity)
 * @param gamma gamma = sqrt(1 + p^2)
 * @param gamma_inv gamma_inv = 1/gamma = 1/sqrt(1 + p^2)
 * @param V_drift bulk fluid velocity (computed from M0*V_drift = M1i with weak division)
 * @param GammaV2 Gamma^2 = 1/(1 - V_drift^2/c^2), Lorentz boost factor squared from bulk fluid velocity
 * @param GammaV_inv Gamma_inv = sqrt(1 - V_drift^2/c^2), inverse Lorentz boost factor from bulk fluid velocity
 * @param fIn Input to updater
 * @param mout Output moment
 */
void
gkyl_dg_updater_moment_advance(struct gkyl_dg_updater_moment *moment,
  const struct gkyl_range *update_phase_rng, const struct gkyl_range *update_conf_rng,
  const struct gkyl_array *p_over_gamma, const struct gkyl_array *gamma, 
  const struct gkyl_array *gamma_inv, const struct gkyl_array *V_drift, 
  const struct gkyl_array *GammaV2, const struct gkyl_array *GammaV_inv, 
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT mout);

void
gkyl_dg_updater_moment_advance_cu(struct gkyl_dg_updater_moment *moment,
  const struct gkyl_range *update_phase_rng, const struct gkyl_range *update_conf_rng,
  const struct gkyl_array *p_over_gamma, const struct gkyl_array *gamma, 
  const struct gkyl_array *gamma_inv, const struct gkyl_array *V_drift, 
  const struct gkyl_array *GammaV2, const struct gkyl_array *GammaV_inv, 
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT mout);

/**
 * Return total time spent in computing moments
 *
 * @param moment Updater object
 * @return timers
 */
struct gkyl_dg_updater_moment_tm gkyl_dg_updater_moment_get_tm(const struct gkyl_dg_updater_moment *moment);

/**
 * Delete updater.
 *
 * @param moment Updater to delete.
 */
void gkyl_dg_updater_moment_release(struct gkyl_dg_updater_moment* moment);