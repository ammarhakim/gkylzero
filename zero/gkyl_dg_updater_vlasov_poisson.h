#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_vlasov_poisson gkyl_dg_updater_vlasov_poisson;

// return type for vlasov timers
struct gkyl_dg_updater_vlasov_poisson_tm {
  double vlasov_poisson_tm; // time for vlasov updates
};

/**
 * Create new updater to update Vlasov-Poisson equations using hyper dg.
 * Supports Vlasov-Poisson with and without external potentials (phi_ext,A_ext).
 *
 * @param grid Grid object.
 * @param cbasis Configuration space basis functions.
 * @param pbasis Phase space basis function.
 * @param conf_range Configuration space range.
 * @param vel_range Velocity space range.
 * @param phase_range Phase space range.
 * @param is_zero_flux_dir True in directions with (lower and upper) zero flux BCs.
 * @param model_id Enum identifier for model type (e.g., General Geometry, see gkyl_eqn_type.h).
 * @param field_id Enum identifier for field type (e.g., Poisson, Poisson with external potentials, see gkyl_eqn_type.h).
 * @param aux_inp Void pointer to auxiliary fields. Void to be flexible to different auxfields structs.
 * @param use_gpu Boolean to determine whether struct objects are on host or device.
 * 
 * @return New vlasov_poisson updater object.
 */
gkyl_dg_updater_vlasov_poisson* gkyl_dg_updater_vlasov_poisson_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  const bool *is_zero_flux_dir, enum gkyl_vpmodel_id model_id, enum gkyl_vpfield_id field_id, void *aux_inp, bool use_gpu);

/**
 * Acquire Vlasov-Poisson equation object
 *
 * @param up Vlasov-Poisson updater object
 * 
 * @return Vlasov-Poisson equation object
 */
struct gkyl_dg_eqn* 
gkyl_dg_updater_vlasov_poisson_acquire_eqn(const gkyl_dg_updater_vlasov_poisson* up);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up Vlasov-Poisson updater object.
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater.
 * @param cflrate CFL scalar rate (frequency) array.
 * @param rhs RHS output.
 */
void gkyl_dg_updater_vlasov_poisson_advance(gkyl_dg_updater_vlasov_poisson *up,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in vlasov equation
 *
 * @param up Updater object.
 * @return timers.
 */
struct gkyl_dg_updater_vlasov_poisson_tm gkyl_dg_updater_vlasov_poisson_get_tm(const gkyl_dg_updater_vlasov_poisson *up);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_updater_vlasov_poisson_release(gkyl_dg_updater_vlasov_poisson* up);