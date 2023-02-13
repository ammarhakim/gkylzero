#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_updater_vlasov.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_bflux_vlasov gkyl_dg_updater_bflux_vlasov;

// return type for drag and diffusion timers
struct gkyl_dg_updater_bflux_vlasov_tm {
  double bflux_tm; // time for bflux updates
};

/**
 * Create new updater to update vlasov equations using hyper dg.
 * Supports Vlasov-Maxwell, Vlasov-Poisson (with and without vector potential A)
 * and special relativistic Vlasov-Maxwell
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis function
 * @param conf_range Config space range
 * @param vel_range Velocity space range
 * @param phase_range Phase space range
 * @param model_id Enum identifier for model type (e.g., SR, General Geometry, PKPM, see gkyl_eqn_type.h)
 * @param field_id Enum identifier for field type (e.g., Maxwell's, Poisson, see gkyl_eqn_type.h)
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 * 
 * @return New vlasov updater object
 */
gkyl_dg_updater_bflux_vlasov* gkyl_dg_updater_bflux_vlasov_new(const struct gkyl_rect_grid *grid, 
  int cdim, const gkyl_dg_updater_vlasov *vlasov, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param vlasov vlasov updater object
 * @param update_rng Range on which to compute.
 * @param aux1 Auxiliary field 1 (usually field, i.e., q/m*EM or q/m*phi)
 * @param aux2 Auxiliary field 2
 * @param aux3 Auxiliary field 3
 * @param aux4 Auxiliary field 4
 * @param aux5 Auxiliary field 5
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_bflux_vlasov_advance(gkyl_dg_updater_bflux_vlasov *bflux,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT rhs);

void gkyl_dg_updater_bflux_vlasov_advance_cu(gkyl_dg_updater_bflux_vlasov *bflux,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in vlasov equation
 *
 * @param vlasov Updater object
 * @return timers
 */
struct gkyl_dg_updater_bflux_vlasov_tm gkyl_dg_updater_bflux_vlasov_get_tm(const gkyl_dg_updater_bflux_vlasov *bflux);

/**
 * Delete updater.
 *
 * @param bflux Updater to delete.
 */
void gkyl_dg_updater_bflux_vlasov_release(gkyl_dg_updater_bflux_vlasov* bflux);
