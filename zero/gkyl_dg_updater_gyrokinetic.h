#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_gyrokinetic gkyl_dg_updater_gyrokinetic;

// Boundary condition types.
enum gkyl_gkeqn_id {
  GKYL_GK_DEFAULT=0,    // Default GK equation (kernels).
//  GKYL_GK_SOFTMIRROR,
};

/**
 * Create new updater to update gyrokinetic equations using hyper dg.
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis function
 * @param conf_range Config space range
 * @param is_zero_flux_dir True in directions with (lower and upper) zero flux BCs.
 * @param gkyl_gkeqn_id Enum identifier for gyrokinetic equation type.
 * 
 * @return New gyrokinetic updater object
 */
gkyl_dg_updater_gyrokinetic* gkyl_dg_updater_gyrokinetic_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const bool *is_zero_flux_dir,
  enum gkyl_gkeqn_id eqn_id, double charge, double mass, bool use_gpu);

/**
 * Acquire gyrokinetic equation object
 *
 * @param up gyrokinetic updater object
 * 
 * @return gyrokinetic equation object
 */
struct gkyl_dg_eqn* 
gkyl_dg_updater_gyrokinetic_acquire_eqn(const gkyl_dg_updater_gyrokinetic* up);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up gyrokinetic updater object
 * @param field_id Enum identifier for field type (see gkyl_eqn_type.h)
 * @param update_rng Range on which to compute.
 * @param bmag Magnetic field magnitude.
 * @param jacobtot_inv Reciprocal of the total (conf * guiding center) Jacobian.
 * @param cmag Multiplicative factor in Clebsch-like form of the magnetic field.
 * @param b_i Covariant components of hat{b}.
 * @param phi Electrostatic potential.
 * @param apar Parallel component of magnetic vector potential.
 * @param apardot Time rate of change of apar.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_gyrokinetic_advance(gkyl_dg_updater_gyrokinetic *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *bmag, const struct gkyl_array *jacobtot_inv,
  const struct gkyl_array *cmag, const struct gkyl_array *b_i,
  const struct gkyl_array *phi, const struct gkyl_array *apar,
  const struct gkyl_array *apardot, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_updater_gyrokinetic_release(gkyl_dg_updater_gyrokinetic* up);
