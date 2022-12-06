#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_collisions gkyl_dg_updater_collisions;

// return type for drag and diffusion timers
struct gkyl_dg_updater_lbo_gyrokinetic_tm {
  double drag_tm, diff_tm; // time for drag and diffusion updates
};

/**
 * Create new updater to update lbo equations using hyper dg.
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis function
 * @param conf_range Config space range
 * @param mass Species mass
 * @return New LBO updater object
 */
struct gkyl_dg_updater_collisions* 
gkyl_dg_updater_lbo_gyrokinetic_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, double mass, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param lbo LBO updater object
 * @param update_rng Range on which to compute.
 * @param bmag Magnitude of magnetic field
 * @param nu_sum Sum of coll freq
 * @param nu_prim_moms Sum of coll freq*u and freq*vtsq
 * @param m2self 2nd velocity moment of this species.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_lbo_gyrokinetic_advance(struct gkyl_dg_updater_collisions *lbo,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *bmag_inv,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_prim_moms,
  const struct gkyl_array *m2self, 
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

void gkyl_dg_updater_lbo_gyrokinetic_advance_cu(struct gkyl_dg_updater_collisions *lbo,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *bmag_inv,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_prim_moms,
  const struct gkyl_array *m2self, 
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in drag and diffusion terms
 *
 * @param lbo Updater object
 * @return timers
 */
struct gkyl_dg_updater_lbo_gyrokinetic_tm gkyl_dg_updater_lbo_gyrokinetic_get_tm(const struct gkyl_dg_updater_collisions *coll);

/**
 * Delete updater.
 *
 * @param lbo Updater to delete.
 */
void gkyl_dg_updater_lbo_gyrokinetic_release(struct gkyl_dg_updater_collisions* coll);
