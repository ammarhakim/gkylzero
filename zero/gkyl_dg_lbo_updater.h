#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov_lbo.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_lbo_updater gkyl_dg_lbo_updater;

/**
 * Create new updater to update lbo equations using hyper dg.
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis function
 * @param conf_range Config space range
 * @return New LBO updater object
 */
gkyl_dg_lbo_updater* gkyl_dg_lbo_updater_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, const struct gkyl_range *conf_range);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param lbo LBO updater object
 * @param update_rng Range on which to compute.
 * @param nu_sum Sum of coll freq
 * @param nu_u Sum of coll freq*u
 * @param nu_vthsq Sum of coll freq*vth
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_lbo_updater_advance(gkyl_dg_lbo_updater *lbo, struct gkyl_range update_rng,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_u, const struct gkyl_array *nu_vthsq,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs);
  
/**
 * Delete updater.
 *
 * @param lbo Updater to delete.
 */
void gkyl_dg_lbo_updater_release(gkyl_dg_lbo_updater* lbo);
