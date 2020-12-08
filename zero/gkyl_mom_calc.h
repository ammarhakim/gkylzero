#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_mom_calc gkyl_mom_calc;

/**
 * Create new updater to compute moments of distribution
 * function. Free using gkyl_mom_calc_new_release.
 *
 * @param grid Grid object
 * @param momt Pointer to moment type object
 * @return New updater pointer.
 */
gkyl_mom_calc* gkyl_mom_calc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt);

/**
 * Compute moment of distribution function. The update_rng MUST be a
 * sub-range of the range on which the distribution function is
 * defined. That is, it must be either the same range as the
 * distribution function range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param calc Moment calculator updater to run
 * @param update_rng Range on which to run updater
 * @param fin Input distribution function array
 * @param mout Output moment array
 */
void gkyl_mom_calc_advance(const gkyl_mom_calc* calc,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *fin, struct gkyl_array *mout);

/**
 * Delete pointer to moment calculator updater.
 *
 * @param calc Updater to delete.
 */
void gkyl_mom_calc_release(gkyl_mom_calc* calc);
