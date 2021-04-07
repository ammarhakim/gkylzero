#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_wave_prop gkyl_wave_prop;

/**
 * Create new updater to update equations using wave-propagation
 * algorithm.
 *
 * @param grid Grid object
 * @param equation Equation object
 * @param num_up_dirs Number of directions to update
 * @param update_dirs List of directions to update (size 'num_up_dirs')
 */
gkyl_wave_prop* gkyl_wave_prop_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_wv_eqn *equation, int num_up_dirs, int update_dirs[]);

/**
 * Compute wave-propagation update. The update_rng MUST be a sub-range
 * of the range on which the array is defined. That is, it must be
 * either the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param wv Updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_wave_prop_advance(const gkyl_wave_prop *wv, const struct gkyl_range *update_rng,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs);
  
/**
 * Delete updater.
 *
 * @param wv Updater to delete.
 */
void gkyl_wave_prop_release(gkyl_wave_prop* wv);
