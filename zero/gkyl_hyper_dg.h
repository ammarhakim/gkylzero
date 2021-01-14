#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_hyper_dg gkyl_hyper_dg;

/**
 * Create new updater to update equations using DG algorithm.
 *
 * @param grid Grid object
 * @param basis Basis functions
 * @param equation Equation object
 * @param num_up_dirs Number of directions to update
 * @param update_dirs List of directions to update (size 'num_up_dirs')
 * @param zero_flux_flags Flags with zero-flux BCs (size 'num_up_dirs')
 * @param update_vol_term Set to 0 to skip volume update
 */
gkyl_hyper_dg* gkyl_hyper_dg_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation,
  int num_up_dirs, int update_dirs[], int zero_flux_flags[], int update_vol_term);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param hdg Hyper DG updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_hyper_dg_advance(gkyl_hyper_dg *hdg, const struct gkyl_range *update_rng,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs);
  
/**
 * Delete updater.
 *
 * @param hdg Updater to delete.
 */
void gkyl_hyper_dg_release(gkyl_hyper_dg* hdg);
