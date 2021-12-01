#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

/** Flags for indicating acting edge of velocity space */
enum gkyl_vel_edge { GKYL_VX_LOWER, GKYL_VY_LOWER, GKYL_VZ_LOWER, GKYL_VX_UPPER, GKYL_VY_UPPER, GKYL_VZ_UPPER };

// Object type
typedef struct gkyl_mom_bcorr gkyl_mom_bcorr;

struct gkyl_mom_bcorr {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_basis; // number of basis functions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  const struct gkyl_mom_type *momt; // equation object
};

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
gkyl_mom_bcorr* gkyl_mom_bcorr_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt);

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
void gkyl_mom_bcorr_advance(gkyl_mom_bcorr *bcorr, struct gkyl_range phase_rng, struct gkyl_range conf_rng,
  const struct gkyl_array *fIn, struct gkyl_array *out);
  
/**
 * Delete updater.
 *
 * @param hdg Updater to delete.
 */
void gkyl_mom_bcorr_release(gkyl_mom_bcorr* bcorr);
