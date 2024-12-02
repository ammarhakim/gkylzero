#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_array_average gkyl_array_average;

  // The symbol used for the coordinates are not related
  // to 2x or 3x runs, x is always the first dim, y the second, z the third.
  // E.g.: in a 2x run, y here translates into the z coordinate
enum gkyl_array_average_op {
    GKYL_ARRAY_AVERAGE_OP = 0,    // average over all dimensions
    GKYL_ARRAY_AVERAGE_OP_X,      // 1st dimension will remain
    GKYL_ARRAY_AVERAGE_OP_Y,      // 2nd dimension will remain
    GKYL_ARRAY_AVERAGE_OP_Z,      // 3rd dimension will remain
    GKYL_ARRAY_AVERAGE_OP_XY,     // 1st and 2nd dimensions will remain
    GKYL_ARRAY_AVERAGE_OP_XZ,     // 1st and 3rd dimensions will remain
    GKYL_ARRAY_AVERAGE_OP_YZ      // 2nd and 3rd dimensions will remain
};

/**
 * Create a new updater that computes the average of a gkyl_array, with an option to
 * set the dimensions to average
 *
 * @param grid computational grid to compute the surface element
 * @param basis gkyl basis, indicate the dimensionality ndim and polynomial type and order
 * @param op enum type to describe which kind of average has to be perform
 * @param use_gpu Indicate whether to perform integral on the device.
 */
struct gkyl_array_average*
gkyl_array_average_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis, 
  enum gkyl_array_average_op op, bool use_gpu);

/**
 * Compute the array average.
 *
 * @param up array_integrate updater.
 * @param full_rng range of the input array (all dimensions)
 * @param sub_rng range of the output array (only non-avg dimensions)
 * @param fin Input gkyl_array
 * @param avgout Output gkyl_array
 */
void gkyl_array_average_advance(gkyl_array_average *up, 
  const struct gkyl_range *full_rng, const struct gkyl_range *sub_rng,
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT avgout);

/**
 * Release memory associated with this updater.
 *
 * @param up array_integrate updater.
 */
void gkyl_array_average_release(gkyl_array_average *up);
