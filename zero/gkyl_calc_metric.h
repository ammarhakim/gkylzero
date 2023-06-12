#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_calc_metric gkyl_calc_metric;

/**
 * Create new updater to compute the metric coefficients
 *
 * @param cbasis Basis object (configuration space).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_metric* gkyl_calc_metric_new(const struct gkyl_basis *cbasis,
  struct gkyl_rect_grid grid, bool use_gpu);

/**
 * Advance calc_metric (compute the metric coefficients).
 *
 * @param up calc_metric updater object.
 * @param crange Config-space range.
 * @param XYZ field containing DG rep of cartesian coordinates
 * @param gFld output field where metric coefficients will be placed
 */

void gkyl_calc_metric_advance(const gkyl_calc_metric *up, const struct gkyl_range *crange,
    struct gkyl_array *XYZ, struct gkyl_array *gFld);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_metric_release(gkyl_calc_metric* up);
