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
 * @param grid configuration space grid.
 * @param bcs 0 for periodic 1 for non-periodic. Determines whether mapc2p extends to ghost cells
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_metric* gkyl_calc_metric_new(const struct gkyl_basis *cbasis,
  const struct gkyl_rect_grid *grid, const int *bcs, bool use_gpu);

/**
 * Use finite differences to calculate metric coefficients at nodes
 *
 * @param up calc_metric updater object.
 * @param nrange nodal range.
 * @param mc2p_nodal_fd nodal array containing cartesian coordinates at nodes and nearby nodes used for FD
 * @param gFld output field where metric coefficients will be placed
 */
void gkyl_calc_metric_advance(gkyl_calc_metric *up, struct gkyl_range *nrange, struct gkyl_array *mc2p_nodal_fd, double *dzc, struct gkyl_array *gFld, struct gkyl_range *update_range);
//void gkyl_calc_metric_advance(const gkyl_calc_metric *up, const struct gkyl_range *crange, struct gkyl_array *XYZ, struct gkyl_array *gFld);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_metric_release(gkyl_calc_metric* up);
