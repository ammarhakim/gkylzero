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
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_metric* gkyl_calc_metric_new(const struct gkyl_basis *cbasis,
  const struct gkyl_rect_grid *grid, bool use_gpu);

/**
 * Use finite differences to calculate metric coefficients and tangent vectors at nodes
 * Then convert to modal
 *
 * @param up calc_metric updater object.
 * @param nrange nodal range.
 * @param mc2p_nodal_fd nodal array containing cartesian coordinates at nodes and nearby nodes used for FD
 * @param gFld output field where metric modal coefficients will be placed
 * @param tanvecFld output field where tangent vector modal coefficients will be placed
 * @param update range. Modal range over which metric coefficients and tangent vectors will be calculated
 */
void gkyl_calc_metric_advance(gkyl_calc_metric *up, struct gkyl_range *nrange, struct gkyl_array *mc2p_nodal_fd, double *dzc, struct gkyl_array *gFld, struct gkyl_array *tanvecFld, struct gkyl_array *dualFld, const struct gkyl_range *update_range);


/**
 * Use finite differences to calculate metric coefficients and jacobian at nodes
 * Using the explicit in Eq. 66-73 of the GK coordinates document
 * Then convert to modal
 *
 * @param up calc_metric updater object.
 * @param nrange nodal range.
 * @param mc2p_nodal_fd nodal array containing cylindrical coordinates at nodes and nearby nodes used for FD
 * @param gFld output field where metric modal coefficients will be placed
 * @param jFld output field where jacobian modal coefficients will be placed
 * @param update range. Modal range over which metric coefficients and tangent vectors will be calculated
 */
void gkyl_calc_metric_advance_rz(gkyl_calc_metric *up, struct gkyl_range *nrange, struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *dphidtheta_nodal, struct gkyl_array *bmag_nodal, double *dzc, struct gkyl_array *gFld, struct gkyl_array *jFld, const struct gkyl_range *update_range);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_metric_release(gkyl_calc_metric* up);
