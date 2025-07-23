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
 * @param global, _ext computational global Config-space range and extended range.
 * @param local, _ext computational local Config-space range and extended range.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_metric* gkyl_calc_metric_new(const struct gkyl_basis *cbasis,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *global, const struct gkyl_range *global_ext, 
  const struct gkyl_range *local, const struct gkyl_range *local_ext, bool use_gpu);

/**
 * Use finite differences to calculate metric coefficients and tangent vectors at nodes
 * Then convert to modal
 *
 * @param up calc_metric updater object.
 * @param nrange nodal range.
 * @param mc2p_nodal_fd nodal array containing cartesian coordinates at nodes and nearby nodes (epsilon and 2 epsilon away) used for FD
 *    At each array location 39 values are stored.
 *    The arrangement is as follows: X_c, Y_c, Z_c, 
 *    X_L1, Y_L1, Z_L1, X_R1, Y_R1, Z_R1,
 *    X_L2, Y_L2, Z_L2, X_R2, Y_R2, Z_R2,
 *    X_L3, Y_L3, Z_L3, X_R3, Y_R3, Z_R3,
 *    X_LL1, Y_LL1, Z_LL1, X_RR1, Y_RR1, Z_RR1,
 *    X_LL2, Y_LL2, Z_LL2, X_RR2, Y_RR2, Z_RR2,
 *    X_LL3, Y_LL3, Z_LL3, X_RR3, Y_RR3, Z_RR3
 *    where L#/R# indicates a node shifted to the left/right by epsilon in coordinate #
 *    and LL#/RR# indicates a node shifted to the left/right by 2 epsilon in coordinate #
 * @param dzc epsilons used for FD
 * @param gFld output field where metric modal coefficients will be placed
 * @param tanvecFld output field where tangent vector modal coefficients will be placed
 * @param dualFld output field where dual vector modal coefficients will be placed
 * @param dualmagFld output field where magnitude of the dual vectors will be placed
 * @param normFld output field where dual vector modal coefficients will be placed
 * @param update range. Modal range over which metric coefficients and tangent vectors will be calculated
 */
void gkyl_calc_metric_advance(gkyl_calc_metric *up, struct gkyl_range *nrange, struct gkyl_array *mc2p_nodal_fd, 
  double *dzc, struct gkyl_array *gFld, struct gkyl_array *tanvecFld, struct gkyl_array *dualFld, 
  struct gkyl_array *dualmagFld, struct gkyl_array* normFld, const struct gkyl_range *update_range);


/**
 * Use finite differences to calculate metric coefficients and jacobian at nodes
 * Using the explicit in Eq. 66-73 of the GK coordinates document
 * Then convert to modal
 *
 * @param up calc_metric updater object.
 * @param nrange nodal range.
 * @param mc2p_nodal_fd nodal array containing cylindrical coordinates at nodes and nearby nodes used for FD
 * @param ddtheta_nodal input nodal array containing dR/dtheta, dZ/dtheta, and dphi/dtheta at nodes
 * @param bmag_nodal input nodal array containing B at nodes
 * @param dzc epsilons used for FD
 * @param gFld output field where metric modal coefficients will be placed
 * @param tanvecFld output field where tangent vector modal coefficients will be placed
 * @param dualFld output field where dual vector modal coefficients will be placed
 * @param dualmagFld output field where magnitude of the dual vectors will be placed
 * @param normFld output field where dual vector modal coefficients will be placed
 * @param jFld output field where jacobian modal coefficients will be placed
 * @param bcartFld output field where cartesian compnents of b modal coefficients will be placed
 * @param update range. Modal range over which metric coefficients and tangent vectors will be calculated
 */
void gkyl_calc_metric_advance_rz(gkyl_calc_metric *up, struct gkyl_range *nrange, struct gkyl_array *mc2p_nodal_fd, 
  struct gkyl_array *ddtheta_nodal, struct gkyl_array *bmag_nodal, double *dzc, struct gkyl_array *gFld, 
  struct gkyl_array *tanvecFld, struct gkyl_array *dualFld, struct gkyl_array *dualmagFld, struct gkyl_array *normFld, 
  struct gkyl_array *jFld, struct gkyl_array* bcartFld, const struct gkyl_range *update_range);

/**
 * Use finite differences to calculate metric coefficients and jacobian at nodes
 * Using the explicit in Eq. 66-73 of the GK coordinates document
 * This function uses the coordinate definition alpha = phi
 * Then convert to modal
 *
 * @param up calc_metric updater object.
 * @param nrange nodal range.
 * @param mc2p_nodal_fd nodal array containing cylindrical coordinates at nodes and nearby nodes used for FD
 * @param ddtheta_nodal input nodal array containing dphi/dtheta, dR/dtheta, and dZ/dtheta at nodes
 * @param dzc epsilons used for FD
 * @param gFld output field where covariant metric modal coefficients will be placed
 * @param grFld output field where contravariant metric modal coefficients will be placed
 * @param update range. Modal range over which metric coefficients and tangent vectors will be calculated
 */
void gkyl_calc_metric_advance_rz_neut( gkyl_calc_metric *up, struct gkyl_range *nrange,
    struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *ddtheta_nodal, double *dzc, 
    struct gkyl_array *gFld, struct gkyl_array *grFld, const struct gkyl_range *update_range);

/**
 * Use finite differences to calculate metric coefficients and jacobian at nodes
 * Special analytical simplifications for mirrors
 * Then convert to modal
 *
 * @param up calc_metric updater object.
 * @param nrange nodal range.
 * @param mc2p_nodal_fd nodal array containing cylindrical coordinates at nodes and nearby nodes used for FD
 * @param ddtheta_nodal input nodal array containing dR/dtheta, dZ/dtheta, and dphi/dtheta at nodes
 * @param bmag_nodal input nodal array containing B at nodes
 * @param dzc epsilons used for FD
 * @param gFld output field where metric modal coefficients will be placed
 * @param tanvecFld output field where tangent vector modal coefficients will be placed
 * @param dualFld output field where dual vector modal coefficients will be placed
 * @param dualmagFld output field where magnitude of the dual vectors will be placed
 * @param normFld output field where dual vector modal coefficients will be placed
 * @param jFld output field where jacobian modal coefficients will be placed
 * @param bcartFld output field where cartesian compnents of b modal coefficients will be placed
 * @param update range. Modal range over which metric coefficients and tangent vectors will be calculated
 */
void gkyl_calc_metric_advance_mirror(gkyl_calc_metric *up, struct gkyl_range *nrange, struct gkyl_array *mc2p_nodal_fd,
   struct gkyl_array *ddtheta_nodal, struct gkyl_array *bmag_nodal, double *dzc, struct gkyl_array *gFld, 
   struct gkyl_array *tanvecFld, struct gkyl_array *dualFld, struct gkyl_array *dualmagFld, struct gkyl_array *normFld,
    struct gkyl_array *jFld, struct gkyl_array* bcartFld, const struct gkyl_range *update_range);

/**
 * Calculate cartesian components of bhat
 * @param up calc_metric updater object.
 * @param nrange nodal range.
 * @param biFld input field containing b_i DG expansion
 * @param dualFld input field containing dual vectors DG expansion
 * @param bcartFld output field containing DG expansion of cartesian components of bhat
 */
void gkyl_calc_metric_advance_bcart(gkyl_calc_metric *up, struct gkyl_range *nrange, struct gkyl_array *biFld,
   struct gkyl_array *dualFld, struct gkyl_array *bcartFld, const struct gkyl_range *update_range);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_metric_release(gkyl_calc_metric* up);
