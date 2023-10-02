#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_geo_gyrokinetic.h>

// Object type
typedef struct gkyl_calc_bmag gkyl_calc_bmag;
typedef struct bmag_ctx bmag_ctx;

/**
 * Create new updater to compute the metric coefficients
 *
 * @param cbasis Basis object (configuration space).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_bmag* 
gkyl_calc_bmag_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_rect_grid *cgrid, const struct gkyl_rect_grid *pgrid, const gkyl_geo_gyrokinetic *app, const struct gkyl_geo_gyrokinetic_geo_inp *ginp, bool use_gpu);


/**
 * Advance calc_metric (compute the metric coefficients).
 *
 * @param up calc_metric updater object.
 * @param crange Config-space range.
 * @param XYZ field containing DG rep of cartesian coordinates
 * @param gFld output field where metric coefficients will be placed
 */

void gkyl_calc_bmag_advance(const gkyl_calc_bmag *up, const struct gkyl_range *crange, const struct gkyl_range *crange_ext,
     const struct gkyl_range *prange, const struct gkyl_range *prange_ext, struct gkyl_array *psidg, struct gkyl_array *psibyrdg, struct gkyl_array *psibyr2dg, struct gkyl_array *bphidg, struct gkyl_array* bmag_compdg, struct gkyl_array* mapc2p);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_bmag_release(gkyl_calc_bmag* up);
