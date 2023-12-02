#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_calc_derived_geo gkyl_calc_derived_geo;

/**
 * Create new updater to compute the derived_geo coefficients
 *
 * @param cbasis Basis object (configuration space).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_derived_geo* gkyl_calc_derived_geo_new(const struct gkyl_basis *cbasis,
  const struct gkyl_rect_grid *grid, bool use_gpu);

/**
 * Advance calc_derived_geo (compute the derived_geo coefficients).
 *
 * @param up calc_derived_geo updater object.
 * @param crange Config-space range.
 * @param gFld field containing DG rep of the metric coefficients
 * @param jFld output field where jacobian will be placed
 */

void gkyl_calc_derived_geo_advance(const gkyl_calc_derived_geo *up, const struct gkyl_range *crange,
    struct gkyl_array *gFld, struct gkyl_array *bmagFld, struct gkyl_array *jFld, struct gkyl_array *jinvFld,
    struct gkyl_array *grFld, struct gkyl_array *biFld, struct gkyl_array *cmagFld, struct gkyl_array *jtotFld, struct gkyl_array *jtotinvFld, struct gkyl_array *bmaginvFld, struct gkyl_array *bmaginvsqFld, struct gkyl_array *gxxJFld,  struct gkyl_array *gxyJFld, struct gkyl_array *gyyJFld);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_derived_geo_release(gkyl_calc_derived_geo* up);
