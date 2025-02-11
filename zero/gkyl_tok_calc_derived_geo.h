#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_tok_calc_derived_geo gkyl_tok_calc_derived_geo;

/**
 * Create new updater to compute the derived_geo coefficients
 *
 * @param cbasis Basis object (configuration space).
 * @param grid computational grid
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_tok_calc_derived_geo* gkyl_tok_calc_derived_geo_new(const struct gkyl_basis *cbasis,
  const struct gkyl_rect_grid *grid, bool use_gpu);

/**
 * Advance tok_calc_derived_geo (compute the derived_geo coefficients).
 *
 * @param up tok_calc_derived_geo updater object.
 * @param crange Config-space range.
 * @param gFld input field containing DG rep of the metric coefficients
 * @param bmagFld input field containing DG rep of B = |B|
 * @param jFld input field with DG rep of jacobian (J)
 * @param jinvFld output field with DG rep of 1/J
 * @param grFld output field with DG rep of  g^ij
 * @param jtotinvFld output field with DG rep of JB
 * @param bamginvFld output field with DG rep of 1/B
 * @param bmaginvsqFld output field with DG rep of 1/B^2
 * @param gxxJFld output field with DG rep of Jg^xx
 * @param gxyJFld output field with DG rep of Jg^xy
 * @param gyyJFld output field with DG rep of Jg^yy
 * @param gxzJFld output field with DG rep of Jg^xz
 * @param eps2Fld output field with DG rep of eps2 = Jg^33 - J/g_33
 */

void gkyl_tok_calc_derived_geo_advance(const gkyl_tok_calc_derived_geo *up, const struct gkyl_range *crange,
    struct gkyl_array *gFld, struct gkyl_array *bmagFld, struct gkyl_array *jFld, struct gkyl_array *jinvFld,
    struct gkyl_array *grFld, struct gkyl_array *biFld, struct gkyl_array *cmagFld, struct gkyl_array *jtotFld, struct gkyl_array *jtotinvFld, struct gkyl_array *bmaginvFld, struct gkyl_array *bmaginvsqFld, struct gkyl_array *gxxJFld,  struct gkyl_array *gxyJFld, struct gkyl_array *gyyJFld, struct gkyl_array *gxzJFld, struct gkyl_array *eps2Fld);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_tok_calc_derived_geo_release(gkyl_tok_calc_derived_geo* up);
