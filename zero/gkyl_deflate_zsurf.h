#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_deflate_zsurf gkyl_deflate_zsurf;

/**
 * Create new updater to deflate a 2d (x,z) modal expansion to a 1d (x) modal expansion
 *  at z=+/1 surfaces given a z index
 * @param cbasis Basis object (configuration space).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */


gkyl_deflate_zsurf* gkyl_deflate_zsurf_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *deflated_cbasis,
  const struct gkyl_rect_grid *grid, const struct gkyl_rect_grid *deflated_grid, int edge, bool use_gpu);

/**
 * Advance deflate_zsurf (compute the derived_zsurf coefficients).
 *
 * @param up deflate_zsurf updater object.
 * @param crange Config-space range.
 * @param fld 2d (x,z) field containing DG rep 
 * @param deflated_fld 1d (x) field containing DG rep
 */

void gkyl_deflate_zsurf_advance(const gkyl_deflate_zsurf *up, int zidx, const struct gkyl_range *range, const struct gkyl_range* deflated_range, const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_deflate_zsurf_release(gkyl_deflate_zsurf* up);
