#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_efit gkyl_efit;

/**
 * Create new updater to compute the derived_geo coefficients
 *
 * @param rzbasis Basis object 
 * @param rz grid to be filled from efit
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_efit* gkyl_efit_new(char *filepath, const struct gkyl_basis *rzbasis,
  struct gkyl_rect_grid *rzgrid, struct gkyl_range *rzlocal, struct gkyl_range *rzlocal_ext, bool use_gpu);

void gkyl_efit_release(gkyl_efit* up);
