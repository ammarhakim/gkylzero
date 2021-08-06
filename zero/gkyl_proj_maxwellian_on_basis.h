#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_proj_maxwellian_on_basis gkyl_proj_maxwellian_on_basis;

/**
 * Create new updater to project Maxwellian on basis functions. Free
 * using gkyl_proj_maxwellian_on_basis_release method.
 *
 * @param grid Grid object
 * @param conf_basis Conf-space basis functions
 * @param phase_basis Phase-space basis functions
 * @param num_quad Number of quadrature nodes
 * @return New updater pointer.
 */
gkyl_proj_maxwellian_on_basis* gkyl_proj_maxwellian_on_basis_new(
  const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  int num_quad);

/**
 * Compute projection on basis. The update_rng MUST be a sub-range of
 * the range on which the array is defined. That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param pob Project on basis updater to run
 * @param update_rng Range on which to run projection.
 * @param out Output array
 */
void gkyl_proj_maxwellian_on_basis_advance(const gkyl_proj_maxwellian_on_basis *pob,
  const struct gkyl_range *update_rng, struct gkyl_array *out);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_maxwellian_on_basis_release(gkyl_proj_maxwellian_on_basis* pob);
