#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_correct_maxwellian gkyl_correct_maxwellian;

/**
 * Create new updater to correct a Maxwellian to match specified
 * moments.
 *
 * @param grid Grid on which updater lives
 * @param conf_basis Configuration space basis functions
 * @param conf_local_ext_cells Number of cells in local extended configuration space
 */
gkyl_correct_maxwellian *gkyl_correct_maxwellian_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, long conf_local_ext_ncells);

/**
 * Delete updater.
 *
 * @param cmax Updater to delete.
 */
void gkyl_correct_maxwellian_release(gkyl_correct_maxwellian* cmax);
