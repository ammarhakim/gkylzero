#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_correct_MJ gkyl_correct_MJ;

/**
 * Create new updater to correct a Maxwellian to match specified
 * moments.
 *
 * @param grid Grid on which updater lives
 * @param conf_basis Conf space basis functions
 * @param phase_basis Phase space basis functions
 * @param conf_local_cells Number of cells in local config-space
 * @param conf_local_ext_cells Number of cells in local extended config-space
 */
gkyl_correct_MJ *gkyl_correct_MJ_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *conf_basis,
  const struct gkyl_basis *phase_basis,
  long conf_local_ncells, long conf_local_ext_ncells);

/**
 * Fix the Maxwellian so that it's moments match desired moments.
 *
 * @param cMJ Maxwell correction updater
 * @param fout Distribution function to fix (modified in-place)
 * @param m0 Desired number density
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 */
void gkyl_correct_MJ_fix(gkyl_correct_MJ *cMJ,
  struct gkyl_array *fout,
  const struct gkyl_array *m0,
  const struct gkyl_array *m1i,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local);

/**
 * Delete updater.
 *
 * @param cMJ Updater to delete.
 */
void gkyl_correct_MJ_release(gkyl_correct_MJ* cMJ);
