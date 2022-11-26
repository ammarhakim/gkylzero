#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_MJ_moments gkyl_MJ_moments;

/**
 * Create new updater to correct a Maxwellian to match specified
 * moments.
 *
 * @param grid Grid on which updater lives
 * @param conf_basis Conf space basis functions
 * @param phase_basis Phase space basis functions
 * @param conf_range configuration space range
 * @param vel_range velocity space range
 * @param conf_local_cells Number of cells in local config-space
 * @param conf_local_ext_cells Number of cells in local extended config-space
 */
gkyl_MJ_moments *gkyl_MJ_moments_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *conf_basis,
  const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_range,
  const struct gkyl_range *vel_range,
  long conf_local_ncells, long conf_local_ext_ncells, bool use_gpu);

/**
 * Fix the Maxwellian so that it's moments match desired moments.
 *
 * @param cMJ Maxwell correction updater 
 * @param p_over_gamma velocity array
 * @param gamma array
 * @param gamma_inv array
 * @param fout Distribution function to fix (modified in-place)
 * @param m0 Desired number density
 * @param m1i Desired velocity
 * @param m2 Desired Temperature
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 */
void gkyl_MJ_moments_fix(gkyl_MJ_moments *cMJ, const struct gkyl_array *p_over_gamma,
  const struct gkyl_array *gamma, const struct gkyl_array *gamma_inv,
  struct gkyl_array *fout,
  struct gkyl_array *m0,
  struct gkyl_array *m1i,
  struct gkyl_array *m2,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local);

/**
 * Delete updater.
 *
 * @param cMJ Updater to delete.
 */
void gkyl_MJ_moments_release(gkyl_MJ_moments* cMJ);
