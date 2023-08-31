#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_correct_maxwellian_gyrokinetic gkyl_correct_maxwellian_gyrokinetic;

/**
 * Create new updater to correct a Maxwellian to match specified
 * moments.
 *
 * @param grid Grid on which updater lives
 * @param conf_basis Conf space basis functions
 * @param phase_basis Phase space basis functions
 * @param conf_local Local configuration space range
 * @param conf_local_ext Local extended configuration space range
 * @param bmag Magnetic field
 * @param jacob_tot Jacobian to project the maxwellian distribution
 * @param mass Mass of the species
 * @param use_gpu Bool to determine if on GPU
 */
gkyl_correct_maxwellian_gyrokinetic *gkyl_correct_maxwellian_gyrokinetic_new(
  const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *conf_local, const struct gkyl_range *conf_local_ext,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, double mass, bool use_gpu);

/**
 * Fix the Maxwellian so that it's moments match desired moments.
 *
 * @param cmax Maxwellian-fix updater
 * @param fM Distribution function to fix (modified in-place)
 * @param moms_in Input moments
 * @param err_max Tolerance of error in M1 and M2
 * @param iter_max Maximum number of iteration
 * @param conf_local Local configuration space range
 * @param conf_local_ext Local extended configuration space range
 * @param phase_local Local phase-space range
 */
void gkyl_correct_maxwellian_gyrokinetic_fix(gkyl_correct_maxwellian_gyrokinetic *cmax,
  struct gkyl_array *fM, const struct gkyl_array *moms_in, double err_max, int iter_max,
  const struct gkyl_range *conf_local, const struct gkyl_range *conf_local_ext, 
  const struct gkyl_range *phase_local);

/**
 * Delete updater.
 *
 * @param cmax Updater to delete.
 */
void gkyl_correct_maxwellian_gyrokinetic_release(gkyl_correct_maxwellian_gyrokinetic* cmax);
