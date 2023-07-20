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
 * @param conf_basis Conf space basis functions
 * @param phase_basis Phase space basis functions
 * @param conf_local Local configuration space range
 * @param conf_local_ext Local extended configuration space range
 * @param mass Mass of the species
 * @param use_gpu Bool to determine if on GPU
 */
gkyl_correct_maxwellian *gkyl_correct_maxwellian_new(
  const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *conf_local, const struct gkyl_range *conf_local_ext,
  const struct gkyl_array *bmag, double mass, bool use_gpu);

/**
 * Fix the Maxwellian so that it's moments match desired moments.
 *
 * @param cmax Maxwellian-fix updater
 * @param fM Distribution function to fix (modified in-place)
 * @param m0_corr Desired number density
 * @param m1_corr Desired moment density
 * @param m2_corr Desired kinetic energy density
 * @param conf_local Local configuration space range
 * @param conf_local_ext Local extended configuration space range
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 * @param poly_order Polynomial order of basis functions
 * @param use_gpu Bool to determine if on gpu 
 */
void gkyl_correct_maxwellian_gyrokinetic(gkyl_correct_maxwellian *cmax,
  struct gkyl_array *fM, 
  const struct gkyl_array *m0_corr, const struct gkyl_array *m1_corr, const struct gkyl_array *m2_corr,
  const struct gkyl_array *jacob_tot, const struct gkyl_array *bmag, double mass,
  double err_max, int iter_max,
  const struct gkyl_range *conf_local, const struct gkyl_range *conf_local_ext, 
  const struct gkyl_range *phase_local, 
  int poly_order, bool use_gpu);

/**
 * Delete updater.
 *
 * @param cmax Updater to delete.
 */
void gkyl_correct_maxwellian_release(gkyl_correct_maxwellian* cmax);
