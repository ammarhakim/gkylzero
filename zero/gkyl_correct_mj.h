#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_correct_mj gkyl_correct_mj;

/**
 * Create new updater to correct a Maxwellian to match specified
 * moments.
 *
 * @param grid Grid on which updater lives
 * @param conf_basis Conf space basis functions
 * @param phase_basis Phase space basis functions
 * @param conf_range configuration space range
 * @param conf_range_ext Number of cells in local extended config-
 * @param vel_range velocity space range
 * @param p_over_gamma sr quantitiy: velocity
 * @param gamma sr quantitiy: gamma
 * @param gamma_inv sr quantitiy: 1/gamma
 * @param use_gpu bool for gpu useage
 */
gkyl_correct_mj *gkyl_correct_mj_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_range_ext, const struct gkyl_range *vel_range, 
  const struct gkyl_array *p_over_gamma, const struct gkyl_array *gamma, const struct gkyl_array *gamma_inv, 
  bool use_gpu);

/**
 * Fix the Maxwell-Juttner so that it's moments match desired n_stationary moment.
 *
 * @param cmj MJ correction updater
 * @param fout Distribution function to fix (modified in-place)
 * @param n_stationary Desired lab-frame number density
 * @param vbi specified velocity of f in the stationary frame
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 */
void gkyl_correct_mj_fix_n_stationary(gkyl_correct_mj *cmj, 
  struct gkyl_array *fout,
  const struct gkyl_array *n_stationary, const struct gkyl_array *vbi,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local);

/**
 * Fix the Maxwell-Juttner so that it's moments match desired moments.
 * NOTE: If this algorithm fails, the returns the original distribution function
 * with only the n_stationary moment corrected (i.e. runs: gkyl_correct_mj_fix_n_stationary())
 *
 * @param cmj MJ correction updater
 * @param fout Distribution function to fix (modified in-place)
 * @param n_stationary Desired lab-frame number density
 * @param vbi Desired lab-frame velocity
 * @param T_stationary Desired lab-frame temperature
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 * @param poly_order 
 */
void gkyl_correct_mj_fix(gkyl_correct_mj *cmj,
  struct gkyl_array *fout,
  const struct gkyl_array *n_stationary, const struct gkyl_array *vbi, const struct gkyl_array *T_stationary,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
  int poly_order);

/**
 * Delete updater.
 *
 * @param cmj Updater to delete.
 */
void gkyl_correct_mj_release(gkyl_correct_mj *cmj);
