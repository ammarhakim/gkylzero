#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_correct_vlasov_lte gkyl_correct_vlasov_lte;

/**
 * Create new updater to correct a Maxwellian/Maxwell-Juttner to match specified
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
 * @param model_id either GKYL_MODEL_SR (MJ) or GKYL_MODEL_DEFAULT (Maxwellian)
 * @param mass mass of the species
 * @param use_gpu bool for gpu useage
 */
gkyl_correct_vlasov_lte *gkyl_correct_vlasov_lte_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_range_ext, const struct gkyl_range *vel_range, 
  const struct gkyl_array *p_over_gamma, const struct gkyl_array *gamma, const struct gkyl_array *gamma_inv, 
  enum gkyl_model_id model_id, double mass, bool use_gpu);

/**
 * Fix the Maxwellian/Maxwell-Juttner so that it's moments match desired n/n_stationary moment.
 *
 * @param c_corr MJ correction updater
 * @param fout Distribution function to fix (modified in-place)
 * @param n_stationary Desired lab-frame number density
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 */
void gkyl_correct_density_moment_vlasov_lte(gkyl_correct_vlasov_lte *c_corr, 
  struct gkyl_array *fout,
  const struct gkyl_array *n_stationary, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local);

/**
 * Fix the Maxwellian/Maxwell-Juttner so that it's moments match desired moments.
 * NOTE: If this algorithm fails, the returns the original distribution function
 * with only the n/n_stationary moment corrected (i.e. runs: gkyl_correct_density_moment_vlasov_lte())
 *
 * @param c_corr MJ correction updater
 * @param fout Distribution function to fix (modified in-place)
 * @param n_stationary Desired lab-frame number density
 * @param vbi Desired lab-frame velocity
 * @param T_stationary Desired lab-frame temperature
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 * @param poly_order 
 */
void gkyl_correct_all_moments_vlasov_lte(gkyl_correct_vlasov_lte *c_corr,
  struct gkyl_array *fout,
  const struct gkyl_array *n_stationary, const struct gkyl_array *vbi, const struct gkyl_array *T_stationary,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
  int poly_order);

/**
 * Delete updater.
 *
 * @param c_corr Updater to delete.
 */
void gkyl_correct_vlasov_lte_release(gkyl_correct_vlasov_lte *c_corr);
