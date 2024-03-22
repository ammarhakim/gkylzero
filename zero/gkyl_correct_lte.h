#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_correct_vlasov_lte gkyl_correct_vlasov_lte;

// Correction status
struct gkyl_correct_vlasov_lte_status {
  bool iter_converged; // true if iterations converged
  int num_iter; // number of iterations for the correction
};  

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
 * Fix the LTE (local thermodynamic equlibrium) distribution function
 * (Maxwellian for non-relativistic/Maxwell-Juttner for relativisticy)
 * so that its 0th moment matches desired stationary-frame density moment.
 *
 * @param c_corr LTE distribution function moment correction updater
 * @param f_lte LTE distribution function to fix (modified in-place)
 * @param moms_target Target stationary-frame moments (n, V_drift, T/m)
 *                    Target density is the 0th component
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 * @return Status of correction
 */
struct gkyl_correct_vlasov_lte_status gkyl_correct_density_moment_vlasov_lte(gkyl_correct_vlasov_lte *c_corr, 
  struct gkyl_array *f_lte, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local);

/**
 * Fix the LTE (local thermodynamic equlibrium) distribution function
 * (Maxwellian for non-relativistic/Maxwell-Juttner for relativisticy)
 * so that *all* its stationary-frame moments (n, V_drift, T/m) match target
moments.
 * NOTE: If this algorithm fails, the returns the original distribution function
 * with only the desired stationary-frame density moment corrected
 * (i.e. runs: gkyl_correct_density_moment_vlasov_lte)
 *
 * @param c_corr LTE distribution function moment correction updater
 * @param f_lte LTE distribution function to fix (modified in-place)
 * @param moms_target Target stationary-frame moments (n, V_drift, T/m)
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 * @return Status of correction
 */
struct gkyl_correct_vlasov_lte_status gkyl_correct_all_moments_vlasov_lte(gkyl_correct_vlasov_lte *c_corr,
  struct gkyl_array *f_lte, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local);

/**
 * Delete updater.
 *
 * @param c_corr Updater to delete.
 */
void gkyl_correct_vlasov_lte_release(gkyl_correct_vlasov_lte *c_corr);
