#pragma once

#include <gkyl_fpo_vlasov_coeffs_correct_priv.h>

/**
 * Create new updater to compute the necessary corrections to drag and diffusion 
 * coefficients to enforce momentum and energy conservation in the FPO.
 *
*/ 
struct gkyl_fpo_coeffs_correct*
gkyl_fpo_coeffs_correct_new(bool use_gpu, const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_rng);

/**
 * Solve linear system to calculate corrections to drag and diffusion coefficients,
 * and accumulate them onto the existing drag and diffusion coefficient arrays.
 *
*/
void gkyl_fpo_coeffs_correct_advance(struct gkyl_fpo_coeffs_correct *up,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  const struct gkyl_array *fpo_moms, const struct gkyl_array *boundary_corrections,
  const struct gkyl_array *moms, struct gkyl_array *drag_diff_coeff_corrs,
  struct gkyl_array *drag_coeff, struct gkyl_array *drag_coeff_surf,
  struct gkyl_array *diff_coeff, struct gkyl_array *diff_coeff_surf);

