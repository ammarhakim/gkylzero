#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h> 

// Object type
typedef struct gkyl_lte_moments gkyl_lte_moments;

// input packaged as a struct
struct gkyl_lte_moments_vlasov_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid on which to compute moments
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions
  const struct gkyl_basis *phase_basis; // Phase-space basis functions
  const struct gkyl_range *conf_range; // Configuration-space range
  const struct gkyl_range *conf_range_ext; // Extended configuration-space range (for internal memory allocations)
  const struct gkyl_range *vel_range; // velocity space range
  const struct gkyl_array *p_over_gamma; // SR quantitiy: velocity
  const struct gkyl_array *gamma; // SR quantitiy: gamma = sqrt(1 + p^2)
  const struct gkyl_array *gamma_inv; // SR quantitiy: 1/gamma = 1/sqrt(1 + p^2)
  enum gkyl_model_id model_id; // Enum identifier for model type (e.g., SR, see gkyl_eqn_type.h)
  double mass; // Mass factor 
  bool use_gpu; // bool for gpu useage
};


/**
 * Create new updater to compute the moments for a Maxwellian distribution.
 * Updater is agnostic to relativistic vs. non-relativistic and always returns
 * (n, V_drift, T/m) where n and T/m are the stationary frame density and temperature/mass
 * (the frame moving at velocity V_drift).
 * 
 * @param inp Input parameters
 * @return New updater pointer.
 */
struct gkyl_lte_moments*
gkyl_lte_moments_inew(const struct gkyl_lte_moments_vlasov_inp *inp);

/**
 * Compute the moments of an arbitrary distribution function for an equivalent Maxwellian.
 * Computes (n, V_drift, T/m) where n and T/m are the stationary frame density and temperature/mass
 *
 * @param maxwell_moms Maxwellian moments updater
 * @param phase_local Phase-space range on which to compute moments.
 * @param conf_local Configuration-space range on which to compute moments.
 * @param fin Input distribution function
 * @param moms Output maxwellian moments (n, V_drift, T/m)
 */
void gkyl_lte_moments_advance(struct gkyl_lte_moments *maxwell_moms, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local, 
  const struct gkyl_array *fin, struct gkyl_array *moms);

/**
 * Delete updater.
 *
 * @param cmj Updater to delete.
 */
void gkyl_lte_moments_release(gkyl_lte_moments* maxwell_moms);
