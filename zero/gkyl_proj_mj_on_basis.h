#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h> 

// Object type
typedef struct gkyl_proj_mj_on_basis gkyl_proj_mj_on_basis;

/**
 * Create new updater to project the relativistic LTE distribution
 * function on basis functions: the Maxwell-Juttner distribution function. 
 * Free using gkyl_proj_mj_on_basis_release method.
 *
 * @param grid Phase-space grid object
 * @param conf_basis Configuration-space basis functions
 * @param phase_basis Phase-space basis functions
 * @param num_quad Number of quadrature nodes
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_proj_mj_on_basis* gkyl_proj_mj_on_basis_new(
  const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  int num_quad, bool use_gpu);

/**
 * Compute projection of relativistic LTE distribution on basis. 
 * This method takes the LTE moments as a single array moms = (n, V_drift, T/m)
 * to compute the projection of the Maxwell-Juttner distribution function on basis functions.
 * Note: these moments are *defined in the stationary frame moving at V_drift*.
 *
 * @param proj_on_basis Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Configuration-space range
 * @param moms_lte LTE moments for computing Maxwell-Juttner distribution (n, V_drift, T/m)
 *                 Note: LTE moments are defined in stationary frame (frame moving at V_drift)
 * @param f_mj Output Maxwell-Juttner distribution
 */
void gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(const gkyl_proj_mj_on_basis *proj_on_basis,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms_lte, struct gkyl_array *f_mj);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_mj_on_basis_release(gkyl_proj_mj_on_basis* pob);
