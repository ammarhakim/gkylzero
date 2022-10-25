#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_proj_MJ_on_basis gkyl_proj_MJ_on_basis;

/**
 * Create new updater to project MJ on basis functions. Free
 * using gkyl_proj_MJ_on_basis_release method.
 *
 * @param grid Grid object
 * @param conf_basis Conf-space basis functions
 * @param phase_basis Phase-space basis functions
 * @param num_quad Number of quadrature nodes
 * @return New updater pointer.
 */
gkyl_proj_MJ_on_basis* gkyl_proj_MJ_on_basis_new(
  const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  int num_quad);

/**
 * Compute projection of MJ on basis. This method takes
 * lab-frame moments to compute the projection of MJ on basis
 * functions.
 *
 * @param pob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param n Number density moment in the stationary fluid frame
 * @param v velocity of the fluid relative to the lab-frame
 * @param T Temperature in the fluid rest frame

 * @param f_MJ Output MJ
 */
void gkyl_proj_MJ_on_basis_fluid_stationary_frame_mom(const gkyl_proj_MJ_on_basis *pob,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *n_stationary_frame, const struct gkyl_array *vel_stationary_frame, const struct gkyl_array *T_stationary_frame,
  struct gkyl_array *f_MJ);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_MJ_on_basis_release(gkyl_proj_MJ_on_basis* pob);
