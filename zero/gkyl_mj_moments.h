#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h> 

// Object type
typedef struct gkyl_mj_moments gkyl_mj_moments;

/**
 * Create new updater to compute the moments.
 *
 * @param grid Grid on which updater lives
 * @param conf_basis Conf space basis functions
 * @param phase_basis Phase space basis functions
 * @param conf_range configuration space range
 * @param vel_range velocity space range
 * @param conf_local_cells Number of cells in local config-space
 * @param conf_local_ext_cells Number of cells in local extended config-space
 * @param p_over_gamma sr quantitiy: velocity
 * @param gamma sr quantitiy: gamma
 * @param gamma_inv sr quantitiy: 1/gamma
 * @param use_gpu bool for gpu useage
 */
gkyl_mj_moments *gkyl_mj_moments_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range,
  long conf_local_ncells, long conf_local_ext_ncells, 
  const struct gkyl_array *p_over_gamma, const struct gkyl_array *gamma, const struct gkyl_array *gamma_inv, 
  bool use_gpu);

/**
 * Update the n_stationary, vb, and T_stationary moments (n, vb, T) moments of an arbitary sr
 * distribution so they are ready as inputs to the mj routine
 *
 * @param cmj Maxwell-Juttner correction updater
 * @param fin Distribution function
 * @param n_stationary Desired number density
 * @param vb Desired velocity
 * @param T_stationary Desired Temperature
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 */
void gkyl_mj_moments_advance(gkyl_mj_moments *cmj, 
  const struct gkyl_array *fin,
  struct gkyl_array *n_stationary, struct gkyl_array *vb, struct gkyl_array *T_stationary,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local);

/**
 * Delete updater.
 *
 * @param cmj Updater to delete.
 */
void gkyl_mj_moments_release(gkyl_mj_moments* cmj);
