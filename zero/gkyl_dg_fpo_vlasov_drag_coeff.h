#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

/**
 * Compute the drag coefficient a_i for the drag term of the FPO 
 * as the velocity space gradient of an arbitrary input potential H.
 * Calculates all 3 vector components of a_i. 
 *
 * @param grid Grid (for getting cell spacing)
 * @param cbasis Basis functions in configuration space
 * @param pbasis Basis functions in phase space
 * @param range Range to calculate gradient
 * @param fpo_h Input potential
 * @param fpo_drag_coeff Output array of drag coefficient
*/
void gkyl_calc_fpo_drag_coeff_recovery(const struct gkyl_rect_grid *grid, 
  struct gkyl_basis pbasis, const struct gkyl_range *range, const struct gkyl_range *conf_range, 
  const struct gkyl_array *gamma, const struct gkyl_array *fpo_h, const struct gkyl_array* fpo_dhdv_surf, struct gkyl_array *fpo_drag_coeff, struct gkyl_array *fpo_drag_coeff_surf,
  struct gkyl_array *sgn_drag_coeff_surf, struct gkyl_array* const_sgn_drag_coeff_surf);
