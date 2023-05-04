#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

/**
 * Compute the diffusion tensor D_ij for the diffusion term of the FPO
 * from an input potential G. For now, only calculates the diagonal terms.
 * 
 * @param grid Grid (for getting cell spacing)
 * @param cbasis Basis functions in configuration space
 * @param pbasis Basis functions in phase space
 * @param range Range to calculate gradient
 * @param fpo_g Input potential
 * @param fpo_diff_coeff Output array of drag coefficient

*/
void gkyl_calc_fpo_diff_coeff_recovery(const struct gkyl_rect_grid *grid, 
    struct gkyl_basis pbasis, const struct gkyl_range *range,
    const struct gkyl_array *fpo_g, struct gkyl_array *fpo_diff_coeff);
