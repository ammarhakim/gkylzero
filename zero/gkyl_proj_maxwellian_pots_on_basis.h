#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_proj_maxwellian_pots_on_basis gkyl_proj_maxwellian_pots_on_basis;

/**
 * Create a new updater to project H and G potentials of a Mxwellian
 * on basis functions. Free after use with gkyl_proj_maxwellian_pots_on_basis_release.
 *
 * @param grid Grid object
 * @param conf_basis Configuration space basis functions
 * @param phase_basis Phase space basis functions
 * @param num_quad Number of quadrature nodes
 * @return New updater pointer
*/
gkyl_proj_maxwellian_pots_on_basis* gkyl_proj_maxwellian_pots_on_basis_new(
    const struct gkyl_rect_grid *grid,
    const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
    int num_quad);

/**
 * @param mpob Project on basis updater to run
 * @param phase_range Phase space range
 * @param conf_range Configuration space range
 * @param moms Velocity momends (m0, m1i, m2)
 * @param max_pots Output potentials
*/
void gkyl_proj_maxwellian_pots_on_basis_lab_mom(const gkyl_proj_maxwellian_pots_on_basis *up,
    const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
    const struct gkyl_array *moms, const struct gkyl_array *gamma, const double mass, struct gkyl_array *max_pots);

/**
 * Delete updater.
 *
 *@param mpob Updater to delete.
*/
void gkyl_proj_maxwellian_pots_on_basis_release(gkyl_proj_maxwellian_pots_on_basis *mpob);
