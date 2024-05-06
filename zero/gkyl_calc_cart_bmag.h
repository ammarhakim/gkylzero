#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_calc_cart_bmag gkyl_calc_cart_bmag;

/**
 * Create new updater to compute the bmag on the compuational grid 
 *
 * @param basis (physical RZ basis)
 * @param grid physical RZ grid
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_cart_bmag*
gkyl_calc_cart_bmag_new(const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid, const struct gkyl_range *local, const struct gkyl_range *local_ext, bool use_gpu);


/**
 * Advance calc_cart_bmag (compute bmag given dg fields Psi Psi/R and Psi/R^2 on the RZ grid 
 *
 * @param up calc_cart_bmag updater object.
 * @param local, local_ext physical RZ range and extended range.
 * @param psidg, psibyrdg, psibyr2dg: DG Psi(R,Z), Psi(R,Z)/R, Psi(R,Z)/R^2 on the RZ grid
 */
void
gkyl_calc_cart_bmag_advance(const gkyl_calc_cart_bmag *up, const struct gkyl_array *psidg, const struct gkyl_array *psibyrdg, const struct gkyl_array *psibyr2dg);


/**
 * Get the cartesian components of the magnetic field at a specificif cartesian xyz coordinate
 * @param up gkyl_calc_cart_bmag updater
 * @param xcart XYZ coordinates
 * @param bcart B_x, B_Y, B_Z to be populated
 */
void
gkyl_eval_cart_bmag(const gkyl_calc_cart_bmag *up, double xcart[3], double bcart[3]);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_cart_bmag_release(gkyl_calc_cart_bmag* up);
