#pragma once

// Object type
typedef struct gkyl_gkgeom gkyl_gkgeom;

// Inputs to create a new GK geometry creation object
struct gkyl_gkgeom_inp {
  const struct gkyl_rect_grid *rzgrid; // RZ grid on which psi(R,Z) is defined
  const struct gkyl_basis *rzbasis; // basis functions for R,Z grid

  const struct gkyl_array *psiRZ; // psi(R,Z) DG representation
};

/**
 * Create new updater to compute the geometry (mapc2p) needed in GK
 * simulations.
 *
 * @param inp Input parameters
 * @param New GK geometry updater
 */
gkyl_gkgeom *gkyl_gkgeom_new(const struct gkyl_gkgeom_inp *inp);

/**
 * Delete updater.
 *
 * @param geo Geometry obejct to delete
 */
void gkyl_gkgeom_release(gkyl_gkgeom *geo);
