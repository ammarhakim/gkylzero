#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_skin_surf_from_ghost gkyl_skin_surf_from_ghost;

/**
 * Create a new updater to copy the ghost node values to the skin node values
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param basis Basis on which coefficients in array are expanded (a device pointer if use_gpu=true).
 * @param skin_r Skin range.
 * @param ghost_r Ghost range.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_skin_surf_from_ghost* gkyl_skin_surf_from_ghost_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_basis *basis, const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  bool use_gpu);

/**
 * Apply the sheath BC with the skin_surf_from_ghost object.
 *
 * @param up BC updater.
 * @param phi Electrostatic potential.
 * @param phi_wall Wall potential.
 * @param distf Distribution function array to apply BC to.
 * @param conf_r Configuration space range (to index phi).
 */
void gkyl_skin_surf_from_ghost_advance(const struct gkyl_skin_surf_from_ghost *up, struct gkyl_array *phi);

/**
 * Free memory associated with skin_surf_from_ghost updater.
 *
 * @param up BC updater.
 */
void gkyl_skin_surf_from_ghost_release(struct gkyl_skin_surf_from_ghost *up);
