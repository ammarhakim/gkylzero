#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_skin_surf_from_ghost gkyl_skin_surf_from_ghost;

/**
 * Create a new updater to copy ghost cell values to skin cells (boundary cells).
 *
 * @param dir Direction along which to apply the skin-ghost boundary copy.
 * @param edge Lower or upper edge where the boundary condition is applied (see gkyl_edge_loc).
 * @param basis Basis on which the array coefficients are expanded.
 * @param skin_r Range representing the skin (boundary) region.
 * @param ghost_r Range representing the ghost (outer boundary) region.
 * @param use_gpu Boolean flag to indicate whether GPU computation should be used.
 * @return Pointer to the newly created updater.
 */
struct gkyl_skin_surf_from_ghost* gkyl_skin_surf_from_ghost_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_basis basis, const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  bool use_gpu);

/**
 * Enforce that the value of the skin cell at the node facing the ghost cell is equal to the ghost value
 * at that same node. The other nodal value of the skin cell remains unchanged.
 *
 * @param up Pointer to the boundary condition updater.
 * @param field Array representing the field values to update (currently works only in configuration space).
 */
void gkyl_skin_surf_from_ghost_advance(const struct gkyl_skin_surf_from_ghost *up, struct gkyl_array *field);

/**
 * Free the memory associated with the skin_surf_from_ghost updater.
 *
 * @param up Pointer to the boundary condition updater to be released.
 */
void gkyl_skin_surf_from_ghost_release(struct gkyl_skin_surf_from_ghost *up);