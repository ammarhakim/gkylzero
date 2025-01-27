#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_skin_surf_from_ghost gkyl_skin_surf_from_ghost;

/**
 * Create a new updater to copy ghost cell values to skin cells (boundary cells).
 * Specifically, the updater enforces that the skin cell value at the interface (e.g., 
 * upper edge) equals the ghost cell value, achieving continuity across the interface 
 * while maintaining the other edge nodal value.
 * Mathematically, the updater alters the DG coefficient of the skin cell so that
 *  f(x_s) = f_g(x_s) ∀ x_s ∈ ∂Ω
 * where f is the manipulated field, f_g is the ghost field, and ∂Ω is the 
 * boundary of the domain where the matching is enforced.
 * 
 * Example for GKYL_UPPER_EDGE in 1D:
 *   vls            vus  vlg                vug
 *    |___skin cell___|   |____ghost cell____|
 *   legend: f(x_s) = vus the upper nodal skin cell value, 
 *           f_g(x_s) = vlg the lower nodal ghost cell value.
 *  In this situation, ∂Ω is only one point and the updater will alter the DG coefficient so that 
 *    vus = vlg
 *  keeping vls unchanged.
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
 * at that same node. The other nodal value of the skin cell remains unchanged (see new routine description above).
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