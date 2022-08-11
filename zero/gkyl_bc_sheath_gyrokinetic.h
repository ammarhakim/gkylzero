#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_bc_sheath_gyrokinetic gkyl_bc_sheath_gyrokinetic;

/**
 * Create a new updater to apply conducting sheath BCs in gyrokinetics.
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param local_range_ext Local extended range.
 * @param num_ghosts Number of ghosts in each dimension.
 * @param basis Basis on which coefficients in array are expanded.
 * @param grid cartesian grid dynamic field is defined on.
 * @param cdim Configuration space dimensions.
 * @param q2Dm charge-to-mass ratio times 2.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_sheath_gyrokinetic* gkyl_bc_sheath_gyrokinetic_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *local_range_ext, const int *num_ghosts, const struct gkyl_basis *basis,
  const struct gkyl_rect_grid *grid, int cdim, double q2Dm, bool use_gpu);

/**
 * Create new updater to apply basic BCs to a field
 * in a gkyl_array. Basic BCs are those in which the
 * ghost cell depends solely on the skin cell next to it
 * via a function of type array_copy_func_t (e.g. absorb, reflect).
 *
 * @param up BC updater.
 * @param phi Electrostatic potential.
 * @param phi_wall Wall potential.
 * @param distf Distribution function array to apply BC to.
 */
void gkyl_bc_sheath_gyrokinetic_advance(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_array *distf);

/**
 * Free memory associated with bc_sheath_gyrokinetic updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_sheath_gyrokinetic_release(struct gkyl_bc_sheath_gyrokinetic *up);
