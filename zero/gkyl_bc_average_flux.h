#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>

// Object type
typedef struct gkyl_bc_average_flux gkyl_bc_average_flux;

/**
 * Create a new updater to apply conducting sheath BCs in gyrokinetics.
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param basis Basis on which coefficients in array are expanded (a device pointer if use_gpu=true).
 * @param skin_r Skin range.
 * @param ghost_r Ghost range.
 * @param vel_map Velocity space mapping object.
 * @param cdim Configuration space dimensions.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_average_flux* gkyl_bc_average_flux_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_basis *basis, const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  const struct gkyl_velocity_map *vel_map, int cdim, bool use_gpu);

/**
 * Apply the sheath BC with the bc_average_flux object.
 *
 * @param up BC updater.
 * @param buff Buffer to copy from. This stores the data at R=0
 * @param distf Distribution function to update.
 */
void gkyl_bc_average_flux_advance(const struct gkyl_bc_average_flux *up, const struct gkyl_array *buff,
   struct gkyl_array *distf);

/**
 * Free memory associated with bc_average_flux updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_average_flux_release(struct gkyl_bc_average_flux *up);
