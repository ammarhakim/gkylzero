#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>

// BC types in this updater.
enum gkyl_bc_emission_spectrum_type {
  GKYL_BC_CHUNG_EVERHART = 0,
  GKYL_BC_GAUSSIAN = 1};

// Object type
typedef struct gkyl_bc_emission_spectrum gkyl_bc_emission_spectrum;

/**
 * Create a new updater to apply conducting sheath BCs in gyrokinetics.
 *
 * @param grid cartesian grid dynamic field is defined on.
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (emission_spectrum gkyl_edge_loc).
 * @param local_conf_range_ext Local extended configuration range.
 * @param local_range_ext Local extended phase range.
 * @param num_ghosts Number of ghosts in each dimension.
 * @param bctype BC type (see gkyl_bc_emission_spectrum_type).
 * @param cbasis Configuration basis on which coefficients in array are expanded
 * @param basis Basis on which coefficients in array are expanded (a device pointer if use_gpu=true).
 * @param cdim Configuration space dimensions.
 * @param vdim Velocity space dimensions.
 * @param bc_param Parameters used for calculating BCs.
 * @param gain Array of secondary electron gain values at the cell centers.
 * @param elastic Array of elastic backscattering gain values at the cell centers.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_emission_spectrum* gkyl_bc_emission_spectrum_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *ghost_r,
  enum gkyl_bc_emission_spectrum_type bctype, const struct gkyl_basis *cbasis, const struct gkyl_basis *basis,
  int cdim, int vdim, bool use_gpu);

double gkyl_bc_emission_spectrum_advance_cross(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_array *f_self, struct gkyl_array *f_other, struct gkyl_array *f_proj,
  double *bc_param, double flux, struct gkyl_rect_grid *grid, double *gain, const struct gkyl_range *other_r);

/**
 * Advance boundary conditions. Fill buffer array based on boundary conditions and copy
 * contents to ghost cells of input f_arr
 *
 * @param up BC updater.
 * @param buff_arr Buffer array, big enough for ghost cells at this boundary.
 * @param f_arr Field array to apply BC to.
 */
void gkyl_bc_emission_spectrum_advance(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_array *buff_arr, struct gkyl_array *f_arr, struct gkyl_array *f_proj, double k);

/**
 * Free memory associated with bc_emission_spectrum updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_emission_spectrum_release(struct gkyl_bc_emission_spectrum *up);
