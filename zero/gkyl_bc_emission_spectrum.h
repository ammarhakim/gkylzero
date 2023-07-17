#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>

// BC types in this updater.
enum gkyl_bc_emission_spectrum_type {
  GKYL_BC_CHUNG_EVERHART = 0,
  GKYL_BC_GAUSSIAN = 1};

enum gkyl_bc_emission_spectrum_gamma_type {
  GKYL_BC_FURMAN_PIVI = 0};

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
  enum gkyl_bc_emission_spectrum_type bctype,
  enum gkyl_bc_emission_spectrum_gamma_type gammatype,
  double *bc_param, double *sey_param, int cdim, int vdim, bool use_gpu);

void gkyl_bc_emission_spectrum_advance(const struct gkyl_bc_emission_spectrum *up,
  const struct gkyl_array *f_skin, const struct gkyl_array *f_proj, struct gkyl_array *f_buff,
  struct gkyl_array *weight, struct gkyl_array *k,
  const struct gkyl_array *flux, struct gkyl_rect_grid *grid, struct gkyl_array *gamma,
  const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r, const struct gkyl_range *conf_r);

void gkyl_bc_emission_spectrum_sey_calc(const struct gkyl_bc_emission_spectrum *up, struct gkyl_array *gamma, struct gkyl_rect_grid *grid, const struct gkyl_range *ghost_r);

/**
 * Free memory associated with bc_emission_spectrum updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_emission_spectrum_release(struct gkyl_bc_emission_spectrum *up);
