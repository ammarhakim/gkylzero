#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>

// BC types in this updater.
enum gkyl_bc_emission_elastic_type {
  GKYL_BS_FURMAN_PIVI = 0,
  GKYL_BS_CAZAUX = 1,
  GKYL_BS_CONSTANT = 2};

struct gkyl_bc_emission_elastic_furman_pivi {
  double mass;
  double charge;
  double P1_inf;
  double P1_hat;
  double E_hat;
  double W;
  double p;
};

struct gkyl_bc_emission_elastic_cazaux {
  double mass;
  double charge;
  double E_f;
  double phi;
};

struct gkyl_bc_emission_elastic_constant{
  double delta;
};

// Object type
typedef struct gkyl_bc_emission_elastic gkyl_bc_emission_elastic;

/**
 * Create a new updater to apply emitting wall spectrum boundary conditions.
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (emission_spectrum gkyl_edge_loc).
 * @param bctype BC spectrum type (see gkyl_bc_emission_spectrum_type).
 * @param gammatype SE yield type (see gkyl_bc_emission_spectrum_type).
 * @param bc_param Parameters used for calculating BC spectrum.
 * @param sey_param Parameters used for calculating SE yield.
 * @param cdim Configuration space dimensions.
 * @param vdim Velocity space dimensions.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_emission_elastic* gkyl_bc_emission_elastic_new(enum gkyl_bc_emission_elastic_type elastic_type,
  void *elastic_param, struct gkyl_array *elastic_yield, int dir, enum gkyl_edge_loc edge,
  int cdim, int vdim, struct gkyl_rect_grid *grid, struct gkyl_range *emit_buff_r, int poly_order,
  struct gkyl_basis *basis, bool use_gpu);

/**
 * @param up BC updater
 * @param f_skin Skin cell distribution
 * @param f_proj Projected spectrum distribution
 * @param f_buff Distribution buffer array
 * @param weight Weighting coefficients
 * @param k Normalization factor
 * @param flux Flux into boundary
 * @param grid Domain grid
 * @param gamma SE yield values on incoming ghost space
 * @param skin_r Incoming skin space range
 * @param ghost_r Incoming ghost space range
 * @param conf_r Configuration space range
 * @param buff_r Buffer array range
 */
void gkyl_bc_emission_elastic_advance(const struct gkyl_bc_emission_elastic *up,
  struct gkyl_range *emit_skin_r, struct gkyl_array *buff_arr, struct gkyl_array *f_skin,
  struct gkyl_array *f_emit, struct gkyl_array *elastic_yield, struct gkyl_basis *basis);

/**
 * Free memory associated with bc_emission_elastic updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_emission_elastic_release(struct gkyl_bc_emission_elastic *up);
