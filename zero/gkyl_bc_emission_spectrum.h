#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>

// BC types in this updater.
enum gkyl_bc_emission_spectrum_norm_type {
  GKYL_BC_CHUNG_EVERHART = 0,
  GKYL_BC_GAUSSIAN = 1,
  GKYL_BC_MAXWELLIAN = 2};

enum gkyl_bc_emission_spectrum_yield_type {
  GKYL_BC_FURMAN_PIVI = 0,
  GKYL_BC_SCHOU = 1,
  GKYL_BC_CONSTANT = 2};

struct gkyl_bc_emission_spectrum_norm_gaussian {
  double mass;
  double charge;
  double E_0;
  double tau;
};

struct gkyl_bc_emission_spectrum_norm_chung_everhart {
  double mass;
  double charge;
  double phi;
};

struct gkyl_bc_emission_spectrum_norm_maxwellian {
  double mass;
  double charge;
  double vt;
};

struct gkyl_bc_emission_spectrum_yield_furman_pivi {
  double mass;
  double charge;
  double deltahat_ts;
  double Ehat_ts;
  double t1;
  double t2;
  double t3;
  double t4;
  double s;
};

struct gkyl_bc_emission_spectrum_yield_schou {
  double mass;
  double charge;
  double int_wall;
  double a2;
  double a3;
  double a4;
  double a5;
  double nw;
};

struct gkyl_bc_emission_spectrum_yield_constant{
  double delta;
};

// Object type
typedef struct gkyl_bc_emission_spectrum gkyl_bc_emission_spectrum;

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
struct gkyl_bc_emission_spectrum* gkyl_bc_emission_spectrum_new(enum gkyl_bc_emission_spectrum_norm_type norm_type,
  enum gkyl_bc_emission_spectrum_yield_type yield_type, void *norm_params, void *yield_params,
  struct gkyl_array *yield, struct gkyl_array *spectrum, int dir, enum gkyl_edge_loc edge,
  int cdim, int vdim, struct gkyl_range *impact_skin_r, struct gkyl_range *impact_ghost_r,
  struct gkyl_rect_grid *grid, bool use_gpu);

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
void gkyl_bc_emission_spectrum_advance(const struct gkyl_bc_emission_spectrum *up,
  const struct gkyl_array *f_skin, const struct gkyl_array *f_proj,
  struct gkyl_array *f_buff, struct gkyl_rect_grid *grid);

/**
 * @param up BC updater
 * @param grid Domain grid
 * @param gamma SE yield values on incoming ghost space
 * @param ghost_r Incoming ghost space range
 */
void gkyl_bc_emission_spectrum_sey_calc(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_array *gamma, struct gkyl_rect_grid *grid, const struct gkyl_range *ghost_r);

void gkyl_bc_emission_pos_neg_ranges(struct gkyl_range *pos, struct gkyl_range *neg,
  int dir, const struct gkyl_range *parent, const int *nghost);

/**
 * Free memory associated with bc_emission_spectrum updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_emission_spectrum_release(struct gkyl_bc_emission_spectrum *up);
