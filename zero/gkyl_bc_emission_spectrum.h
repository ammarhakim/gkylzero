#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>

// BC types in this updater.
enum gkyl_bc_emission_spectrum_norm_type {
  GKYL_SEE_CHUNG_EVERHART = 0,
  GKYL_SEE_GAUSSIAN = 1,
  GKYL_SEE_MAXWELLIAN = 2};

enum gkyl_bc_emission_spectrum_yield_type {
  GKYL_SEE_FURMAN_PIVI = 0,
  GKYL_SEE_SCHOU = 1,
  GKYL_SEE_CONSTANT = 2};

// BC normalization factor structs
struct gkyl_bc_emission_spectrum_norm_gaussian {
  int cdim;
  int vdim;
  double mass;
  double charge;
  double E_0;
  double tau;
};

struct gkyl_bc_emission_spectrum_norm_chung_everhart {
  int cdim;
  int vdim;
  double mass;
  double charge;
  double phi;
};

struct gkyl_bc_emission_spectrum_norm_maxwellian {
  int cdim;
  int vdim;
  double mass;
  double charge;
  double vt;
};

// BC SEY equation structs
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
 * @param norm_type Type of spectrum to use for calculation of normalization factor
 * @param yield_type Type of yield equation for SEY calculation
 * @param norm_param Parameters of normalization factor calculation
 * @param yield_param Parameters of yield calculation
 * @param yield Array of calculated yield values at cell centers
 * @param spectrum Emission spectrum projected onto basis
 * @param dir Direction in which to apply BC
 * @param edge Lower or upper edge at which to apply BC (emission_spectrum gkyl_edge_loc)
 * @param cdim Configuration space dimensions
 * @param vdim Velocity space dimensions
 * @param impact_buff_r Range over the impacting species buffer array
 * @param emit_buff_r Range over the emitting species buffer array
 * @param grid Impacting species boundary grid
 * @param poly_order Polynomial order of basis functions.
 * @param basis Basis functions
 * @param proj_buffer Host array to temporarily store projection of emission spectrum
 * @param use_gpu Boolean to indicate whether to use the GPU
 * @return New updater pointer
 */
struct gkyl_bc_emission_spectrum* gkyl_bc_emission_spectrum_new(enum gkyl_bc_emission_spectrum_norm_type norm_type,
  enum gkyl_bc_emission_spectrum_yield_type yield_type, void *norm_param, void *yield_param,
  struct gkyl_array *yield, struct gkyl_array *spectrum, int dir, enum gkyl_edge_loc edge,
  int cdim, int vdim, struct gkyl_range *impact_buff_r,  struct gkyl_range *emit_buff_r,
  struct gkyl_rect_grid *grid, int poly_order, struct gkyl_basis *basis, struct gkyl_array *proj_buffer, bool use_gpu);

/**
 * @param up BC updater
 * @param impact_buff_r Range over the impacting species buffer array
 * @param impact_cbuff_r Configuration space range over the impacting species buffer array
 * @param emit_buff_r Range over the emitting species buffer array
 * @param bflux Boundary flux df/dt in ghost cells
 * @param f_emit Emitted distribution
 * @param yield Array of calculated yield values at cell centers
 * @param spectrum Emission spectrum projected onto basis
 * @param weight Weighting coefficients
 * @param flux Flux into boundary
 * @param k Normalization factor
 */
void gkyl_bc_emission_spectrum_advance(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_range *impact_buff_r, struct gkyl_range *impact_cbuff_r,
  struct gkyl_range *emit_buff_r, struct gkyl_array *bflux, struct gkyl_array *f_emit,
  struct gkyl_array *yield, struct gkyl_array *spectrum, struct gkyl_array *weight,
  struct gkyl_array *flux, struct gkyl_array *k);

/**
 * Loop over impacting species velocity space and calculate SEY at cell centers
 *
 * @param up BC updater
 * @param yield Array of calculated yield values at cell centers
 * @param grid Impacting species boundary grid
 * @param impact_buff_r Range over the impacting species buffer array
 */
void gkyl_bc_emission_spectrum_sey_calc(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_array *yield, struct gkyl_rect_grid *grid, const struct gkyl_range *impact_buffer_r);

/**
 * Create range over velocity space into wall
 *
 * @param flux_r Output range over impacting velocity space
 * @param dir Direction in which to apply BC
 * @param parent Input range over all of velocity space
 * @param nghost Number of ghost cells
 * @param edge Lower or upper edge at which to apply BC (emission_spectrum gkyl_edge_loc)
 */
void gkyl_bc_emission_flux_ranges(struct gkyl_range *impact_buff_r, int dir,
  const struct gkyl_range *parent, const int *nghost, enum gkyl_edge_loc edge);

/**
 * Free memory associated with bc_emission_spectrum updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_emission_spectrum_release(struct gkyl_bc_emission_spectrum *up);
