#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fv_proj.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

struct gkyl_moment_em_coupling_data {
  enum gkyl_eqn_type type; // Equation type.

  double charge; // Species charge.
  double mass; // Species mass.

  double k0; // Closure parameter (for 10-moment equations onl; defaults to 0.0).
};

struct gkyl_moment_em_coupling_inp {
  const struct gkyl_rect_grid *grid; // Grid over which the equations are solved.
  int nfluids; // Number of fluid species.
  struct gkyl_moment_em_coupling_data param[GKYL_MAX_SPECIES]; // Data for each fluid species.
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.

  double t_ramp_E; // Ramp-up time for the linear ramp function for initializing external electric fields.
  double t_ramp_curr; // Ramp-up time for the linear ramp function for initializing applied currents.

  bool has_collision; // Run with collisions switched on.
  bool use_rel; // Assume special relativistic fluid species.
  
  // Matrix of scaling factors for collision frequencies. Should be symmetric (i.e. nu_base_sr = nu_base_rs).
  // These are defined such that nu_sr = nu_base_sr / rho_s, and nu_rs = nu_base_rs / rho_r.
  double nu_base[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];

  bool use_explicit_em_coupling; // Use the explicit source-solver for handling moment-EM coupling (not operational yet).

  bool has_nT_sources; // Run with number density and temperature sources.
};

// Moment-EM coupling object.
typedef struct gkyl_moment_em_coupling gkyl_moment_em_coupling;

/**
* Create a new moment-EM coupling object, used for integrating the electromagnetic source terms that appear in the multi-fluid equations.
*
* @param inp Input parameters for the moment-EM coupling object.
* @return Moment-EM coupling object.
*/
gkyl_moment_em_coupling* gkyl_moment_em_coupling_new(struct gkyl_moment_em_coupling_inp inp);

/**
* Integrate the electromagnetic source terms in the multi-fluid equation system using an implicit forcing solver (specifically the time-centered
* Crank-Nicolson/implicit Runge-Kutta method). The gkyl_range object to be updated should be a (non-strict) subrange of the range over which the
* fluid variable array is defined (i.e. either the fluid variable array gkyl_range object itself, or a gkyl_range object initialized using the
* gkyl_sub_range_init method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param update_range Range object over which to integrate the electromagnetic sources using an implicit time-centered method.
* @param fluid Array of fluid variables (array size = nfluids).
* @param app_accel Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param p_rhs Array of RHS/source terms to be applied to the pressure tensor (for the case of 10-moment gradient-based closure only).
* @param em Array of electromagnetic variables.
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
* @param nT_sources Array of number density and temperature source terms.
*/
void gkyl_moment_em_coupling_implicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], const struct gkyl_array *p_rhs[GKYL_MAX_SPECIES],
  struct gkyl_array* em, const struct gkyl_array* app_current, const struct gkyl_array* ext_em, const struct gkyl_array* nT_sources[GKYL_MAX_SPECIES]);

/**
* Integrate the electromagnetic source terms in the multi-fluid equation system using an explicit forcing solver (specifically either the strong
* stability-preserving Runge-Kutta method or the simple forward Euler method). The gkyl_range object to be updated should be a (non-strict) subrange
* of the range over which the fluid variable array is defined (i.e. either the fluid variable array gkyl_range object itself, or a gkyl_range object
* initialized using the gkyl_sub_range_init method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param update_range Range object over which to integrate the electromagnetic sources using an explicit time integration method.
* @param fluid Array of fluid variables (array size = nfluids).
* @param app_accel Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param p_rhs Array of RHS/source terms to be applied to the pressure tensor (for the case of 10-moment gradient-based closure only).
* @param em Array of electromagnetic variables.
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param app_current1 Array of stage-1 current terms to be applied to the fluid equations (for stage-1 of external current driving).
* @param app_current2 Array of stage-2 current terms to be applied to the fluid equations (for stage-2 of external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
* @param nT_sources Array of number density and temperature source terms.
*/
void gkyl_moment_em_coupling_explicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], const struct gkyl_array* p_rhs[GKYL_MAX_SPECIES],
  struct gkyl_array* em, const struct gkyl_array *app_current, const struct gkyl_array* app_current1, const struct gkyl_array* app_current2,
  const struct gkyl_array* ext_em, const struct gkyl_array* nT_sources[GKYL_MAX_SPECIES], gkyl_fv_proj* proj_app_curr, int nstrang);

/**
* Delete moment-EM coupling object.
*
* @param mom_em Moment-EM coupling object to delete.
*/
void gkyl_moment_em_coupling_release(gkyl_moment_em_coupling* mom_em);