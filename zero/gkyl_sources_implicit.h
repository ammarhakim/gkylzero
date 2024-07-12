#pragma once

#include <gkyl_moment_em_coupling_priv.h>

/**
* Integrate the electromagnetic source terms of a charged multi-fluid equation system within a single cell, using an implicit forcing solver
* (specifically the time-centered Crank-Nicolson/implicit Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
* @param app_accel_s Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param em Array of electromagnetic variables.
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
*/
void implicit_em_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double *app_accel_s[GKYL_MAX_SPECIES], double* em, const double* app_current, const double* ext_em);

/**
* Integrate the momentum source terms of a neutral multi-fluid equation system within a single cell, using an implicit forcing solver (specifically
* the time-centered Crank-Nicolson/implicit Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
* @param app_accel_s Array of acceleration terms to be applied to the fluid equations (for external forces).
*/
void implicit_neut_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES]);

/**
* Integrate the collisional source terms of a multi-fluid equation system within a single cell, using an implicit forcing solver (specifically the
* time-centered Crank-Nicolson/implicit Runge-Kutta method, with a direct matrix inversion).
*
* @param mom_em Moment-EM coupling object.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
*/
void implicit_collision_source_update(const gkyl_moment_em_coupling* mom_em, double dt, double* fluid_s[GKYL_MAX_SPECIES]);

/**
* Integrate the electromagnetic source terms in the multi-fluid equation system within each cell, using an implicit forcing solver (specifically
* the time-centered Crank-Nicolson/implicit Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
* @param app_accel_s Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param p_rhs_s Array of RHS/source terms to be applied to the pressure tensor (for the case of 10-moment gradient-based closure only).
* @param em Array of electromagnetic variables.
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
* @param nT_sources_s Array of number density and temperature source terms.
*/
void implicit_source_coupling_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], const double* p_rhs_s[GKYL_MAX_SPECIES], double* em, const double* app_current,
  const double* ext_em, const double* nT_sources_s[GKYL_MAX_SPECIES]);