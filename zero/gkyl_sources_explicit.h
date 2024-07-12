#pragma once

#include <gkyl_moment_em_coupling_priv.h>

/**
* Integrate the number density and temperature source terms for a single fluid species within a single cell, using an explicit forcing solver
* (specifically the simple forward-Euler method).
*
* @param mass Mass of the fluid species.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables for a single fluid species (before source update).
* @param fluid_new Array of new fluid variables for a single fluid species (after source update).
* @param nT_sources Array of number density and temperature source terms for a single fluid species.
*/
void explicit_nT_source_update_euler(const double mass, const double dt, double* fluid_old, double* fluid_new, const double* nT_sources);

/**
* Integrate the number density and temperature source terms in the multi-fluid equation system within a single cell, using an explicit forcing
* solver (specifically the simple forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
* @param nT_sources_s Array of number density and temperature source terms.
*/
void explicit_nT_source_update(const gkyl_moment_em_coupling* mom_em, const double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* nT_sources_s[GKYL_MAX_SPECIES]);

/**
* Integrate the electromagnetic source terms in the multi-fluid equation system within each cell, using an explicit forcing solver (specifically
* a simple forward-Euler method combined with a Higuera-Cary update).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation tine.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
* @param app_accel_s Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param em Array of electromagnetic variables.
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param app_current1 Array of stage-1 current terms to be applied to the fluid equations (for stage-1 of external current driving).
* @param app_current2 Array of stage-2 current terms to be applied to the fluid equations (for stage-2 of external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
* @param nstrang Indictator of which step in the Strang splitting we are currently considering.
*/
void explicit_source_coupling_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], double* em, const double* app_current, const double* app_current1, const double* app_current2,
  const double* ext_em, int nstrang);