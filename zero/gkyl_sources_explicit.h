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