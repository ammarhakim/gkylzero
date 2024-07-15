#pragma once

#include <string.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_sources_implicit.h>
#include <gkyl_sources_explicit.h>

struct gkyl_moment_em_coupling {
  struct gkyl_rect_grid grid; // Grid over which the equations are solved.
  int ndim; // Number of grid dimensions.
  int nfluids; // Number of fluid species.
  struct gkyl_moment_em_coupling_data param[GKYL_MAX_SPECIES]; // Data for each fluid species.

  bool is_charged_species; // Does the fluid species carry a charge?
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.

  bool ramp_app_E; // Use a linear ramp function for initializing external electric fields.
  double t_ramp_E; // Ramp-up time for the linear ramp function for initializing external electric fields.
  bool ramp_app_curr; // Use a linear ramp function for initializing applied currents.
  double t_ramp_curr; // Ramp-up time for the linear ramp function for initializing applied currents.

  bool has_collision; // Run with collisions switched on.
  bool use_rel; // Assume special relativistic fluid species.

  // Matrix of scaling factors for collision frequencies. Should be symmetric (i.e. nu_base_sr = nu_base_rs).
  // These are defined such that nu_sr = nu_base_sr / rho_s, and nu_rs = nu_base_rs / rho_r.
  double nu_base[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];

  bool use_explicit_em_coupling; // Use the explicit source-solver for handling moment-EM coupling (not operational yet).

  bool has_nT_sources; // Run with number density and temperature sources.
  bool has_frictional_sources; // Run with frictional sources.
};