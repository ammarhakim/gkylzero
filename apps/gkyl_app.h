#pragma once

#include <stdbool.h>
#include <gkyl_array_rio.h>

// Update status
struct gkyl_update_status {
  bool success; // status of update
  double dt_actual; // actual time-step taken
  double dt_suggested; // suggested stable time-step
};

// Boundary conditions on particles
enum gkyl_species_bc_type {
  GKYL_SPECIES_COPY = 0, // copy BCs
  GKYL_SPECIES_REFLECT, // perfect reflector
  GKYL_SPECIES_ABSORB, // Absorbing BCs
  GKYL_SPECIES_NO_SLIP, // no-slip boundary conditions
  GKYL_SPECIES_WEDGE, // specialized "wedge" BCs for RZ-theta
  GKYL_SPECIES_FUNC, // Function boundary conditions
  GKYL_SPECIES_FIXED_FUNC, // Fixed function, time-independent, boundary conditions
  GKYL_SPECIES_ZERO_FLUX, // Zero flux BCs; must be applied on both lower and upper BC
  GKYL_SPECIES_GK_SHEATH, // Gyrokinetic sheath BCs
  GKYL_SPECIES_RECYCLE, // Recycling BCs
  GKYL_SPECIES_GK_IWL, // Gyrokinetic inner wall limited.
};

// Boundary conditions on fields
enum gkyl_field_bc_type {
  GKYL_FIELD_COPY = 0, // copy BCs
  GKYL_FIELD_PEC_WALL, // Maxwell's perfect electrical conductor (zero normal B and zero tangent E)
  GKYL_FIELD_SYM_WALL, // Maxwell's symmetry BC (zero normal E and zero tangent B)
  GKYL_FIELD_RESERVOIR, // Reservoir Maxwell's BCs for heat flux problem
  GKYL_FIELD_WEDGE, // specialized "wedge" BCs for RZ-theta
  GKYL_FIELD_FUNC, // Function boundary conditions
};
