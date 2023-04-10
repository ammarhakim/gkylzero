#pragma once

#include <stdbool.h>

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
};

// Boundary conditions on fields
enum gkyl_field_bc_type {
  GKYL_FIELD_COPY = 0, // copy BCs
  GKYL_FIELD_PEC_WALL, // perfect electrical conductor (PEC) BCs
  GKYL_FIELD_RESERVOIR, // Reservoir Maxwell's BCs for heat flux problem
  GKYL_FIELD_WEDGE, // specialized "wedge" BCs for RZ-theta
};
