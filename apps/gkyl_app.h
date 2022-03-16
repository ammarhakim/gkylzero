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
  GKYL_SPECIES_WALL, // perfect reflector
};

// Boundary conditions on fields
enum gkyl_field_bc_type {
  GKYL_FIELD_COPY = 0, // copy BCs
  GKYL_FIELD_PEC_WALL, // perfect electrical conductor (PEC) BCs
};
