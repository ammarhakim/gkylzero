#pragma once

// Identifiers for various equation systems
enum gkyl_eqn_type {
  GKYL_EULER, // Euler equations
  GKYL_ISO_EULER, // Isothermal Euler equations
  GKYL_TEN_MOMENT, // Ten-moment (with pressure tensor)
  GKYL_MAXWELL, // Maxwell equations
};
