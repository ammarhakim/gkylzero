#pragma once

// Identifiers for various equation systems
enum gkyl_eqn_type {
  GKYL_EULER, // Euler equations
  GKYL_ISO_EULER, // Isothermal Euler equations
  GKYL_TEN_MOMENT, // Ten-moment (with pressure tensor)
  GKYL_MAXWELL, // Maxwell equations
};

// Identifiers for specific Vlasov equation type
enum gkyl_field_id {
  GKYL_EM, // Vlasov-Maxwell
  GKYL_NO_FIELD, // Neutrals
  GKYL_PHI_ONLY, // Vlasov-Poisson (only phi)
  GKYL_PHI_A, // Vlasov-Poisson (phi + A)
};
