#pragma once

// Identifiers for various equation systems
enum gkyl_eqn_type {
  GKYL_EQN_EULER,     // Euler equations
  GKYL_EQN_SR_EULER,  // SR Euler equations
  GKYL_EQN_ISO_EULER, // Isothermal Euler equations
  GKYL_EQN_COLDFLUID, // Cold fluid equations
  GKYL_EQN_COLDFLUID_SR, // Relativistic Cold fluid equations
  GKYL_EQN_TEN_MOMENT, // Ten-moment (with pressure tensor)
  GKYL_EQN_MAXWELL, // Maxwell equations
  GKYL_EQN_MHD,  // Ideal MHD equations
  GKYL_EQN_BURGERS, // Burgers equations
  GKYL_EQN_ADVECTION, // Scalar advection equation
  GKYL_EQN_GR_EULER, // General relativistic Euler equations with ideal gas equation of state.
  GKYL_EQN_GR_ULTRA_REL_EULER, // General relativistic Euler equations with ultra-relativistic equation of state.
  GKYL_EQN_GR_MAXWELL, // General relativistic Maxwell equations.
  GKYL_EQN_GR_MAXWELL_TETRAD, // General relativistic Maxwell equations in the tetrad basis.
  GKYL_EQN_REACTIVE_EULER, // Reactive Euler equations.
  GKYL_EQN_EULER_MIXTURE, // Euler mixture equations.
  GKYL_EQN_ISO_EULER_MIXTURE, // Isothermal Euler mixture equations.
};

// Identifiers for specific field object types
enum gkyl_field_id {
  GKYL_FIELD_E_B = 0, // Maxwell (E, B). This is default
  GKYL_FIELD_PHI = 1, // Poisson (only phi)
  GKYL_FIELD_PHI_A = 2, // Poisson with static B = curl(A) (phi, A)
  GKYL_FIELD_NULL = 3, // no field is present
};

// Identifiers for subsidary models
// These are used to distinguish things like special relativistic from non-relativistic
enum gkyl_model_id {
  GKYL_MODEL_DEFAULT = 0, // No subsidiary model specified
  GKYL_MODEL_SR = 1,
  GKYL_MODEL_GEN_GEO = 2,
  GKYL_MODEL_CANONICAL_PB = 3,
};

// Identifiers for specific collision object types
enum gkyl_collision_id {
  GKYL_NO_COLLISIONS = 0, // No collisions. This is default
  GKYL_BGK_COLLISIONS, // BGK Collision operator
  GKYL_LBO_COLLISIONS, // LBO Collision operator
  GKYL_FPO_COLLISIONS, // FPO Collision operator
};

// Identifiers for specific source object types
enum gkyl_source_id {
  GKYL_NO_SOURCE = 0, // No source. This is default
  GKYL_FUNC_SOURCE, // Function source
  GKYL_PROJ_SOURCE, // Source given by projection object determined by gkyl_projection_id
  GKYL_BFLUX_SOURCE // Source which scales to boundary fluxes
};

// Identifiers for specific projection object types
enum gkyl_projection_id {
  GKYL_PROJ_FUNC = 0, // Function projection. This is default
  GKYL_PROJ_MAXWELLIAN_PRIM, // Maxwellian projection from primitive moments (n, u, T)
  GKYL_PROJ_MAXWELLIAN_LAB, // Maxwellian projection from lab moments (M0, M1, M2)
  GKYL_PROJ_BIMAXWELLIAN, // Bi-Maxwellian projection
  GKYL_PROJ_VLASOV_LTE, // LTE (Local thermodynamic equilibrium) projection for Vlasov
                        // (Maxwellian for non-relativistic, Maxwell-Juttner for relativistic)
};

// Identifiers for specific radiation object types
enum gkyl_radiation_id {
  GKYL_NO_RADIATION = 0, // No radiation. This is default
  GKYL_GK_RADIATION, // Radiation in gyrokinetic equations.
  GKYL_VM_COMPTON_RADIATION, // Vlasov simple Compton radiation model. 
};

// type of quadrature to use
enum gkyl_quad_type {
  GKYL_GAUSS_QUAD = 0, // Gauss-Legendre quadrature
  GKYL_GAUSS_LOBATTO_QUAD, // Gauss-Lobatto quadrature
};

/** Flags for indicating acting edge of velocity space */
enum gkyl_vel_edge { 
  GKYL_VX_LOWER, GKYL_VY_LOWER, GKYL_VZ_LOWER, 
  GKYL_VX_UPPER, GKYL_VY_UPPER, GKYL_VZ_UPPER 
};
