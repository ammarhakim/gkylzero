#pragma once

#include <gkyl_gk_block_geom.h>
#include <gkyl_gyrokinetic.h>

typedef struct gkyl_gyrokinetic_multib_app gkyl_gyrokinetic_multib_app;

// Boundary conditions on fields and particles
enum gkyl_gyrokinetic_bc_type {
  GKYL_BC_GK_SKIP, // Do not apply any BCs 
  GKYL_BC_GK_SPECIES_REFLECT, // perfect reflector
  GKYL_BC_GK_SPECIES_ABSORB, // Absorbing BCs
  GKYL_BC_GK_SPECIES_FUNC, // Function boundary conditions
  GKYL_BC_GK_SPECIES_FIXED_FUNC, // Fixed function, time-independent, boundary conditions
  GKYL_BC_GK_SPECIES_EMISSION, // Emission spectrum BCs
  GKYL_BC_GK_SPECIES_ZERO_FLUX, // Zero flux BCs; must be applied on both lower and upper BC
  GKYL_BC_GK_SPECIES_GK_SHEATH, // Gyrokinetic sheath BCs
  GKYL_BC_GK_SPECIES_RECYCLE, // Recycling BCs
  GKYL_BC_GK_SPECIES_GK_IWL, // Gyrokinetic inner wall limited.

  GKYL_BC_GK_FIELD_DIRICHLET, // Dirichlet boundary conditions
  GKYL_BC_GK_FIELD_NEUMANN, // Nemann boundary conditions
};

// BC for blocks
struct gkyl_gyrokinetic_block_physical_bcs {
  int bidx; // block index
  int dir;  // direction in which BC is specified
  enum gkyl_edge_loc edge; // which edge this BC is for
  enum gkyl_gyrokinetic_bc_type bc_type; // BC type flag.
  void (*aux_profile)(double t, const double *xn, double *fout, void *ctx); // Auxiliary function (e.g. wall potential).
  void *aux_ctx; // Context for aux_profile.
  double aux_parameter; // Parameter for aux_profile (maybe redundant).
  struct gkyl_gyrokinetic_projection projection; // Projection object input (e.g. for FIXED_FUNC).
};

// Species input per block
struct gkyl_gyrokinetic_multib_species_pb {
  int block_id; // block ID

  // Initial conditions using projection routine.
  struct gkyl_gyrokinetic_projection projection;

  double polarization_density;

  // Source to include.
  struct gkyl_gyrokinetic_source source;
};

// Species input
struct gkyl_gyrokinetic_multib_species {
  char name[128]; // Species name.

  enum gkyl_gkmodel_id gkmodel_id;
  double charge, mass; // Charge and mass.

  double lower[3], upper[3]; // Lower, upper bounds of velocity-space.
  int cells[3]; // Velocity-space cells.
  struct gkyl_mapc2p_inp mapc2p; // Velocity mapping.


  bool no_by; // Boolean for whether we are using specialized GK kernels with no b_y.
              // These more computationally efficient kernels are for slab or mirror 
              // calculations where there is no toroidal field. 

  int num_diag_moments; // Number of diagnostic moments.
  enum gkyl_distribution_moments diag_moments[12]; // List of diagnostic moments.
  int num_integrated_diag_moments; // Number of integrated diagnostic moments.
  enum gkyl_distribution_moments integrated_diag_moments[12]; // List of integrated diagnostic moments.
  bool time_rate_diagnostics; // Whether to ouput df/dt diagnostics.

  struct gkyl_phase_diagnostics_inp boundary_flux_diagnostics;

  // Collisions to include.
  struct gkyl_gyrokinetic_collisions collisions;

  // Diffusion coupling to include.
  struct gkyl_gyrokinetic_diffusion diffusion;

  // Radiation to include.
  struct gkyl_gyrokinetic_radiation radiation;

  // Reactions between plasma species to include.
  struct gkyl_gyrokinetic_react react;
  // Reactions with neutral species to include.
  struct gkyl_gyrokinetic_react react_neut;


  bool duplicate_across_blocks; // Set to true if all blocks are identical.
  // Species inputs per-block: only one is needed if duplicate_across_blocks = true.
  const struct gkyl_gyrokinetic_multib_species_pb *blocks;

  // Physical boundary conditions.
  int num_physical_bcs;
  const struct gkyl_gyrokinetic_block_physical_bcs *bcs;
};

// Neutral species input per block
struct gkyl_gyrokinetic_multib_neut_species_pb {
  int block_id; // block ID

  // Initial conditions using projection routine.
  struct gkyl_gyrokinetic_projection projection;

  // Source to include.
  struct gkyl_gyrokinetic_source source;
};

// Neutral species input
struct gkyl_gyrokinetic_multib_neut_species {
  char name[128]; // Species name.

  double mass; // Mass.
  double lower[3], upper[3]; // Lower, upper bounds of velocity-space.
  int cells[3]; // Velocity-space cells.

  struct gkyl_mapc2p_inp mapc2p;

  bool is_static; // Set to true if neutral species does not change in time.

  int num_diag_moments; // Number of diagnostic moments.
  enum gkyl_distribution_moments diag_moments[12]; // List of diagnostic moments.

  // Reactions with plasma species to include.
  struct gkyl_gyrokinetic_react react_neut;

  bool duplicate_across_blocks; // set to true if all blocks are identical  
  // species inputs per-block: only one is needed if duplicate_across_blocks = true
  const struct gkyl_gyrokinetic_multib_neut_species_pb *blocks;

  // Physical boundary conditions
  int num_physical_bcs;
  const struct gkyl_gyrokinetic_block_physical_bcs *bcs;
};

// Field input per block 
struct gkyl_gyrokinetic_multib_field_pb {
  int block_id; // block ID

  double polarization_bmag; // B  to use in linearized polarization weight
  double kperpSq; // kperp^2 parameter for 1D field equations

  bool time_rate_diagnostics; // Writes the time rate of change of field energy.

  void *phi_wall_lo_ctx; // context for biased wall potential on lower wall
  // pointer to biased wall potential on lower wall function
  void (*phi_wall_lo)(double t, const double *xn, double *phi_wall_lo_out, void *ctx);
  bool phi_wall_lo_evolve; // set to true if biased wall potential on lower wall function is time dependent  

  void *phi_wall_up_ctx; // context for biased wall potential on upper wall
  // pointer to biased wall potential on upper wall function
  void (*phi_wall_up)(double t, const double *xn, double *phi_wall_up_out, void *ctx);
  bool phi_wall_up_evolve; // set to true if biased wall potential on upper wall function is time dependent  
};

// Field input 
struct gkyl_gyrokinetic_multib_field {
  enum gkyl_gkfield_id gkfield_id;
  double kperpSq; // kperp^2 parameter for 1D field equations

  // parameters for adiabatic electrons simulations
  double electron_mass, electron_charge, electron_density, electron_temp;

  //struct gkyl_poisson_bc poisson_bcs;

  bool duplicate_across_blocks; // set to true if all blocks are identical  
  // field inputs per-block: only one is needed if duplicate_across_blocks = true
  const struct gkyl_gyrokinetic_multib_field_pb *blocks;

  // Physical boundary conditions
  int num_physical_bcs;
  const struct gkyl_gyrokinetic_block_physical_bcs *bcs;

  bool time_rate_diagnostics; // Writes the time rate of change of field energy.
};

// Top-level app parameters: this
struct gkyl_gyrokinetic_multib {
  char name[128]; // name of app

  int cdim, vdim; // conf, velocity space dimensions
  int poly_order; // polynomial order
  enum gkyl_basis_type basis_type; // type of basis functions to use
  bool use_gpu; // Flag to indicate if solver should use GPUs

  // geometry and topology of all blocks in simulation
  struct gkyl_gk_block_geom *gk_block_geom;

  double cfl_frac; // CFL fraction to use (default 1.0)

  bool enforce_positivity; // Positivity enforcement via shift in f.

  int num_species; // number of species
  // species inputs
  struct gkyl_gyrokinetic_multib_species species[GKYL_MAX_SPECIES];

  int num_neut_species; // number of neutral species  
  // neutral species inputs
  struct gkyl_gyrokinetic_multib_neut_species neut_species[GKYL_MAX_SPECIES];

  bool skip_field; // Skip field update -> phi = 0 for all time  
  // field inputs
  struct gkyl_gyrokinetic_multib_field field;

  // communicator to use.  
  struct gkyl_comm *comm;  
};

/**
 * Construct a new gk multi-block app.
 *
 * @param mbinp App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New multi-block gk app object.
 */
gkyl_gyrokinetic_multib_app* gkyl_gyrokinetic_multib_app_new(const struct gkyl_gyrokinetic_multib *mbinp);

/**
 * Initialize species by projecting initial conditions on
 * basis functions, and the field by solving the field equations.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_gyrokinetic_multib_app_apply_ic(gkyl_gyrokinetic_multib_app* app, double t0);

/**
 * Initialize species by projecting initial conditions on basis
 * functions. Species index (sidx) is the same index used to specify
 * the species in the gkyl_gk object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_gyrokinetic_multib_app_apply_ic_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0);

/**
 * Initialize neutral species by projecting initial conditions on basis
 * functions. Neutral species index (sidx) is the same index used to specify
 * the neutral species in the gkyl_gk object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of neutral species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_gyrokinetic_multib_app_apply_ic_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0);

/**
 * Initialize field from file
 *
 * @param app App object
 * @param fname file to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_file_field(gkyl_gyrokinetic_multib_app *app, const char *fname);

/**
 * Initialize gyrokinetic species from file
 *
 * @param app App object
 * @param sidx gk species index
 * @param fname file to read
 */
struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname);

/**
 * Initialize neutral species from file
 *
 * @param app App object
 * @param sidx neut species index
 * @param fname file to read
 */
struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname);

/**
 * Initialize the gyrokinetic app from a specific frame.
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_read_from_frame(gkyl_gyrokinetic_multib_app *app, int frame);

/**
 * Initialize field from frame
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_field(gkyl_gyrokinetic_multib_app *app, int frame);

/**
 * Initialize gyrokinetic species from file
 *
 * @param app App object
 * @param sidx gk species index
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame);

/**
 * Initialize neutral species from file
 *
 * @param app App object
 * @param sidx neut species index
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame);

/**
 * Write output to console: this is mainly for diagnostic messages the
 * driver code wants to write to console. It accounts for parallel
 * output by not messing up the console with messages from each rank.
 *
 * @param app App object
 * @param fp File pointer for open file for output
 * @param fmt Format string for console output
 * @param argp Objects to write
 */
void gkyl_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, ...);

/**
 * Write block topology to file.
 * 
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_write_topo(const gkyl_gyrokinetic_multib_app* app);

/**
 * Write geometry file.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_write_geometry(gkyl_gyrokinetic_multib_app *app);

/**
 * Write field data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_field(gkyl_gyrokinetic_multib_app* app, double tm, int frame);

/**
 * Calculate integrated field energy
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_calc_field_energy(gkyl_gyrokinetic_multib_app* app, double tm);

/**
 * Write field energy to file. Field energy data is appended to the
 * same file.
 * 
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_write_field_energy(gkyl_gyrokinetic_multib_app* app);

/**
 * Write species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write neutral species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for neutral species to file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Calculate integrated diagnostic moments for a plasma species.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_species_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm);

/**
 * Calculate integrated diagnostic moments for a neutral species.
 *
 * @param app App object.
 * @param sidx Index of neutral species to initialize.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_neut_species_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm);

/**
 * Calculate integrated diagnostic moments of the boundary fluxes for a plasma species.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm);

/**
 * Write integrated diagnostic moments for charged species to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 */
void gkyl_gyrokinetic_multib_app_write_species_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write integrated diagnostic moments for neutral species to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species to initialize.
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write integrated diagnostic moments of the boundary fluxes for charged
 * species to file. Integrated moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 */
void gkyl_gyrokinetic_multib_app_write_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write species source to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_source(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write neutral species source to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_source(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for species source to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_source_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for neutral species source to file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_source_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Calculate integrated diagnostic moments for a plasma species source.
 *
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_species_source_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm);

/**
 * Calculate integrated diagnostic moments for a neutral species source.
 *
 * @param app App object.
 * @param sidx Index of neutral species to write.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_neut_species_source_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm);

/**
 * Write integrated diagnostic moments for charged species source to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 */
void gkyl_gyrokinetic_multib_app_write_species_source_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write integrated diagnostic moments for neutral species source to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species to write.
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_source_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write species integrated correct Maxwellian status of the to file. 
 * Correct Maxwellian status is appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 */
void gkyl_gyrokinetic_multib_app_write_species_lte_max_corr_status(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write neutral species integrated correct Maxwellian status of the to file. 
 * Correct Maxwellian status is appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_lte_max_corr_status(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write LBO collisional moments for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_lbo_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Write BGK cross moments for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_bgk_cross_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Write radiation drag coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_rad_drag(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Write radiation emissivity of each species that species sidx collides with
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_rad_emissivity(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Calculate integrated diagnostic moments of the radiation model.
 *
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_species_rad_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm);

/**
 * Write integrated moments of radiation rhs for radiating species 
 * 
 * @param app App object.
 * @param sidx Index of species from which to write radiation.
 */
void gkyl_gyrokinetic_multib_app_write_species_rad_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write iz react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_iz_react(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write iz react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_iz_react_neut(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write recomb react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_recomb_react(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write recomb react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_recomb_react_neut(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write cx react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_cx_react_neut(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write the phase-space diagnostics for a charged species.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_phase(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write the phase-space diagnostics for a neutral species.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_phase(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write the conf-space diagnostics for a charged species.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_conf(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write the conf-space diagnostics for a neutral species.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_conf(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for all species (including sources) to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_mom(gkyl_gyrokinetic_multib_app *app, double tm, int frame);

/**
 * Calculate integrated diagnostic moments for all species (including sources).
 *
 * @param app App object.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_integrated_mom(gkyl_gyrokinetic_multib_app* app, double tm);

/**
 * Write integrated diagnostic moments for all species (including sources)
 * to file. Integrated moments are appended to the same file.
 * 
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_write_integrated_mom(gkyl_gyrokinetic_multib_app *app);

/**
 * Write phase space diagnostics to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_phase(gkyl_gyrokinetic_multib_app* app, double tm, int frame);

/**
 * Write configuration space diagnostics to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_conf(gkyl_gyrokinetic_multib_app* app, double tm, int frame);

/**
 * Write both conf and phase-space diagnostics to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write(gkyl_gyrokinetic_multib_app* app, double tm, int frame);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_stat_write(gkyl_gyrokinetic_multib_app* app);

/**
 * Print timing of solver components to iostream.
 *
 * @param app App object.
 * @param iostream Where to write timers to (e.g. stdout, stderr);
 */
void
gkyl_gyrokinetic_multib_app_print_timings(gkyl_gyrokinetic_multib_app* app, FILE *iostream);

/**
 * Record the time step (in private dynvector).
 *
 * @param app App object.
 * @param tm Time stamp.
 * @param dt Time step to record (e.g. provided by app's status object).
 */
void
gkyl_gyrokinetic_multib_app_save_dt(gkyl_gyrokinetic_multib_app* app, double tm, double dt);

/**
 * Write the time step over time.
 *
 * @param app App object.
 */
void
gkyl_gyrokinetic_multib_app_write_dt(gkyl_gyrokinetic_multib_app* app);

/**
 * Read geometry file.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_read_geometry(gkyl_gyrokinetic_multib_app *app);

/**
 * Write output to console: this is mainly for diagnostic messages the
 * driver code wants to write to console. It accounts for parallel
 * output by not messing up the console with messages from each rank.
 *
 * @param app App object
 * @param fp File pointer for open file for output
 * @param fmt Format string for console output
 * @param argp Objects to write
 */
void gkyl_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, ...);

/**
 * Advance simulation by a suggested time-step 'dt'. The dt may be too
 * large in which case method will attempt to take a smaller time-step
 * and also return it as the 'dt_actual' field of the status
 * object. If the suggested time-step 'dt' is smaller than the largest
 * stable time-step the method will use the smaller value instead,
 * returning the larger time-step in the 'dt_suggested' field of the
 * status object. If the method fails to find any stable time-step
 * then the 'success' flag will be set to 0. At that point the calling
 * code must abort the simulation as this signals a catastrophic
 * failure and the simulation can't be safely continued.
 * 
 * @param app App object.
 * @param dt Suggested time-step to advance simulation
 * @return Status of update.
 */
struct gkyl_update_status gkyl_gyrokinetic_multib_update(gkyl_gyrokinetic_multib_app* app, double dt);

/**
 * Return simulation statistics.
 * 
 * @return Return statistics object.
 */
struct gkyl_gyrokinetic_stat gkyl_gyrokinetic_multib_app_stat(gkyl_gyrokinetic_multib_app* app);

/**
 * Run the RHS for the species update. This is used to compute kernel
 * timers and is not otherwise a useful function for a full
 * simulation.
 *
 * @param app App object.
 * @param update_vol_term Set to 1 to update vol term also, 0 otherwise
 */
void gkyl_gyrokinetic_multib_app_species_ktm_rhs(gkyl_gyrokinetic_multib_app* app, int update_vol_term);

/**
 * Free gk app.
 *
 * @param app App to release.
 */
void gkyl_gyrokinetic_multib_app_release(gkyl_gyrokinetic_multib_app* app);
