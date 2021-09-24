#pragma once

#include <gkyl_app.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_util.h>

#include <stdbool.h>

// Parameters for Vlasov species
struct gkyl_vlasov_species {
  char name[128]; // species name
  double charge, mass; // charge and mass
  double lower[3], upper[3]; // lower, upper bounds of velocity-space
  int cells[3]; // velocity-space cells

  bool evolve; // evolve species? 1-yes, 0-no

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  int num_diag_moments; // number of diagnostic moments
  char diag_moments[16][16]; // list of diagnostic moments

  double nu; // collision frequency
  enum gkyl_collision_id collision_id; // type of collisions (see gkyl_eqn_type.h)
};

// Parameter for EM field
struct gkyl_vlasov_field {
  enum gkyl_field_id field_id; // type of field (see gkyl_eqn_type.h)
  bool evolve; // evolve field? 1-yes, 0-no
  
  double epsilon0, mu0;
  double elcErrorSpeedFactor, mgnErrorSpeedFactor;

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);
};

// Top-level app parameters
struct gkyl_vm {
  char name[128]; // name of app: used as output prefix

  int cdim, vdim; // conf, velocity space dimensions
  double lower[3], upper[3]; // lower, upper bounds of config-space
  int cells[3]; // config-space cells
  int poly_order; // polynomial order
  enum gkyl_basis_type basis_type; // type of basis functions to use

  double cfl_frac; // CFL fraction to use (default 1.0)

  bool use_gpu; // Flag to indicate if solver should use GPUs

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  int num_species; // number of species
  struct gkyl_vlasov_species species[GKYL_MAX_SPECIES]; // species objects

  bool skip_field; // Skip field update or no field specified
  struct gkyl_vlasov_field field; // field object
};

// Simulation statistics
struct gkyl_vlasov_stat {
  long nup; // calls to update
  long nfeuler; // calls to forward-Euler method
    
  long nstage_2_fail; // number of failed RK stage-2s
  long nstage_3_fail; // number of failed RK stage-3s

  double stage_2_dt_diff[2]; // [min,max] rel-diff for stage-2 failure
  double stage_3_dt_diff[2]; // [min,max] rel-diff for stage-3 failure
    
  double total_tm; // time for simulation (not including ICs)
  double init_species_tm; // time to initialize all species
  double init_field_tm; // time to initialize fields

  double species_rhs_tm; // time to compute species RHS
  double field_rhs_tm; // time to compute field RHS
  double current_tm; // time to compute currents and accumulation

  long nmom; // calls to moment calculation
  double mom_tm; // time to compute moments
};

// Object representing Vlasov app
typedef struct gkyl_vlasov_app gkyl_vlasov_app;

/**
 * Construct a new Vlasov app.
 *
 * @param vm App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New vlasov app object.
 */
gkyl_vlasov_app* gkyl_vlasov_app_new(struct gkyl_vm vm);

/**
 * Initialize species and field by projecting initial conditions on
 * basis functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_vlasov_app_apply_ic(gkyl_vlasov_app* app, double t0);

/**
 * Initialize field by projecting initial conditions on basis
 * functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions
 */
void gkyl_vlasov_app_apply_ic_field(gkyl_vlasov_app* app, double t0);

/**
 * Initialize species by projecting initial conditions on basis
 * functions. Species index (sidx) is the same index used to specify
 * the species in the gkyl_vm object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_vlasov_app_apply_ic_species(gkyl_vlasov_app* app, int sidx, double t0);

/**
 * Calculate diagnostic moments.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_calc_mom(gkyl_vlasov_app* app);

/**
 * Write field and species data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write(gkyl_vlasov_app* app, double tm, int frame);

/**
 * Write field data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_field(gkyl_vlasov_app* app, double tm, int frame);

/**
 * Write species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_species(gkyl_vlasov_app* app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for species to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_mom(gkyl_vlasov_app* app, double tm, int frame);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_stat_write(const gkyl_vlasov_app* app);

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
struct gkyl_update_status gkyl_vlasov_update(gkyl_vlasov_app* app, double dt);

/**
 * Return simulation statistics.
 * 
 * @return Return statistics object.
 */
struct gkyl_vlasov_stat gkyl_vlasov_app_stat(gkyl_vlasov_app* app);

/**
 * Run the RHS for the species update. This is used to compute kernel
 * timers and is not otherwise a useful function for a full
 * simulation.
 *
 * @param app App object.
 * @param update_vol_term Set to 1 to update vol term also, 0 otherwise
 */
void gkyl_vlasov_app_species_ktm_rhs(gkyl_vlasov_app* app, int update_vol_term);

/**
 * Free Vlasov app.
 *
 * @param app App to release.
 */
void gkyl_vlasov_app_release(gkyl_vlasov_app* app);
