#pragma once

#include <gkyl_app.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wave_prop.h>

#include <time.h>

// Boundary conditions on fields and fluids
enum gkyl_moment_bc_type {
  GKYL_MOMENT_COPY = 0, // copy BCs for fluid and field
  GKYL_MOMENT_SPECIES_WALL, // perfect reflector for moments
  GKYL_MOMENT_FIELD_COND, // perfect conductor for fields
};

// Parameters for moment species
struct gkyl_moment_species {
  char name[128]; // species name
  double charge, mass; // charge and mass
  double k0; // closure parameter 
  enum gkyl_wave_limiter limiter; // limiter to use
  const struct gkyl_wv_eqn *equation; // equation object

  int evolve; // evolve species? 1-yes, 0-no

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  // boundary conditions
  enum gkyl_moment_bc_type bcx[2], bcy[2], bcz[2];
};

// Parameter for EM field
struct gkyl_moment_field {
  double epsilon0, mu0;
  double elc_error_speed_fact, mag_error_speed_fact;

  enum gkyl_wave_limiter limiter; // limiter to use

  int evolve; // evolve field? 1-yes, 0-no

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  // boundary conditions
  enum gkyl_moment_bc_type bcx[2], bcy[2], bcz[2];
};

// Choices of schemes to use in the fluid solver 
enum gkyl_moment_fluid_scheme {
  GKYL_MOMENT_FLUID_WAVE_PROP = 0, // default
  GKYL_MOMENT_FLUID_KEP
};

// Top-level app parameters
struct gkyl_moment {
  char name[128]; // name of app: used as output prefix

  int ndim; // space dimensions
  double lower[3], upper[3]; // lower, upper bounds
  int cells[3]; // config-space cells

  double cfl_frac; // CFL fraction to use

  enum gkyl_moment_fluid_scheme fluid_scheme; // scheme to update fluid equations

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  int num_species; // number of species
  struct gkyl_moment_species species[GKYL_MAX_SPECIES]; // species objects
  struct gkyl_moment_field field; // field object
};

// Simulation statistics
struct gkyl_moment_stat {
  long nup; // calls to update
  long nfail; // number of failed time-steps

  double total_tm; // time for simulation (not including ICs)
  double species_tm; // time to compute species updates
  double field_tm; // time to compute field updates
  double sources_tm; // time to compute source terms
};

// Object representing moments app
typedef struct gkyl_moment_app gkyl_moment_app;

/**
 * Construct a new moments app.
 *
 * @param vm App inputs. See struct docs.
 * @return New moment app object.
 */
gkyl_moment_app* gkyl_moment_app_new(struct gkyl_moment mom);

/**
 * Compute maximum estimated stable dt wtih current app state. Call
 * after app initialized and after initial conditions set.
 *
 * @param app App object.
 * @retuen maximum estimated stable dt
 */
double gkyl_moment_app_max_dt(gkyl_moment_app* app);

/**
 * Initialize species and field.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_moment_app_apply_ic(gkyl_moment_app* app, double t0);

/**
 * Initialize field.
 *
 * @param app App object.
 * @param t0 Time for initial conditions
 */
void gkyl_moment_app_apply_ic_field(gkyl_moment_app* app, double t0);

/**
 * Initialize species.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_moment_app_apply_ic_species(gkyl_moment_app* app, int sidx, double t0);

/**
 * Write field and species data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_moment_app_write(const gkyl_moment_app* app, double tm, int frame);

/**
 * Write field data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_moment_app_write_field(const gkyl_moment_app* app, double tm, int frame);

/**
 * Write species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_moment_app_write_species(const gkyl_moment_app* app, int sidx, double tm, int frame);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_moment_app_stat_write(const gkyl_moment_app* app);

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
struct gkyl_update_status gkyl_moment_update(gkyl_moment_app* app, double dt);

/**
 * Return simulation statistics.
 * 
 * @return Return statistics.
 */
struct gkyl_moment_stat gkyl_moment_app_stat(gkyl_moment_app* app);

/**
 * Free moment app.
 *
 * @param app App to release.
 */
void gkyl_moment_app_release(gkyl_moment_app* app);
