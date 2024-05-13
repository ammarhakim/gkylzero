#pragma once

#include <gkyl_app.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#include <stdbool.h>

// Parameters for projection
struct gkyl_vlasov_poisson_projection {
  enum gkyl_projection_id proj_id; // type of projection (see gkyl_eqn_type.h)

  union {
    struct {
      // pointer and context to initialization function 
      void *ctx_func; 
      void (*func)(double t, const double *xn, double *fout, void *ctx); 
    };
    struct {
      // pointers and contexts to initialization functions for LTE distribution projection
      // (Maxwellian for non-relativistic, Maxwell-Juttner for relativistic)
      void *ctx_density;
      void (*density)(double t, const double *xn, double *fout, void *ctx);
      void *ctx_V_drift;
      void (*V_drift)(double t, const double *xn, double *fout, void *ctx);
      void *ctx_temp;
      void (*temp)(double t, const double *xn, double *fout, void *ctx);

      // boolean if we are correcting all the moments or only density
      bool correct_all_moms;       
    };
  };
};

// Parameters for species collisions
struct gkyl_vlasov_poisson_collisions {
  enum gkyl_collision_id collision_id; // type of collisions (see gkyl_eqn_type.h)

  void *ctx; // context for collision function
  // function for computing self-collision frequency
  void (*self_nu)(double t, const double *xn, double *fout, void *ctx);

  // nputs for Spitzer collisionality
  bool normNu; // Set to true if you want to rescale collision frequency
  double n_ref; // Density used to calculate coulomb logarithm
  double T_ref; // Temperature used to calculate coulomb logarithm
  double nuFrac; // Parameter for rescaling collision frequency from SI values
  double hbar, eps0, eV; // Planck's constant/2 pi, vacuum permittivity, elementary charge

  // boolean if we are correcting all the moments or only density:
  // only used by BGK collisions
  bool correct_all_moms;
  double iter_eps; // error tolerance for moment fixes (density is always exact)
  int max_iter; // maximum number of iteration

  int num_cross_collisions; // number of species to cross-collide with
  char collide_with[GKYL_MAX_SPECIES][128]; // names of species to cross collide with
};

// Parameters for species source
struct gkyl_vlasov_poisson_source {
  enum gkyl_source_id source_id; // type of source
  bool write_source;

  double source_length; // required for boundary flux source
  char source_species[128];
  
  // sources using projection routine
  struct gkyl_vlasov_poisson_projection projection;
};

// Parameters for boundary conditions
struct gkyl_vlasov_poisson_bc {
  enum gkyl_species_bc_type type;
  void *aux_ctx;
  void (*aux_profile)(double t, const double *xn, double *fout, void *ctx);
  double aux_parameter;
};

struct gkyl_vlasov_poisson_bcs {
  struct gkyl_vlasov_poisson_bc lower, upper;
};

// Parameters for Vlasov species
struct gkyl_vlasov_poisson_species {
  char name[128]; // species name

  enum gkyl_vpmodel_id model_id; // type of model, see gkyl_eqn_type.h

  double charge, mass; // charge and mass
  double lower[3], upper[3]; // lower, upper bounds of velocity-space
  int cells[3]; // velocity-space cells

  // initial conditions using projection routine
  struct gkyl_vlasov_poisson_projection projection;

  int num_diag_moments; // number of diagnostic moments
  char diag_moments[16][16]; // list of diagnostic moments

  // collisions to include
  struct gkyl_vlasov_poisson_collisions collisions;

  // source to include
  struct gkyl_vlasov_poisson_source source;

  // boundary conditions
  struct gkyl_vlasov_poisson_bcs bcx, bcy, bcz;
};

// Parameter for VP field.
struct gkyl_vlasov_poisson_field {
  enum gkyl_vpfield_id field_id; // Type of field.

  double permittivity;

  struct gkyl_poisson_bc poisson_bcs;

  void *external_potentials_ctx; // Context for external electromagnetic fields function
  // Function definiting external potentials.
  void (*external_potentials)(double t, const double *xn, double *ext_em_out, void *ctx);

  bool external_potentials_evolve; // Set to true if external potentials are time dependent.
};

// Top-level app parameters
struct gkyl_vp {
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
  struct gkyl_vlasov_poisson_species species[GKYL_MAX_SPECIES]; // species objects

  bool skip_field; // Skip field update or no field specified
  struct gkyl_vlasov_poisson_field field; // field object

  // this should not be set by typical user-facing code but only by
  // higher-level drivers
  bool has_low_inp; // should one use low-level inputs?
  struct gkyl_app_comm_low_inp low_inp; // low-level inputs  
};

// Simulation statistics
struct gkyl_vlasov_poisson_stat {
  bool use_gpu; // did this sim use GPU?
  
  long nup; // calls to update
  long nfeuler; // calls to forward-Euler method
    
  long nstage_2_fail; // number of failed RK stage-2s
  long nstage_3_fail; // number of failed RK stage-3s

  double stage_2_dt_diff[2]; // [min,max] rel-diff for stage-2 failure
  double stage_3_dt_diff[2]; // [min,max] rel-diff for stage-3 failure
    
  double total_tm; // time for simulation (not including ICs)
  double init_species_tm; // time to initialize all species
  double init_field_tm; // time to initialize fields

  double species_rhs_tm; // time to compute species collisionless RHS
  double field_rhs_tm; // time to compute field RHS
    
  double species_coll_mom_tm; // time needed to compute various moments needed in LBO
  double species_lbo_coll_drag_tm[GKYL_MAX_SPECIES]; // time to compute LBO drag terms
  double species_lbo_coll_diff_tm[GKYL_MAX_SPECIES]; // time to compute LBO diffusion terms
  double species_coll_tm; // total time for collision updater (excluded moments)

  double species_bc_tm; // time to compute species BCs
  double field_bc_tm; // time to compute field

  long nspecies_omega_cfl; // number of times CFL-omega all-reduce is called
  double species_omega_cfl_tm; // time spent in all-reduce for omega-cfl

  long nfield_omega_cfl; // number of times CFL-omega for field all-reduce is called
  double field_omega_cfl_tm; // time spent in all-reduce for omega-cfl for field

  long nmom; // calls to moment calculation
  double mom_tm; // time to compute moments

  long ndiag; // calls to diagnostics
  double diag_tm; // time to compute diagnostics

  long nio; // number of calls to IO
  double io_tm; // time to perform IO
};

// Object representing Vlasov app
typedef struct gkyl_vlasov_poisson_app gkyl_vlasov_poisson_app;

/**
 * Construct a new Vlasov app.
 *
 * @param vp App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New vlasov_poisson app object.
 */
gkyl_vlasov_poisson_app* gkyl_vlasov_poisson_app_new(struct gkyl_vp *vp);

/**
 * Initialize species and field by projecting initial conditions on
 * basis functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_vlasov_poisson_app_apply_ic(gkyl_vlasov_poisson_app* app, double t0);

/**
 * Initialize species by projecting initial conditions on basis
 * functions. Species index (sidx) is the same index used to specify
 * the species in the gkyl_vp object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_vlasov_poisson_app_apply_ic_species(gkyl_vlasov_poisson_app* app, int sidx, double t0);

/**
 * Initialize field from file
 *
 * @param app App object
 * @param fname file to read
 */
struct gkyl_app_restart_status
gkyl_vlasov_poisson_app_from_file_field(gkyl_vlasov_poisson_app *app, const char *fname);

/**
 * Initialize Vlasov-Poisson species from file
 *
 * @param app App object
 * @param sidx gk species index
 * @param fname file to read
 */
struct gkyl_app_restart_status
gkyl_vlasov_poisson_app_from_file_species(gkyl_vlasov_poisson_app *app, int sidx,
  const char *fname);

/**
 * Initialize the Vlasov-Poisson app from a specific frame.
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_vlasov_poisson_app_read_from_frame(gkyl_vlasov_poisson_app *app, int frame);

/**
 * Initialize field from frame
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_vlasov_poisson_app_from_frame_field(gkyl_vlasov_poisson_app *app, int frame);

/**
 * Initialize Vlasov-Poisson species from file
 *
 * @param app App object
 * @param sidx gk species index
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_vlasov_poisson_app_from_frame_species(gkyl_vlasov_poisson_app *app, int sidx, int frame);

/**
 * Calculate diagnostic moments.
 *
 * @param app App object.
 */
void gkyl_vlasov_poisson_app_calc_mom(gkyl_vlasov_poisson_app *app);

/**
 * Calculate integrated diagnostic moments.
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_vlasov_poisson_app_calc_integrated_mom(gkyl_vlasov_poisson_app* app, double tm);

/**
 * Calculate integrated field energy
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_vlasov_poisson_app_calc_field_energy(gkyl_vlasov_poisson_app* app, double tm);

/**
 * Write field and species data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_poisson_app_write(gkyl_vlasov_poisson_app* app, double tm, int frame);

/**
 * Write field data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_poisson_app_write_field(gkyl_vlasov_poisson_app* app, double tm, int frame);

/**
 * Write species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_poisson_app_write_species(gkyl_vlasov_poisson_app* app, int sidx, double tm, int frame);

/**
 * Write source species data to file.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_poisson_app_write_source_species(gkyl_vlasov_poisson_app* app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for species to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_poisson_app_write_mom(gkyl_vlasov_poisson_app *app, double tm, int frame);

/**
 * Write diagnostic moments for species source to file.
 *
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_poisson_app_write_source_mom(gkyl_vlasov_poisson_app *app, double tm, int frame);

/**
 * Write integrated diagnostic moments for species to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 */
void gkyl_vlasov_poisson_app_write_integrated_mom(gkyl_vlasov_poisson_app *app);

/**
 * Write integrated diagnostic moments for sources to file. Integrated
 * moments are appended to the same file.
 *
 * @param app App object.
 */
void gkyl_vlasov_poisson_app_write_integrated_source_mom(gkyl_vlasov_poisson_app *app);

/**
 * Write integrated correct lte status of the species distribution function to file. Correct
 * lte status is appended to the same file.
 * 
 * @param app App object.
 */
void gkyl_vlasov_poisson_app_write_lte_corr_status(gkyl_vlasov_poisson_app *app);

/**
 * Write field energy to file. Field energy data is appended to the
 * same file.
 * 
 * @param app App object.
 */
void gkyl_vlasov_poisson_app_write_field_energy(gkyl_vlasov_poisson_app* app);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_vlasov_poisson_app_stat_write(gkyl_vlasov_poisson_app* app);

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
void gkyl_vlasov_poisson_app_cout(const gkyl_vlasov_poisson_app* app, FILE *fp, const char *fmt, ...);

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
struct gkyl_update_status gkyl_vlasov_poisson_update(gkyl_vlasov_poisson_app* app, double dt);

/**
 * Return simulation statistics.
 * 
 * @return Return statistics object.
 */
struct gkyl_vlasov_poisson_stat gkyl_vlasov_poisson_app_stat(gkyl_vlasov_poisson_app* app);

/**
 * Run the RHS for the species update. This is used to compute kernel
 * timers and is not otherwise a useful function for a full
 * simulation.
 *
 * @param app App object.
 * @param update_vol_term Set to 1 to update vol term also, 0 otherwise
 */
void gkyl_vlasov_poisson_app_species_ktm_rhs(gkyl_vlasov_poisson_app* app, int update_vol_term);

/**
 * Free Vlasov app.
 *
 * @param app App to release.
 */
void gkyl_vlasov_poisson_app_release(gkyl_vlasov_poisson_app* app);
