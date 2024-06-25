#pragma once

#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_comms.h>

// Input struct to create a multiblock gyrokinetic app.
struct gkyl_gk_multib {
  char name[128]; // name of app: used as output prefix

  int cdim, vdim; // conf, velocity space dimensions
  int poly_order; // polynomial order
  enum gkyl_basis_type basis_type; // type of basis functions to use

  double cfl_frac; // CFL fraction to use (default 1.0)

  bool use_mpi; // Flag to indicate if solver should use MPI.
  bool use_gpu; // Flag to indicate if solver should use GPUs.

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  int num_species; // number of species
  struct gkyl_gyrokinetic_species species[GKYL_MAX_SPECIES]; // species objects

  int num_neut_species; // number of species
  struct gkyl_gyrokinetic_neut_species neut_species[GKYL_MAX_SPECIES]; // species objects
  
  bool skip_field;
  struct gkyl_gyrokinetic_field field; // field object

  int num_blocks;
  struct gkyl_gk *blocks[GKYL_MAX_BLOCKS];
};

// Object representing gk app
typedef struct gkyl_gyrokinetic_multib_app gkyl_gyrokinetic_multib_app;

/**
 * Construct a new multiblock GK app.
 *
 * @param gk_multib App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New gk app object.
 */
gkyl_gyrokinetic_multib_app* gkyl_gyrokinetic_multib_app_new(struct gkyl_gk_multib *inp);

/**
 * Initialize species and field by projecting initial conditions on
 * basis functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_gyrokinetic_multib_app_apply_ic(gkyl_gyrokinetic_multib_app* mba, double t0);

/**
 * Calculate diagnostic moments.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_calc_mom(gkyl_gyrokinetic_multib_app* mba);

/**
 * Calculate integrated diagnostic moments for plasma species (including sources).
 *
 * @param tm Time at which integrated diagnostics are to be computed
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_calc_integrated_mom(gkyl_gyrokinetic_multib_app* mba, double tm);

/**
 * Write diagnostic moments for species to file.
 *
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_mom(gkyl_gyrokinetic_multib_app* mba, double tm, int frame);

/**
 * Write integrated diagnostic moments for species to file. Integrated
 * moments are appended to the same file.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_write_integrated_mom(gkyl_gyrokinetic_multib_app* mba);

/**
 * Write field and species data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write(gkyl_gyrokinetic_multib_app* mba, double tm, int frame);

/**
 * Initialize the gyrokinetic app from a specific frame.
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_read_from_frame(gkyl_gyrokinetic_multib_app* mba, int frame);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_stat_write(gkyl_gyrokinetic_multib_app* mba);

/**
 * Return simulation statistics.
 *
 * @return Return statistics object.
 */
struct gkyl_gyrokinetic_stat gkyl_gyrokinetic_multib_app_stat(gkyl_gyrokinetic_multib_app* mba);

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
void gkyl_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* mba, FILE *fp, const char *fmt, ...);

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
struct gkyl_update_status gkyl_gyrokinetic_multib_update(gkyl_gyrokinetic_multib_app* mba, double dt);

/**
 * Free a multiblock GK app.
 *
 * @param app App to release.
 */
void gkyl_gyrokinetic_multib_app_release(gkyl_gyrokinetic_multib_app* mba);
