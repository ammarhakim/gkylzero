// Private header for use in MB gyrokinetic app: do not include in
// user-facing header files!
#pragma once

#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_rrobin_decomp.h>
#include <gkyl_multib_comm_conn.h>

// A multib_comm_conn send/recv pair used for communicating between blocks.
struct gkyl_mbcc_sr {
  struct gkyl_multib_comm_conn **recv, **send;
};

// top-level internal App
struct gkyl_gyrokinetic_multib_app {
  char name[128]; // name of app
  struct gkyl_comm *comm; // global communicator to use
  bool use_gpu; // Whether to use the GPU.
  
  // geometry and topology of all blocks in simulation
  struct gkyl_block_geom *block_geom;
  struct gkyl_block_topo *block_topo;
  
  double cfl_frac; // CFL fraction to use
  int num_species; // number of species
  int num_neut_species; // number of neutral species

  bool update_field; // true if there solving Poisson equation
  struct gk_field *field; // Field object. MF 2024/10/20: will replace this
                          // with MB field object.

  char species_name[GKYL_MAX_SPECIES][128]; // name of each species
  char neut_species_name[GKYL_MAX_SPECIES][128]; // name of each neutral species  

  struct gkyl_comm **block_comms; // list of block-communicators

  int num_local_blocks; // total number of blocks on current rank
  int *local_blocks; // local blocks IDs handled by current rank
  struct gkyl_gyrokinetic_app **singleb_apps; // App objects: one per local block

  const struct gkyl_rrobin_decomp *round_robin; // round-robin decomp
  struct gkyl_rect_decomp **decomp; // list of decomps (num_blocks)

  struct gkyl_mbcc_sr *mbcc_sync_conf; // Connections for conf-space sync.
  struct gkyl_mbcc_sr *mbcc_sync_charged; // Connections for charged species phase-space.
  struct gkyl_mbcc_sr *mbcc_sync_neut; // Connections for neut species phase-space.

  double tcurr; // current time
  
  struct gkyl_gyrokinetic_stat stat; // statistics
};

// Meta-data for IO
struct gyrokinetic_multib_output_meta {
  int frame; // frame number
  double stime; // output time
  const char *app_name; // name of App
  const char *topo_file_name; // name of topology file
};

/** Time stepping API */

/**
 * Compute the gyrokinetic fields.
 *
 * @param app Multiblock gyrokinetic app.
 * @param tcurr Current simulation time.
 * @param fin Array of distribution functions (one for each species) .
 */
void gyrokinetic_multib_calc_field(struct gkyl_gyrokinetic_multib_app* app, double tcurr, const struct gkyl_array *fin[]);

/**
 * Compute the gyrokinetic fields and apply boundary conditions.
 *
 * @param app Gyrokinetic app.
 * @param tcurr Current simulation time.
 * @param distf Array of distribution functions (for each charged species).
 * @param distf_neut Array of distribution functions (for each neutral species).
 */
void gyrokinetic_multib_calc_field_and_apply_bc(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array *distf_neut[]);

/**
 * Take time-step using the RK3 method. Also sets the status object
 * which has the actual and suggested dts used. These can be different
 * from the actual time-step.
 *
 * @param app Gyrokinetic app.
 * @param dt0 Suggessted time step.
 */
struct gkyl_update_status gyrokinetic_multib_update_ssp_rk3(struct gkyl_gyrokinetic_multib_app* app, double dt0);
