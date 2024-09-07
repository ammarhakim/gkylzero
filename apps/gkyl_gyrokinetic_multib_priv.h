// Private header for use in MB gyrokinetic app: do not include in
// user-facing header files!
#pragma once

#include <gkyl_comm.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_multib_comm_conn.h>

// top-level internal App
struct gkyl_gyrokinetic_multib_app {
  char name[128]; // name of app
  struct gkyl_comm *comm; // global communicator to use
  
 // geometry and topology of all blocks in simulation
  struct gkyl_block_geom *block_geom;
  struct gkyl_block_topo *block_topo;
  
  double cfl_frac; // CFL fraction to use
  int num_species; // number of species
  int num_neut_species; // number of neutral species
  bool update_field; // true if there solving Poisson equation

  char species_name[GKYL_MAX_SPECIES][128]; // name of each species
  char neut_species_name[GKYL_MAX_SPECIES][128]; // name of each neutral species  

  struct gkyl_comm **block_comms; // list of block-communicators

  int num_local_blocks; // total number of blocks on current rank
  int *local_blocks; // local blocks IDs handled by current rank
  struct gkyl_gyrokinetic_app **singleb_apps; // App objects: one per local block

  // one per local block:
  struct gkyl_multib_comm_conn **send_conn; // conn for inter-block send
  struct gkyl_multib_comm_conn **recv_conn; // conn for inter-block recv

  const struct gkyl_rrobin_decomp *round_robin; // round-robin decomp
  struct gkyl_rect_decomp **decomp; // list of decomps (num_blocks)

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


// Take a single time-step using a SSP-RK3 stepper
static struct gkyl_update_status gyrokinetic_multib_update_ssp_rk3(struct gkyl_gyrokinetic_multib_app* mbapp, double dt0);

// Take a forward Euler step with the suggested time-step dt. This may
// not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by
// stability. The actual time-step and dt_suggested are returned in
// the status object.
void gyrokinetic_multib_forward_euler(struct gkyl_gyrokinetic_multib_app* mbapp, double tcurr, double dt,
  const struct gkyl_array *fin[], struct gkyl_array *fout[],
  const struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[],
  struct gkyl_update_status *st, struct gkyl_update_status *sb_st);
