// Private header for use in MB moment app: do not include in
// user-facing header files!
#pragma once

#include <gkyl_moment_priv.h>
#include <gkyl_comm.h>

// top-level internal App
struct gkyl_moment_multib_app {
  char name[128]; // name of app
  struct gkyl_comm *comm; // global communicator to use
  
 // geometry and topology of all blocks in simulation
  struct gkyl_block_geom *block_geom;
  struct gkyl_block_topo *block_topo;
  
  double cfl_frac; // CFL fraction to use
  int num_species; // number of species

  struct gkyl_comm **block_comms; // list of block-communicators

  bool has_field; // true if there is a field present
  char species_name[GKYL_MAX_SPECIES][128]; // name of each species

  int num_local_blocks; // total number of blocks on current rank
  int *local_blocks; // local blocks IDs handled by current rank
  struct gkyl_moment_app **singleb_apps; // App objects: one per local block

  struct gkyl_rect_decomp **decomp; // list of decomps (num_blocks)

  double tcurr; // current time
  
  struct gkyl_moment_stat stat; // statistics
};

// Meta-data for IO
struct moment_multib_output_meta {
  int frame; // frame number
  double stime; // output time
  const char *app_name; // name of App
  const char *topo_file_name; // name of topology file
};
