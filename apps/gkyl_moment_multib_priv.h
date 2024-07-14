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

  int num_local_blocks; // total number of blocks on current rank
  int *local_blocks; // local blocks IDs handled by current rank
  struct gkyl_moment_app **singleb_apps; // individual App objects, one per-block

  double tcurr; // current time
  
  struct gkyl_moment_stat stat; // statistics
};
