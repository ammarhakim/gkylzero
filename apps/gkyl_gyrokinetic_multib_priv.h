// Private header for use in MB gyrokinetic app: do not include in
// user-facing header files!
#pragma once

#include <gkyl_gyrokinetic_priv.h>

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

  struct gkyl_comm **block_comms; // list of block-communicators

  int num_local_blocks; // total number of blocks on current rank
  int *local_blocks; // local blocks IDs handled by current rank
  struct gkyl_gyrokinetic_app **singleb_apps; // individual App objects, one per-block

  double tcurr; // current time
  
  struct gkyl_gyrokinetic_stat stat; // statistics
};
