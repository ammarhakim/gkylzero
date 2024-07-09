// Private header for use in MB moment app: do not include in
// user-facing header files!
#pragma once

#include <gkyl_moment_priv.h>

// top-level internal App
struct gkyl_moment_multib_app {
  char name[128]; // name of app
  struct gkyl_comm *comm; // communicator to use
  
 // geometry and topology of all blocks in simulation
  struct gkyl_block_geom *block_geom;
  double cfl_frac; // CFL fraction to use
  int num_species; // number of species
  
  int num_local_blocks; // total number of blocks on current rank
  struct gkyl_moment_app **app; // individual App objects, one per-block

  struct gkyl_moment_stat stat; // statistics
};
