// Private header file for multi-block gyrokinetic solver.
#pragma once

#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gyrokinetic_mb.h>

// Top-level multiblock gyrokinetic App object.
struct gkyl_gyrokinetic_mb_app {
  char name[128]; // Name of app.
  struct gkyl_job_pool *job_pool; // Job pool.

  int cdim, vdim; // Conf, velocity space dimensions.
  int poly_order; // Polynomial order.
  double cfl; // CFL number.

  bool use_gpu; // Should we use GPU (if present).

  struct gkyl_block_topo *btopo; // Block topology.

  int num_blocks, num_blocks_local; // Total and local number of blocks.
  int *block_idxs; // List of blocks handled on this rank.
  struct gkyl_gyrokinetic_app **blocks; // Pointers to blocks on this rank.

  struct gkyl_range local, local_ext; // local, local-ext conf-space ranges. Local means individual blokck range.
  struct gkyl_range global, global_ext; // global, global-ext conf-space ranges  . Global means cross-block range.
  struct gkyl_basis confBasis; // conf-space basis

  struct gk_field_mb *field; // mb field object.
  struct gkyl_comm *comm;   // communicator object for conf-space arrays
};

// field data
struct gk_field_mb {
  struct gkyl_gyrokinetic_field info; // data for field

  enum gkyl_gkfield_id gkfield_id;

  struct gkyl_job_pool *job_pool; // Job pool  
  // arrays for global charge density and global smoothed (in z) charge density
  struct gkyl_array *rho_c_global_dg;
  struct gkyl_array *rho_c_global_smooth; 
  struct gkyl_array *phi_fem; // arrays for updates

  struct gkyl_array *phi_host;  // host copy for use IO and initialization

  struct gkyl_range global_sub_range; // sub range of intersection of global range and local range
                                      // for solving subset of Poisson solves with parallelization in z

  struct gkyl_array *weight;
  struct gkyl_fem_parproj *fem_parproj; // FEM smoother for projecting DG functions onto continuous FEM basis
                                        // weight*phi_{fem} = phi_{dg} 
  struct gkyl_deflated_fem_poisson *deflated_fem_poisson; // poisson solver which solves on lines in x or planes in xy
                                                          // - nabla . (epsilon * nabla phi) - kSq * phi = rho
};

