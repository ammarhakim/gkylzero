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
  double tcurr; // Current time.
  double cfl; // CFL number.

  bool use_mpi; // Should we use MPI (if present).
  bool use_gpu; // Should we use GPU (if present).

  struct gkyl_comm *comm_mb;   // Multiblock communicator object.
  struct gkyl_comm **comm_intrab;
  struct gkyl_rect_decomp **decomp_intrab;

  struct gkyl_block_topo *btopo; // Block topology.

  int num_blocks, num_blocks_local; // Total and local number of blocks.
  int *block_idxs; // List of blocks handled on this rank.
  struct gkyl_gyrokinetic_app **blocks; // Pointers to blocks on this rank.

  struct gk_field_mb **field; // mb field object.
  struct gkyl_gyrokinetic_stat stat; // statistics
};

// field data
struct gk_field_mb {
  struct gkyl_gyrokinetic_field info; // data for field

  enum gkyl_gkfield_id gkfield_id;

  // z ranges
  struct gkyl_range crossz, crossz_ext; // crossz, crossz-ext conf-space ranges. Cross-block ranges across z boundaries for the smoother. We need 5 of these in the double null case (2-3-4, 11-12, 1-8-9, 5-6, 1-10)
  struct gkyl_comm *zcomm;   // communicator object for z smoothing

  struct gkyl_job_pool *job_pool; // Job pool  
  // arrays for global charge density and global smoothed (in z) charge density
  struct gkyl_array *rho_c_global_dg;
  struct gkyl_array *rho_c_global_smooth; 
  struct gkyl_array *phi; // arrays for updates

  struct gkyl_array *phi_host;  // host copy for use IO and initialization

  struct gkyl_range crossz_sub_range; // sub range of intersection of global range and local range
                                      // for solving subset of Poisson solves with parallelization in z

  struct gkyl_array *epsilon;
  struct gkyl_fem_parproj *fem_parproj; // FEM smoother for projecting DG functions onto continuous FEM basis
                                        // weight*phi_{fem} = phi_{dg} 
  struct gkyl_deflated_fem_poisson *deflated_fem_poisson; // poisson solver which solves on lines in x or planes in xy
                                                          // - nabla . (epsilon * nabla phi) - kSq * phi = rho
};


struct gk_field_mb* gk_field_mb_new(struct gkyl_gk_mb *gk_mb, struct gkyl_gyrokinetic_mb_app *mb_app, struct gkyl_gyrokinetic_app *app, int bidx);

void gk_field_mb_rhs(gkyl_gyrokinetic_mb_app *mb_app, struct gk_field_mb *field_mb, struct gkyl_gyrokinetic_app *app);

void
gk_field_mb_release(const gkyl_gyrokinetic_mb_app* mb_app, struct gk_field_mb *f);