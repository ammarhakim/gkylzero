// Private header file for multi-block gyrokinetic solver.
#pragma once

#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gyrokinetic_multib.h>

// Top-level multiblock gyrokinetic App object.
struct gkyl_gyrokinetic_multib_app {
  char name[128]; // Name of app.
  struct gkyl_job_pool *job_pool; // Job pool.

  int cdim, vdim; // Conf, velocity space dimensions.
  int poly_order; // Polynomial order.
  double tcurr; // Current time.
  double cfl; // CFL number.

  bool use_mpi; // Should we use MPI (if present).
  bool use_gpu; // Should we use GPU (if present).

  struct gkyl_comm *comm_multib;   // Multiblock communicator object.
  struct gkyl_comm **comm_intrab;
  struct gkyl_rect_decomp **decomp_intrab;

  struct gkyl_block_topo *btopo; // Block topology.

  int num_blocks, num_blocks_local; // Total and local number of blocks.
  int *block_idxs; // List of blocks handled on this rank.
  struct gkyl_gyrokinetic_app **blocks; // Pointers to blocks on this rank.
  int *ranks_per_block; // All ranks in a single block. Num ranks per block is
                        // the same as decomp_intrab[bidx]->ndecomp.
  int *cuts_vol_cum_per_block; // Cumulative num_ranks prior at given rank.

  struct gk_field_multib *field_multib; // mb field object.
  struct gkyl_gyrokinetic_stat stat; // statistics
};

// field data
struct gk_field_multib {
  struct gkyl_gyrokinetic_field info; // data for field

  enum gkyl_gkfield_id gkfield_id;

  // z ranges
  struct gkyl_range crossz, crossz_ext; // crossz, crossz-ext conf-space ranges. Cross-block ranges across z boundaries for the smoother. We need 5 of these in the double null case (2-3-4, 11-12, 1-8-9, 5-6, 1-10)
  int crossz_block_idxs[GKYL_MAX_BLOCKS];
  int numz_blocks;
  struct gkyl_range *cut_ranges; // All ranges connected along z
  int num_cuts; // number of cut ranges
  int *zranks; // List of ranks that own each of these ranges
  int *unique_zranks; // List of unique ranks connected along z
  struct gkyl_comm *zcomm;   // communicator object for z smoothing
  struct gkyl_rect_decomp *zdecomp; // decomp on which zcomm is based

  struct gkyl_job_pool *job_pool; // Job pool  
  // arrays for global charge density and global smoothed (in z) charge density
  struct gkyl_array *rho_c_global_dg;
  struct gkyl_array *rho_c_global_smooth; 
  struct gkyl_array *phi; // arrays for updates

  struct gkyl_array *phi_host;  // host copy for use IO and initialization

  struct gkyl_range crossz_sub_range[GKYL_MAX_BLOCKS]; // sub range of intersection of global range and local range
                                      // for solving subset of Poisson solves with parallelization in z

  struct gkyl_array *epsilon;
  struct gkyl_fem_parproj *fem_parproj; // FEM smoother for projecting DG functions onto continuous FEM basis
                                        // weight*phi_{fem} = phi_{dg} 
  struct gkyl_deflated_fem_poisson *deflated_fem_poisson; // poisson solver which solves on lines in x or planes in xy
                                                          // - nabla . (epsilon * nabla phi) - kSq * phi = rho
};


struct gk_field_multib* gk_field_multib_new(struct gkyl_gk_multib *gk_multib, struct gkyl_gyrokinetic_multib_app *mb_app);

/**
 * Get the number of ranks in a given block of ID @a bidx.
 *
 * @param mba MB gyrokinetic app.
 * @param bidx Block ID.
 *
 * @return number of ranks in that block.
 */
int
gkyl_gyrokinetic_multib_num_ranks_per_block(gkyl_gyrokinetic_multib_app *mba, int bidx);

/**
 * Get the list of rank IDs for @a ranks owning a given block of ID @a bidx.
 * Also returns the number of ranks in that block.
 *
 * @param mba MB gyrokinetic app.
 * @param bidx Block ID.
 * @param ranks list of ranks in that block.
 *
 * @return number of ranks in that block.
 */
int
gkyl_gyrokinetic_multib_ranks_per_block(gkyl_gyrokinetic_multib_app *mba, int bidx, int *ranks);

void gk_field_multib_rhs(gkyl_gyrokinetic_multib_app *mb_app, struct gk_field_multib *field_multib);

void
gk_field_multib_release(const gkyl_gyrokinetic_multib_app* mb_app, struct gk_field_multib *f);
