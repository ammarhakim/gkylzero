// Private header file for multi-block gyrokinetic solver.
#pragma once

#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gyrokinetic_multib.h>

// Collection of a buffer and send/recv states used in sync.
struct gkyl_gyrokinetic_multib_buff_state {
  gkyl_mem_buff send_buff, recv_buff;
  struct gkyl_comm_state *send_state, *recv_state;
};

struct gkyl_gyrokinetic_multib_sync_data {
  struct gkyl_gyrokinetic_multib_buff_state lower[GKYL_MAX_CDIM], upper[GKYL_MAX_CDIM];
};

// Top-level multiblock gyrokinetic App object.
struct gkyl_gyrokinetic_multib_app {
  char name[128]; // Name of app.
  struct gkyl_job_pool *job_pool; // Job pool.

  int cdim, vdim; // Conf, velocity space dimensions.
  int poly_order; // Polynomial order.
  double tcurr; // Current time.
  double cfl; // CFL number.

  int num_species; // Number of charged species.
  int num_neut_species; // Number of neutral species.

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

  // Buffers and statuses used in syncing, 1 per boundary, per block.
  struct gkyl_gyrokinetic_multib_sync_data *syncdat;

  struct gkyl_gyrokinetic_stat stat; // statistics
};

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
