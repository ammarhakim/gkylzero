#pragma once

#include <gkyl_ref_count.h>
#include <gkyl_util.h>

// edge (lower/upper) with orientation (positive/negative)
enum gkyl_oriented_edge {
  GKYL_LOWER_POSITIVE = 1, // this needs to be 1 and not 0
  GKYL_LOWER_NEGATIVE,
  GKYL_UPPER_POSITIVE,
  GKYL_UPPER_NEGATIVE,
  GKYL_PHYSICAL, // edge on physical domain boundary
};

// Connection to a given target edge
struct gkyl_target_edge {
  int bid; // block ID of target edge
  int dir; // logical direction of target edge
  enum gkyl_oriented_edge edge; // edge+orientaton of target edge
};

// Connections of a block in each direction and along each edge (0:
// lower, 1:upper)
struct gkyl_block_connections {
  struct gkyl_target_edge connections[GKYL_MAX_CDIM][2];
};

// Connectivity information for all blocks, describing topology of
// blocks
struct gkyl_block_topo {
  int ndim; // dimension
  int num_blocks; // total number of blocks
  struct gkyl_block_connections *conn; // connection info for each block

  struct gkyl_ref_count ref_count;
};

/**
 * Construct a new empty N-dim block topology with given total number
 * of blocks. The returned block topology needs to be manually
 * constructed for it to be useful.
 *
 * @param ndim Dimension of space
 * @param nblocks Total number of blocks in topology
 * @return New block topology
 */
struct gkyl_block_topo* gkyl_block_topo_new(int ndim, int nblocks);

/**
 * Acquire pointer to block-topology. The pointer must be released
 * using release method.
 *
 * @param btopo Block-topo to which reference is required
 * @return Pointer to acquired block-topo
 */
struct gkyl_block_topo* gkyl_block_topo_acquire(const struct gkyl_block_topo* btopo);

/**
 * Check consistency of block topology: the topology typically has
 * redundant data in it. This method ensures the redundant data is
 * consistent. It also checks if there are any blocks with unspecified
 * connections.
 *
 * @param btopo Block topology to check
 * @return 1 if topology is consistent, 0 otherwise.
 */
int gkyl_block_topo_check_consistency(const struct gkyl_block_topo *btopo);

/**
 * Free block topology.
 *
 * @return Block topology to free
 */
void gkyl_block_topo_release(struct gkyl_block_topo* btopo);
