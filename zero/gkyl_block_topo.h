#pragma once

#include <gkyl_util.h>

// edge (lower/upper) with orientation (positive/negative)
enum gkyl_oriented_edge {
  GKYL_LOWER_POSITIVE,
  GKYL_LOWER_NEGATIVE,
  GKYL_UPPER_POSITIVE,
  GKYL_UPPER_NEGATIVE,
  GKYL_PHYSICAL, // edge on physical domain boundary
};

// Structrue to describe a single connection to a given target block
struct gkyl_target_edge {
    int bid; // block ID of target edge
    int dir; // logical direction of target edge
    enum gkyl_oriented_edge edge; // edge+orientaton of target edge
};

// Connections of a block in each direction and along each edge
struct gkyl_block_connections {
    struct gkyl_target_edge connections[GKYL_MAX_CDIM][2];
};

// Connectivity information for all blocks, describing the topology of
// the blocks
struct gkyl_block_topo {
    int num_blocks; // total number of blocks
    struct gkyl_block_connections *conn; // connection info for each block
};

/**
 * Construct a new empty block topology with given total number of
 * blocks.
 *
 * @param nblocks Total number of blocks in topology
 * @return New block topology
 */
struct gkyl_block_topo* gkyl_block_topo_new(int nblocks);

/**
 * Free block topology.
 */
void gkyl_block_topo_free(struct gkyl_block_topo* btopo);
