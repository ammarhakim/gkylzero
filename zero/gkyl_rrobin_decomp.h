#pragma once

// Round-robin decomposition
struct gkyl_rrobin_decomp {
  int total_ranks; // total number of ranks in decomp
  int nblocks; // number of blocks in decomposition
};  

/**
 * Create a new round-robin decomposition with given @a total_ranks,
 * number of "blocks" and the ranks needed for each block.
 *
 * @param total_ranks Total number of ranks available
 * @param nblocks Number of "blocks" in system
 * @param branks An array of @a nblocks size, each with ranks in each block
 * @return New round-robin decomposition
 */
const struct gkyl_rrobin_decomp* gkyl_rrobin_decomp_new(int total_ranks, int nblocks,
  const int *branks);

/**
 * Release round-robin decomposition.
 *
 * @param rr Round-robin decom to release
 */
void gkyl_rrobin_decomp_release(const struct gkyl_rrobin_decomp *rr);
