#pragma once

#include <gkyl_alloc.h>
#include <gkyl_rect_decomp.h>

struct gkyl_flash_decomp {
    int ndim; // Dimensionality.
    int total_ranks; // Total number of ranks in decomposition.
    int nblocks; // Number of blocks in decomposition.
}; 

/**
 * Initialize a new decomposition from given @a total_ranks,
 * number of blocks and the ranks needed for each block.
 * 
 * @param total_ranks Total number of available ranks.
 * @param nblocks Total number of blocks in system.
 * @param branks Array of size @a nblocks, with ranks in each block.
 * @return const struct gkyl_flash_decomp* 
 */
const struct gkyl_flash_decomp*
gkyl_flash_decomp_new(int total_ranks, int nblocks, int ndim, const int *branks);


/**
 * Get number of ranks in block @a bn.
 *
 * @param decomp Decomp object
 * @param bn Block nummber for which to fetch rank list
 * @return Number of ranks in block @a bn
 */
int gkyl_flash_decomp_nranks(const struct gkyl_flash_decomp *decomp, int bn);

/**
 * Get the list of ranks in block @a bn. The ranks are returned in @a
 * ranks. This must be big enough to store the list.
 *
 * @param decomp Decomp object
 * @param bn Block nummber for which to fetch rank list
 * @param ranks List of ranks in @a bn block
 */
void gkyl_flash_decomp_getranks(const struct gkyl_flash_decomp *decomp, int bn, int ranks[]);

/**
 * Release decomposition.
 *
 * @param decomp Decomp to release
 */
void gkyl_flash_decomp_release(const struct gkyl_flash_decomp *decomp);
