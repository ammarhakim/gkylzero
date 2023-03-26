#pragma once

#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

struct gkyl_rect_decomp {
  int ndim; // dimension of decomposition
  int ndecomp; // number of sub-domains
  struct gkyl_range *ranges; // decomposed ranges
};

/**
 * Create a new decomposition given "cuts" in each direction. The
 * total number of decomposed ranges are product of all cuts.
 *
 * @param ndim Number of dimensions
 * @param cuts Cuts in each direction.
 * @param range Range to decompose
 */
struct gkyl_rect_decomp *gkyl_rect_decomp_new_from_cuts(int ndim,
  const int cuts[], const struct gkyl_range *range);

/**
 * Free decomposition.
 *
 * @param decomp Decomposition to free
 */
void gkyl_rect_decomp_release(struct gkyl_rect_decomp *decomp);

/**
 * Create range and extended ranges from grid and ghost-cell data. The
 * range is a sub-range of the extended range.
 *
 * @param grid Grid to compute ranges for
 * @param nghost Number of ghost-cells in each direction
 * @param ext_range On output, extended range spanning grid+ghost-cells
 * @param range On output, range spanning grid. Sub-range of ext_range.
 */
void gkyl_create_grid_ranges(const struct gkyl_rect_grid *grid,
  const int *nghost, struct gkyl_range *ext_range,
  struct gkyl_range *range);

/**
 * Create range and extended ranges from give range and ghost-cell
 * data. The range is a sub-range of the extended range.
 *
 * @param inrange Input range to use
 * @param nghost Number of ghost-cells in each direction
 * @param ext_range On output, extended range spanning inrange+ghost-cells
 * @param range On output, range same as inrange, but sub-range of ext_range.
 */
void gkyl_create_ranges(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range);
