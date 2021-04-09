#pragma once

#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

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
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range);
