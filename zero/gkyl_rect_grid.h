#pragma once

#include <stdio.h>
#include <gkyl_util.h>

/**
 * Rectangular grid object.
 */
struct gkyl_rect_grid {
  int ndim; // number of dimensions
  double lower[GKYL_MAX_DIM]; // lower-left corner
  double upper[GKYL_MAX_DIM]; // upper-right corner
  int cells[GKYL_MAX_DIM]; // number of cells    
  double dx[GKYL_MAX_DIM]; // cell spacing
  double cellVolume; // cell volume
};

/**
 * Create new grid object.
 *
 * @param grid Grid object to initialize.
 * @param ndim Dimension of grid
 * @param lower Coordinates of lower-left corner of grid
 * @param upper Coordinates of upper-right corner of grid
 * @param cells Number of cells in each direction
 */
void gkyl_rect_grid_init(struct gkyl_rect_grid *grid, int ndim,
  const double *lower, const double *upper, const int *cells);

/**
 * Get cell-center coordinates. idx is the zero-based cell index.
 *
 * @param grid Grid object
 * @param idx Index of cell (lower-left corner has all index (0,0,...) )
 * @param xc On output, cell-center coordinates of cell 'idx'
 */
void gkyl_rect_grid_cell_center(const struct gkyl_rect_grid *grid,
  const int *idx, double *xc);

/**
 * Write grid data to file. File must be opened by caller of this
 * function. Data is written in binary format.
 *
 * @param grid Grid object to write
 * @param fp File handle to write to.
 */
void gkyl_rect_grid_write(const struct gkyl_rect_grid *grid, FILE *fp);
