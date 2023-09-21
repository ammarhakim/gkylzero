#pragma once

#include <stdio.h>
#include <stdbool.h>

#include <gkyl_util.h>
#include <gkyl_range.h>

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
 * Find cell indecies of point
 *
 * @param grid Grid object
 * @param point The point to find the cell incecies at
 * @param pickLower
 * @param knownIdx Any known indecies of where the point is
 * @return Pointer to cell indecies 
 */
int* gkyl_find_cell(struct gkyl_rect_grid *grid, const double *point, bool pickLower,const int **knownIdx);

/**
 * Get cell-center coordinates. Note that idx is a 1-based cell index,
 * i.e. the lower-left corner is (1,1,...).
 *
 * @param grid Grid object
 * @param idx Index of cell (lower-left corner has all index (1,1,...) )
 * @param xc On output, cell-center coordinates of cell 'idx'
 */
GKYL_CU_DH
static inline void
gkyl_rect_grid_cell_center(const struct gkyl_rect_grid *grid,
  const int *idx, double *xc)
{
  for (int i=0; i<grid->ndim; ++i)
    xc[i] = grid->lower[i]+(idx[i]-0.5)*grid->dx[i];
}

/**
 * Get index extents in direction @a dir. The extents are inclusive.
 *
 * @param grid Grid object
 * @param dir Direction in which to get extents
 * @param ext On output, inclusive extents in direction @a dir.
 */
GKYL_CU_DH
static inline void
gkyl_rect_grid_extents(const struct gkyl_rect_grid *grid, int dir, int ext[2])
{
  ext[0] = 1; ext[1] = grid->cells[dir];
}

/**
 * Write grid data to file. File must be opened by caller of this
 * function. Data is written in binary format.
 *
 * @param grid Grid object to write
 * @param fp File handle to write to.
 */
void gkyl_rect_grid_write(const struct gkyl_rect_grid *grid, FILE *fp);

/**
 * Read data from file and initialize grid. File must be opened by
 * caller of this function.
 *
 * @param grid Grid object to initialize
 * @param fp File handle to read fromx.
 * @return True if read succeeded, false otherwise
 */
bool gkyl_rect_grid_read(struct gkyl_rect_grid *grid, FILE *fp);
