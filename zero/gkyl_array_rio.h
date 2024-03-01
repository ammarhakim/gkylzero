#pragma once

#include <stdio.h>

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

/**
 * Write out grid and array data to file in .gkyl format so postgkyl
 * can understand it. 
 * 
 * @param grid Grid object to write
 * @param range Range describing portion of the array to output.
 * @param arr Array object to write
 * @param fname Name of output file (include .gkyl extension)
 * @return Status flag: 0 if write succeeded, 'errno' otherwise
 */
int gkyl_grid_sub_array_write(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_array *arr, const char *fname);

/**
 * Read grid and array data from file. The input array must be
 * pre-allocated and must be big enough to hold the read data.
 * 
 * @param grid Grid object to read
 * @param range Range describing portion of the array.
 * @param arr Array object to read
 * @param fname Name of input file
 * @return Status flag: 0 if write succeeded, 'errno' otherwise
 */
int gkyl_grid_sub_array_read(struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  struct gkyl_array *arr, const char* fname);

/**
 * Read grid and array data from file, creating a new array.
 * 
 * @param grid On outout, grid on which array is defined.
 * @param fname Name of input file
 * @return Newly created array object
 */
struct gkyl_array *gkyl_grid_array_new_from_file(struct gkyl_rect_grid *grid,
  const char* fname);
