#pragma once

#include <stdio.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

/**
 * Write array data to file. Data is written as a binary file.
 *
 * @param arr Array object to write
 * @param fp File handle to write to.
 */
void gkyl_array_write(const struct gkyl_array *arr, FILE *fp);

/**
 * Write part of the array data to file. Data is written as a binary
 * file. The region of the array to write is specified in the range
 * object. This method will fail if range volume is greater than array
 * size.
 *
 * @param range Range describing portion of the array to output.
 * @param arr Array object to write
 * @param fp File handle to write to.
 */
void gkyl_sub_array_write(const struct gkyl_range *range,
  const struct gkyl_array *arr, FILE *fp);

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
int gkyl_grid_array_write(const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  const struct gkyl_array *arr, const char* fname);

/**
 * Print range information to file object.
 *
 * @param range Range object to print
 * @param nm Name of range
 * @param fp File object to print range information
 */
void gkyl_print_range(const struct gkyl_range* range, const char *nm, FILE *fp);
