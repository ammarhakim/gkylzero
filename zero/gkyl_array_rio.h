#pragma once

#include <stdio.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Array header data to write: this is for low-level control and is
// typically not something most user would every encounter
struct gkyl_array_header_info {
  uint64_t file_type; // file type
  enum gkyl_elem_type etype; // element type
  uint64_t esznc; // elem sz * number of components
  uint64_t tot_cells; // total number of cells in grid
};

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
int gkyl_grid_sub_array_write(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_array *arr, const char *fname);

// Same as above method, except takes an open FILE pointer. @a fp must
// be opened with right permissions.
int gkyl_grid_sub_array_write_fp(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_array *arr, FILE *fp);

/**
 * Write out grid and array data header data to file. Note that only
 * HEADER is written and NOT the array data itself.
 *
 * @param grid Grid object to write
 * @param hrd Header data.
 * @return Status flag: 0 if write succeeded, 'errno' otherwise
 */
int gkyl_grid_sub_array_header_write_fp(const struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, FILE *fp);

/**
 * Read data from file and create new array. 
 *
 * @param type Type of data in array
 * @param fp File handle to read from.
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_new_from_file(enum gkyl_elem_type type, FILE *fp);

/**
 * Read part of the array from file. The region of the array to read
 * is specified in the range object. This method will fail if range
 * volume is greater than array size. The input array must be
 * pre-allocated and must be big enough to hold the read data.
 *
 * @param range Range describing portion of the array to read.
 * @param arr Array object to read into
 * @param fp File handle to read from.
 * @param True on successful read, false otherwise
 */
bool gkyl_sub_array_read(const struct gkyl_range *range,
  struct gkyl_array *arr, FILE *fp);

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
int gkyl_grid_sub_array_read(struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char* fname);

/**
 * Read grid and array data from file, creating a new array.
 * 
 * @param grid On outout, grid on which array is defined.
 * @param fname Name of input file
 * @return Newly created array object
 */
struct gkyl_array* gkyl_grid_array_new_from_file(struct gkyl_rect_grid *grid, const char* fname);
