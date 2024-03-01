#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Array header data to write: this is for low-level control and is
// typically not something most users would every encounter
struct gkyl_array_header_info {
  uint64_t file_type; // file type
  enum gkyl_elem_type etype; // element type
  uint64_t esznc; // elem sz * number of components
  uint64_t tot_cells; // total number of cells in grid
  uint64_t meta_size; // size in bytes of meta-data embedded in header
  char *meta;         // meta-data as byte array
  uint64_t nrange; // number of ranges
};

/**
 * Write out grid and array data header data to file. Note that only
 * HEADER is written and NOT the array data itself.
 *
 * @param grid Grid object to write
 * @param hrd Header data.
 * @return Status flag: 0 if write succeeded, 'errno' otherwise
 */
int gkyl_grid_sub_array_header_write_fp(const struct gkyl_rect_grid *grid,
  const struct gkyl_array_header_info *hdr, FILE *fp);

/**
 * Read grid and array data header data from file. Note that only
 * HEADER is read and NOT the array data itself.
 *
 * @param grid Grid object to read
 * @param hrd On output, Header data.
 * @return Status flag: 0 if read succeeded, 'errno' otherwise
 */
int gkyl_grid_sub_array_header_read_fp(struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, FILE *fp);

