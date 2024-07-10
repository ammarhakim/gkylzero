#pragma once

/**

  Format description for raw Gkeyll output file.

  IF THIS FORMAT IF MODIFIED, PLEASE COPY AND THEN CHANGE THE
  DESCRIPTION SO WE HAVE THE OLDER VERSIONS DOCUMENTED HERE. UPDATE
  VERSION BY 1 EACH TIME YOU CHANGE THE FORMAT.

  The format of the gkyl binary output is as follows.

  ## Version 0: Jan 2021. Created by A.H. Note Version 0 has no header
     information

  Data      Type and meaning
  --------------------------
  ndim      uint64_t Dimension of field
  cells     uint64_t[ndim] number of cells in each direction
  lower     float64[ndim] Lower bounds of grid
  upper     float64[ndim] Upper bounds of grid
  esznc     uint64_t Element-size * number of components in field
  size      uint64_t Total number of cells in field
  DATA      size*esznc bytes of data

  ## Version 1: May 9th 2022. Created by A.H

  Data      Type and meaning
  --------------------------
  gkyl0     5 bytes
  version   uint64_t
  file_type uint64_t (See header gkyl_elem_type.h for file types)
  meta_size uint64_t Number of bytes of meta-data
  DATA      meta_size bytes of data. This is in msgpack format

  * For file_type = 1 (field) the above header is followed by

  real_type uint64_t. Indicates real type of data
  ndim      uint64_t Dimension of field
  cells     uint64_t[ndim] number of cells in each direction
  lower     float64[ndim] Lower bounds of grid
  upper     float64[ndim] Upper bounds of grid
  esznc     uint64_t Element-size * number of components in field
  size      uint64_t Total number of cells in field
  DATA      size*esznc bytes of data

  * For file_type = 2 (dynvec) the above header is followed by

  real_type uint64_t. Indicates real type of data
  esznc     uint64_t Element-size * number of components in field
  size      uint64_t Total number of cells in field
  TIME_DATA float64[size] bytes of data
  DATA      size*esznc bytes of data

  * For file_type = 3 (multi-range field) the above header is followed by

  real_type uint64_t. Indicates real type of data
  ndim      uint64_t Dimension of field
  cells     uint64_t[ndim] number of cells in each direction
  lower     float64[ndim] Lower bounds of grid
  upper     float64[ndim] Upper bounds of grid
  esznc     uint64_t Element-size * number of components in field
  size      uint64_t Total number of cells in field
  nrange    uint64_t Number of ranges stored in this file

  For each of the nrange ranges in the field the following data is
  present

  loidx     uint64_t[ndim] Index of lower-left corner of the range
  upidx     uint64_t[ndim] Index of upper-right corner of the range
  size      uint64_t Total number of cells in range
  DATA      size*esznc bytes of data

  Note: the global range in Gkeyll, of which each range is a part,
  is 1-indexed.

  * For file_type = 4 (block topology) there is no additional data 

 */

#include <stdio.h>

// The following utility functions allow to determine the size in
// bytes of the headers for gkyl output files.

size_t gkyl_base_hdr_size(size_t meta_sz);
size_t gkyl_file_type_1_hrd_size(int ndim);
size_t gkyl_file_type_1_partial_hrd_size(int ndim);
size_t gkyl_file_type_2_hrd_size(void);
size_t gkyl_file_type_3_hrd_size(int ndim);
size_t gkyl_file_type_3_range_hrd_size(int ndim);
