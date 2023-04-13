#pragma once

#include <stdio.h>
#include <stdint.h>

// Type of element stored in array
enum gkyl_elem_type { GKYL_INT, GKYL_FLOAT, GKYL_DOUBLE, GKYL_USER };

// Array reduce operators
enum gkyl_array_op { GKYL_MIN, GKYL_MAX, GKYL_SUM };

// code for array datatype for use in IO
static const uint64_t gkyl_array_data_type[] = {
  [GKYL_INT] = 0,
  [GKYL_FLOAT] = 1,
  [GKYL_DOUBLE] = 2,
  [GKYL_USER] = 32,
};

// size in bytes for various data-types
static const size_t gkyl_elem_type_size[] = {
  [GKYL_INT] = sizeof(int),
  [GKYL_FLOAT] = sizeof(float),
  [GKYL_DOUBLE] = sizeof(double),
  [GKYL_USER] = 1,
};

// file types for raw IO of Gkeyll data
enum gkyl_file_type {
  GKYL_FIELD_DATA_FILE,
  GKYL_DYNVEC_DATA_FILE,
  GKYL_MULTI_RANGE_DATA_FILE
};

static const uint64_t gkyl_file_type_int[] = {
  [GKYL_FIELD_DATA_FILE] = 1,
  [GKYL_DYNVEC_DATA_FILE] = 2,
  [GKYL_MULTI_RANGE_DATA_FILE] = 3
};
