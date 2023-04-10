#pragma once

#include <stdio.h>

// Type of element stored in array
enum gkyl_elem_type { GKYL_INT, GKYL_FLOAT, GKYL_DOUBLE, GKYL_USER };

// Array reduce operators
enum gkyl_array_op { GKYL_MIN, GKYL_MAX, GKYL_SUM };

// size in bytes for various data-types
static const size_t gkyl_elem_type_size[] = {
  [GKYL_INT] = sizeof(int),
  [GKYL_FLOAT] = sizeof(float),
  [GKYL_DOUBLE] = sizeof(double),
  [GKYL_USER] = 1,
};
