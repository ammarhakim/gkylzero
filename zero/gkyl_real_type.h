#pragma once

// Supported real types
enum gkyl_real_type { GKYL_FLOAT = 1, GKYL_DOUBLE };

#ifdef GKYL_REAL_TYPE
// choose real type based on user-specified flag
# if GKYL_REAL_TYPE == GKYL_FLOAT
typedef float gkyl_real;
# elif GKYL_REAL_TYPE == GKYL_DOUBLE
typedef double gkyl_real;
# else
#  error "GKYL_REAL_TYPE should be 1 (float) or 2 (double)"
# endif

#else
// default real type is double
typedef double gkyl_real;

#endif
