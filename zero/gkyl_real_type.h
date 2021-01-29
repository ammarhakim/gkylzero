#pragma once

// Supported real types
#define GKYL_FLOAT 1
#define  GKYL_DOUBLE 2

#ifdef GKYL_REAL_TYPE
// choose real type based on user-specified flag

# if GKYL_REAL_TYPE == 1

#define GKYL_REAL_TYPE_IS_FLOAT
typedef float gkyl_real;
static const int gkyl_real_type_id = 1;
static const gkyl_real gkyl_real_type_zero = 0.0f;

# elif GKYL_REAL_TYPE == 2

#define GKYL_REAL_TYPE_IS_DOUBLE

typedef double gkyl_real;
static const int gkyl_real_type_id = 2;
static const gkyl_real gkyl_real_type_zero = 0.0;

# else
#  error "GKYL_REAL_TYPE should be 1 (float) or 2 (double)"
# endif

#else

// default real type is double
#define GKYL_REAL_TYPE_IS_DOUBLE

typedef double gkyl_real;
static const int gkyl_real_type_id = 2;
static const gkyl_real gkyl_real_type_zero = 0.0;

#endif
