#pragma once

#include <stddef.h>
#include <time.h>

// Maximum dimensions supported
#ifndef GKYL_MAX_DIM
# define GKYL_MAX_DIM 7
#endif

// Default alignment boundary (seems agressive or maybe not enough?)
#ifndef GKYL_DEF_ALIGN
# define GKYL_DEF_ALIGN 64
#endif

// This funny looking macro allows getting a pointer to the 'type'
// struct that contains an object 'member' given the 'ptr' to the
// 'member' inside 'type'. See https://en.wikipedia.org/wiki/Offsetof
#define container_of(ptr, type, member)                                 \
    ((type *)((char *)(1 ? (ptr) : &((type *)0)->member) - offsetof(type, member)))

// Select type-specific compare function
#define gkyl_compare(a, b, eps)                 \
    _Generic((a),                               \
      float: gkyl_compare_flt,                  \
      double: gkyl_compare_dbl)                 \
    (a, b, eps)

/**
 * Print error message to stderr and exit.
 *
 * @param msg Error message.
 */
void gkyl_exit(const char* msg);

/**
 * Compares two float numbers 'a' and 'b' to check if they are
 * sufficiently close by, where 'eps' is the relative tolerance.
 */
int gkyl_compare_flt(float a, float b, float eps);

/**
 * Compares two double numbers 'a' and 'b' to check if they are
 * sufficiently close by, where 'eps' is the relative tolerance.
 */
int gkyl_compare_dbl(double a, double b, double eps);

/**
 * Gets wall-clock time in secs/nanoseconds.
 * 
 * @return Status flag.
 */
int gkyl_wall_clock(struct timespec *ts);

/**
 * Difference between two timespec objects.
 *
 * @param tstart Start time
 * @param tend End time 
 * @return Time object representing difference
 */
struct timespec gkyl_time_diff(struct timespec *tstart, struct timespec *tend);

/**
 * Compute in secs time stored in timespec object.
 *
 * @param tm Timespec object
 * @return Time in seconds
 */
double gkyl_time_sec(struct timespec tm);
