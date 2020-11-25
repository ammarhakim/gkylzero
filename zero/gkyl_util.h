#ifndef GKYL_UTIL_H
#define GKYL_UTIL_H

#include <stddef.h>

// This funny looking macro allows getting a pointer to the 'type'
// struct that contains an object 'member' given the 'ptr' to the
// 'member' inside 'type'. See https://en.wikipedia.org/wiki/Offsetof
#define container_of(ptr, type, member)                                 \
    ((type *)((char *)(1 ? (ptr) : &((type *)0)->member) - offsetof(type, member)))

// Select type-specific compare function
#define gkyl_compare(a, b, eps)                 \
    _Generic((a),                               \
      float: gkyl_compare_float,                \
      double: gkyl_compare_double)              \
    (a,b,eps)

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
int gkyl_compare_float(float a, float b, float eps);

/**
 * Compares two double numbers 'a' and 'b' to check if they are
 * sufficiently close by, where 'eps' is the relative tolerance.
 */
int gkyl_compare_double(double a, double b, double eps);

#endif // GKYL_UTIL_H
