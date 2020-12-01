#pragma once

#include <gkyl_array.h>

/** Generic accumulate macro */
#define gkyl_array_accumulate(out, a, inp) \
    _Generic((a),                          \
      float: gkyl_array_accumulate_flt,    \
      double: gkyl_array_accumulate_dbl)   \
    (out, a, inp)

/** Generic set macro */
#define gkyl_array_set(out, a, inp) \
    _Generic((a),                   \
      float: gkyl_array_set_flt,    \
      double: gkyl_array_set_dbl)   \
    (out, a, inp)

/**
 * Compute out = out + a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiple input array
 * @param inp Input array
 */
struct gkyl_array* gkyl_array_accumulate_dbl(struct gkyl_array *out,
  double a, const struct gkyl_array *inp);

struct gkyl_array* gkyl_array_accumulate_flt(struct gkyl_array *out,
  float a, const struct gkyl_array *inp);

/**
 * Set out = a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiple input array
 * @param inp Input array
 */
struct gkyl_array* gkyl_array_set_dbl(struct gkyl_array *out,
  double a, const struct gkyl_array *inp);

struct gkyl_array* gkyl_array_set_flt(struct gkyl_array *out,
  double a, const struct gkyl_array *inp);


