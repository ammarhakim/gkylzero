#pragma once

#include <gkyl_array.h>

/** Generic clear macro */
#define gkyl_array_clear(out, a)       \
    _Generic((a),                      \
      float: gkyl_array_clear_float,   \
      double: gkyl_array_clear_double) \
    (out, a)

/** Generic accumulate macro */
#define gkyl_array_accumulate(out, a, inp)      \
    _Generic((a),                               \
      float: gkyl_array_accumulate_float,       \
      double: gkyl_array_accumulate_double)     \
    (out, a, inp)

/** Generic set macro */
#define gkyl_array_set(out, a, inp)             \
    _Generic((a),                               \
      float: gkyl_array_set_float,              \
      double: gkyl_array_set_double)            \
    (out, a, inp)

/** Generic scale macro (uses set) */
#define gkyl_array_scale(out, a)                \
    _Generic((a),                               \
      float: gkyl_array_set_float,              \
      double: gkyl_array_set_double)            \
    (out, a, out)

/** Generic uniop macro */
#define gkyl_array_uniop(op, a, out, b, inp)    \
    _Generic((a),                               \
      float: gkyl_array_uniop_float,            \
      double: gkyl_array_uniop_double)          \
    (op, a, out, b, inp)

/**
 * Clear out = val. Returns out.
 *
 * @param out Output array
 * @param val Factor to set 
 * @return out array
 */
struct gkyl_array* gkyl_array_clear_double(struct gkyl_array *out, double val);

struct gkyl_array* gkyl_array_clear_float(struct gkyl_array *out, float val);

/**
 * Compute out = out + a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 */
struct gkyl_array* gkyl_array_accumulate_double(struct gkyl_array *out,
  double a, const struct gkyl_array *inp);

struct gkyl_array* gkyl_array_accumulate_float(struct gkyl_array *out,
  float a, const struct gkyl_array *inp);

/**
 * Set out = a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 */
struct gkyl_array* gkyl_array_set_double(struct gkyl_array *out,
  double a, const struct gkyl_array *inp);

struct gkyl_array* gkyl_array_set_float(struct gkyl_array *out,
  float a, const struct gkyl_array *inp);

/**
 * Generic unary operator. Computes out = a*out + b*OP[inp], where
 * "OP" represents an operator that takes a single scalar input.
 *
 * @param op Operator to apply
 * @param a Multiplicative factor on out
 * @param out Output array
 * @param b Multiplicative factor on inp
 * @param inp Input array
 * @return out array
 */
struct gkyl_array* gkyl_array_uniop_double(const char *op, double a,
  struct gkyl_array *out, double b, const struct gkyl_array *inp);

struct gkyl_array* gkyl_array_uniop_float(const char *op, float a,
  struct gkyl_array *out, float b, const struct gkyl_array *inp);
