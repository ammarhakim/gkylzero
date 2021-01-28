#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>

/** Generic clear macro */
#define gkyl_array_clear(out, a)       \
    _Generic((a),                      \
      float: gkyl_array_clear_float,   \
      double: gkyl_array_clear_double) \
    (out, a)

/** Generic clear macro */
#define gkyl_array_clear_range(out, a, range)   \
    _Generic((a),                               \
      float: gkyl_array_clear_range_float,      \
      double: gkyl_array_clear_range_double)    \
    (out, a, range)

/** Generic accumulate macro */
#define gkyl_array_accumulate(out, a, inp)      \
    _Generic((a),                               \
      float: gkyl_array_accumulate_float,       \
      double: gkyl_array_accumulate_double)     \
    (out, a, inp)

/** Generic accumulate macro */
#define gkyl_array_accumulate_range(out, a, inp, range) \
    _Generic((a),                                       \
      float: gkyl_array_accumulate_range_float,         \
      double: gkyl_array_accumulate_range_double)       \
    (out, a, inp, range)

/** Generic set macro */
#define gkyl_array_set(out, a, inp)             \
    _Generic((a),                               \
      float: gkyl_array_set_float,              \
      double: gkyl_array_set_double)            \
    (out, a, inp)

/** Generic set macro */
#define gkyl_array_set_range(out, a, inp, range)        \
    _Generic((a),                                       \
      float: gkyl_array_set_range_float,                \
      double: gkyl_array_set_range_double)              \
    (out, a, inp, range)

/** Generic scale macro (uses set) */
#define gkyl_array_scale(out, a)                \
    _Generic((a),                               \
      float: gkyl_array_set_float,              \
      double: gkyl_array_set_double)            \
    (out, a, out)

/** Generic scale macro (uses set) */
#define gkyl_array_scale_range(out, a, range)   \
    _Generic((a),                               \
      float: gkyl_array_set_range_float,        \
      double: gkyl_array_set_range_double)      \
    (out, a, out, range)

/** Generic reduce macro */
#define gkyl_array_reduce(arr, op, out)                      \
    _Generic((out),                                          \
      float *: gkyl_array_reduce_float,                      \
      double *: gkyl_array_reduce_double)                    \
    (arr, op, out)

/** Generic reduce macro */
#define gkyl_array_reduce_range(res, arr, op, range)         \
    _Generic((res),                                          \
      float *: gkyl_array_reduce_range_float,                \
      double *: gkyl_array_reduce_range_double)              \
    (res, arr, op, range)

// Array reduce operators
enum gkyl_array_op { GKYL_MIN, GKYL_MAX };

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
 * Clear out = val. Returns out.
 *
 * @param out Output array
 * @param val Factor to set 
 * @return out array
 */
struct gkyl_array* gkyl_array_clear_range_double(struct gkyl_array *out, double val,
  const struct gkyl_range *range);
struct gkyl_array* gkyl_array_clear_range_float(struct gkyl_array *out, float val,
  const struct gkyl_range *range);

/**
 * Compute out = out + a*inp over a range of indices.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @param range Range specifying region to accumulate
 * @return out array
 */
struct gkyl_array* gkyl_array_accumulate_range_double(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, const struct gkyl_range *range);
struct gkyl_array* gkyl_array_accumulate_range_float(struct gkyl_array *out,
  float a, const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Set out = a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 * @param range Range specifying region to set
 */
struct gkyl_array* gkyl_array_set_range_double(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, const struct gkyl_range *range);
struct gkyl_array* gkyl_array_set_range_float(struct gkyl_array *out,
  float a, const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Copy out inp. Returns out.
 *
 * @param out Output array
 * @param inp Input array
 * @return out array
 */
struct gkyl_array* gkyl_array_copy_range(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Perform an "reduce" operation of data in the array.
 *
 * @param arr Array to perform reduction on.
 * @param op Reduction operators
 * @return Reduced result
 */
double gkyl_array_reduce_double(const struct gkyl_array *arr, enum gkyl_array_op op, double *out);
float gkyl_array_reduce_float(const struct gkyl_array *arr, enum gkyl_array_op op, float *out);

/**
 * Perform an "reduce" operation of data in the array. Data is reduced
 * component-wise.
 *
 * @param res On output, reduced values
 * @param arr Array to perform reduction on.
 * @param op Reduction operators
 * @param range Range specifying region
 */
void gkyl_array_reduce_range_double(double *res,
  const struct gkyl_array *arr, enum gkyl_array_op op, const struct gkyl_range *range);
void gkyl_array_reduce_range_float(float *res,
  const struct gkyl_array *arr, enum gkyl_array_op op, const struct gkyl_range *range);

/**
 * Copy region of array into a buffer. The buffer must be preallocated
 * and at least of size arr->size*arr->elemSz bytes.
 *
 * @param arr Array to copy from
 * @param data Output data buffer.
 * @param range Range specifying region to copy
 */
void gkyl_array_copy_to_buffer(void *data, const struct gkyl_array *arr,
  const struct gkyl_range *range);

/**
 * Copy buffer into region of array. The array must be preallocated.
 *
 * @param arr Array to copy from
 * @param data Output data buffer.
 * @param range Range specifying region to copy
 */
void gkyl_array_copy_from_buffer(struct gkyl_array *arr,
  const void *data, const struct gkyl_range *range);
