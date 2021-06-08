#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>

// Array reduce operators
enum gkyl_array_op { GKYL_MIN, GKYL_MAX };

/**
 * Clear out = val. Returns out.
 *
 * @param out Output array
 * @param val Factor to set 
 * @return out array
 */
struct gkyl_array* gkyl_array_clear(struct gkyl_array *out, double val);

/**
 * Compute out = out + a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 */
struct gkyl_array* gkyl_array_accumulate(struct gkyl_array *out,
  double a, const struct gkyl_array *inp);

/**
 * Set out = a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 */
struct gkyl_array* gkyl_array_set(struct gkyl_array *out,
  double a, const struct gkyl_array *inp);

/**
 * Scale out = a*out. Returns out.
 *
 * @param out Output array
 * @param a Factor to scale
 * @return out array
 */
struct gkyl_array* gkyl_array_scale(struct gkyl_array *out, double a);

/**
 * Clear out = val. Returns out.
 *
 * @param out Output array
 * @param val Factor to set 
 * @return out array
 */
struct gkyl_array* gkyl_array_clear_range(struct gkyl_array *out, double val,
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
struct gkyl_array* gkyl_array_accumulate_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Set out = a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 * @param range Range specifying region to set
 */
struct gkyl_array* gkyl_array_set_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Scale out = a*ut. Returns out.
 *
 * @param out Output array
 * @param a Factor to scale by
 * @return out array
 * @param range Range specifying region to scale
 */
struct gkyl_array* gkyl_array_scale_range(struct gkyl_array *out,
  double a, const struct gkyl_range *range);

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
double gkyl_array_reduce(const struct gkyl_array *arr, enum gkyl_array_op op, double *out);

/**
 * Perform an "reduce" operation of data in the array. Data is reduced
 * component-wise.
 *
 * @param res On output, reduced values
 * @param arr Array to perform reduction on.
 * @param op Reduction operators
 * @param range Range specifying region
 */
void gkyl_array_reduce_range(double *res,
  const struct gkyl_array *arr, enum gkyl_array_op op, const struct gkyl_range *range);

/**
 * Copy region of array into a buffer. The buffer must be preallocated
 * and at least of size arr->size*arr->elemSz bytes.
 *
 * @param arr Array to copy from
 * @param data Output data buffer.
 * @param range Range specifying region to copy from
 */
void gkyl_array_copy_to_buffer(void *data, const struct gkyl_array *arr,
  const struct gkyl_range *range);

/**
 * Copy buffer into region of array. The array must be preallocated.
 *
 * @param arr Array to copy into
 * @param data Output data buffer.
 * @param range Range specifying region to copy into
 */
void gkyl_array_copy_from_buffer(struct gkyl_array *arr,
  const void *data, const struct gkyl_range *range);

/**
 * Host-side wrappers for array operations
 */
void
gkyl_array_clear_cu(struct gkyl_array* out, double val);

void
gkyl_array_accumulate_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp);

void
gkyl_array_set_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp);

void
gkyl_array_scale_cu(struct gkyl_array* out, double a);

/**
 * Host-side wrappers for range-based array operations
 */
void
gkyl_array_clear_range_cu(int numBlocks, int numThreads,
  struct gkyl_array *out, double val, const struct gkyl_range* range);

void
gkyl_array_accumulate_range_cu(int numBlocks, int numThreads, struct gkyl_array *out,
  double a, const struct gkyl_array* inp, const struct gkyl_range* range);

void
gkyl_array_set_range_cu(int numBlocks, int numThreads, struct gkyl_array *out,
  double a, const struct gkyl_array* inp, const struct gkyl_range* range);

void
gkyl_array_scale_range_cu(int numBlocks, int numThreads, struct gkyl_array *out,
  double a, const struct gkyl_range* range);

void
gkyl_array_copy_range_cu(int numBlocks, int numThreads, struct gkyl_array *out,
  const struct gkyl_array* inp, const struct gkyl_range* range);
