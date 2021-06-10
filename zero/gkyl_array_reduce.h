#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>

struct gkyl_array_reduce_util
{
  double initValue;
//  void *d_temp = NULL;
//  size_t temp_bytes = 0;
  void *d_temp;
  size_t temp_bytes;
  void *out_d;
};

/**
 * Initialize permanent objects needed in reductions (stored in a reduceObject struct).
 *
 * @param inp A gkyl_array (the same one to be reduced or one like it).
 * @param redOb A reduceObject struct with permanent reduction objects assigned here.
 */
void gkyl_array_reduce_init_cu(struct gkyl_array_reduce_util* redutil, struct gkyl_array* inp);

/**
 * Reduce a gkyl_array component-wise.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_max_cu(double *out_d, const struct gkyl_array* inp);

/**
 * Reduce a gkyl_array component-wise over a specified range.
 *
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_range_max_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range range);

/**
 * Free the memory allocated in reduce_init_cu.
 *
 * @param redOb A reduceObject struct with permanent reduction objects (released here).
 */
void gkyl_array_reduce_free_cu(struct gkyl_array_reduce_util* redutil);
