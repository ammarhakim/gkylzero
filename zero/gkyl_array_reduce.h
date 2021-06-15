#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>

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
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_max_cu(double *out_d, const struct gkyl_array* inp, struct gkyl_range range);
