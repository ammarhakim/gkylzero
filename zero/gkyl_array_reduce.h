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
void gkyl_array_reduce_range_max_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range);

/**
 * Reduce a gkyl_array component-wise.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_min_cu(double *out_d, const struct gkyl_array* inp);

/**
 * Reduce a gkyl_array component-wise over a specified range.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_min_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range);

/**
 * Reduce a gkyl_array component-wise.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 */
void gkyl_array_reduce_sum_cu(double *out_d, const struct gkyl_array* inp);

/**
 * Reduce a gkyl_array component-wise over a specified range.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param range A gkyl_range over which to perform the reduction.
 */
void gkyl_array_reduce_range_sum_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range);

/**
 * Max reduce a gkyl_array evaluating the DG field in each cell at Gauss-Legendre nodes
 * and reducing those values. Allows for reducing a single vector component for arrays
 * with multiple scalar fields within it.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param comp Vector component to reduce.
 * @param basis Baisis DG coefficients expand in (device pointer).
 */
void gkyl_array_reducec_dg_max_cu(double *out_d, const struct gkyl_array* inp, int comp, const struct gkyl_basis *basis);

/**
 * Min reduce a gkyl_array evaluating the DG field in each cell at Gauss-Legendre nodes
 * and reducing those values. Allows for reducing a single vector component for arrays
 * with multiple scalar fields within it.
 *
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param inp A gkyl_array to be reduced.
 * @param comp Vector component to reduce.
 * @param basis Baisis DG coefficients expand in (device pointer).
 */
void gkyl_array_reducec_dg_min_cu(double *out_d, const struct gkyl_array* inp, int comp, const struct gkyl_basis *basis);

/**
 * Sum reduce a gkyl_array evaluating the DG field in each cell at Gauss-Legendre nodes
 * and reducing those values. Allows for reducing a single vector component for arrays
 * with multiple scalar fields within it.
 *
 * @param inp A gkyl_array to be reduced.
 * @param out_d A device array with as many elements as the 'inp' has components.
 * @param comp Vector component to reduce.
 * @param basis Baisis DG coefficients expand in (device pointer).
 */
void gkyl_array_reducec_dg_sum_cu(double *out_d, const struct gkyl_array* inp, int comp, const struct gkyl_basis *basis);

