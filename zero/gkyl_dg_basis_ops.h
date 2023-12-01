#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

/**
 * Given values and gradients at the edges of a 1D cell, compute the
 * expansion cofficients of modal cubic (p=3) basis functions. The
 * output @a coeff must be preallocated and have at least 4 elements.
 *
 * @param val val[0] and val[1] are the values at left/right edges
 * @param grad grad[0] and grad[1] are the gradients at left/right edges
 * @param coeff On output, the DG expansion coefficients for p=3 basis.
 */
void gkyl_dg_calc_cubic_1d(const double val[2], const double grad[2], double *coeff);
