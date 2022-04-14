#pragma once

#include <gkyl_wv_eqn.h>

/**
 * Create a new Ten moment equation object.
 * 
 * @param k0 Closure parameter
 * @return Pointer to Ten moment equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_ten_moment_new(double k0);

/**
 * Get closure parameter.
 * 
 * @param wv Ten-moment equation object
 * @return Closure parameter
 */
double gkyl_wv_ten_moment_k0(const struct gkyl_wv_eqn* wv);
