#pragma once

#include <gkyl_wv_eqn.h>

/**
 * Create a new Ten moment equation object.
 *
 * @param k0 Closure parameter
 * @param use_grad_closure Should we use gradient-based closure?
 * @return Pointer to Ten moment equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_ten_moment_new(double k0, bool use_grad_closure);

/**
 * Get closure parameter.
 * 
 * @param wv Ten-moment equation object
 * @return Closure parameter
 */
double gkyl_wv_ten_moment_k0(const struct gkyl_wv_eqn* wv);

/**
 * Should we use grad-closure?
 * 
 * @param wv Ten-moment equation object
 * @return True if grad closure, false otherwise
 */
double gkyl_wv_ten_moment_use_grad_closure(const struct gkyl_wv_eqn* wv);
