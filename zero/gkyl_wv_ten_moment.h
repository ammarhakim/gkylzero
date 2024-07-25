#pragma once

#include <gkyl_wv_eqn.h>

/**
 * Create a new Ten moment equation object.
 * 
 * @param k0 Closure parameter
 * @param use_gpu Boolean to determine whether wave equation object is on host or device
 * @return Pointer to Ten moment equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_ten_moment_new(double k0, bool use_gpu);

/**
 * Create a new Ten moment equation object that lives on NV-GPU.
 * see new() method above for documentation.
 */
struct gkyl_wv_eqn* gkyl_wv_ten_moment_cu_dev_new(double k0);

/**
 * Get closure parameter.
 * 
 * @param wv Ten-moment equation object
 * @return Closure parameter
 */
double gkyl_wv_ten_moment_k0(const struct gkyl_wv_eqn* wv);
