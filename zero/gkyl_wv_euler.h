#pragma once

#include <gkyl_wv_eqn.h>

/**
 * Create a new Euler equation object.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @return Pointer to Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_euler_new(double gas_gamma);

/**
 * Create a new Euler equation object that lives on NV-GPU.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @return Pointer to Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_euler_cu_dev_new(double gas_gamma);

/**
 * Get gas adiabatic constant.
 * 
 * @param wv Euler equation object
 * @return Gas adiabatic constant
 */
double gkyl_wv_euler_gas_gamma(const struct gkyl_wv_eqn* wv);
