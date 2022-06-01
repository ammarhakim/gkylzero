#pragma once

#include <gkyl_wv_eqn.h>

/**
 * Create a new isothermal Euler equation object.
 * 
 * @param vt Thermal velocity
 * @return Pointer to isothermal Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_iso_euler_new(double vt);

/**
 * Create a new isothermal Euler equation object that lives on NV-GPU.
 * 
 * @param vt Thermal velocity
 * @return Pointer to isothermal Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_iso_euler_cu_dev_new(double vt);