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
 * Get thermal velocity.
 * 
 * @param wv isothermal Euler equation object
 * @return thermal velocity
 */
double gkyl_wv_iso_euler_vt(const struct gkyl_wv_eqn* wv);