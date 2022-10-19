#pragma once

#include <gkyl_wv_eqn.h>

/**
 * Create a new advection equation object.
 * 
 * @param c advection speed
 * @return Pointer to Burgers equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_advect_new(double c);
