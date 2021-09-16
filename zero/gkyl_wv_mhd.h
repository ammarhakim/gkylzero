#pragma once

#include <gkyl_wv_eqn.h>

/**
 * Create a new ideal MHD equation object.
 *
 * @param gas_gamma Gas adiabatic constant
 * @return Pointer to mhd equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_mhd_new(double gas_gamma);

/**
 * Get gas adiabatic constant.
 *
 * @param wv mhd equation object
 * @return Get gas adiabatic constant
 */
double gkyl_wv_mhd_gas_gamma(const struct gkyl_wv_eqn* wv);
