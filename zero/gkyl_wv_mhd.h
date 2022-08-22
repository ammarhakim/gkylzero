#pragma once

#include <gkyl_wv_eqn.h>

// flags to indicate which divergence constraint scheme to use
enum gkyl_wv_mhd_div_constraint {
  GKYL_MHD_DIVB_NONE,
  GKYL_MHD_DIVB_EIGHT_WAVES,
  GKYL_MHD_DIVB_GLM
};

/**
 * Create a new ideal MHD equation object.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param divb Divergence constraint method
 * @return Pointer to mhd equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_mhd_new(double gas_gamma, enum gkyl_wv_mhd_div_constraint divb);

/**
 * Get gas adiabatic constant.
 *
 * @param wv mhd equation object
 * @return Get gas adiabatic constant
 */
double gkyl_wv_mhd_gas_gamma(const struct gkyl_wv_eqn* wv);
