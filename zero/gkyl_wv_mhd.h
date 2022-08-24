#pragma once

#include <gkyl_wv_eqn.h>

// flags to indicate which divergence constraint scheme to use
enum gkyl_wv_mhd_div_constraint {
  GKYL_MHD_DIVB_NONE,
  GKYL_MHD_DIVB_EIGHT_WAVES,
  GKYL_MHD_DIVB_GLM
};

struct wv_mhd_inp {
  double gas_gamma; // gas adiabatic constant
  enum gkyl_wv_mhd_div_constraint divergence_constraint; // divB correction
  double glm_ch; // factor to use in GLM scheme
};

/**
 * Create a new ideal MHD equation object.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param divb Divergence constraint method
 * @return Pointer to mhd equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_mhd_new(const struct wv_mhd_inp *inp);

/**
 * Get gas adiabatic constant.
 *
 * @param wv mhd equation object
 * @return Get gas adiabatic constant
 */
double gkyl_wv_mhd_gas_gamma(const struct gkyl_wv_eqn* wv);
