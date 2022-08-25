#pragma once

#include <gkyl_wv_eqn.h>

// flags to indicate which divergence constraint scheme to use
enum gkyl_wv_mhd_div_constraint {
  GKYL_MHD_DIVB_NONE,
  GKYL_MHD_DIVB_EIGHT_WAVES,
  GKYL_MHD_DIVB_GLM
};

struct wv_mhd {
  struct gkyl_wv_eqn eqn; // base object
  double gas_gamma; // gas adiabatic constant
  enum gkyl_wv_mhd_div_constraint divergence_constraint; // divB correction
  double glm_ch; // factor to use in GLM scheme
  double glm_cp; // factor to use in GLM scheme
  double glm_alpha; // Mignone & Tzeferacos, JCP (2010) 229, 2117, Eq (27).
};

struct gkyl_wv_mhd_inp {
  double gas_gamma; // gas adiabatic constant
  enum gkyl_wv_mhd_div_constraint divergence_constraint; // divB correction
  double glm_ch; // factor to use in GLM scheme
  double glm_cp; // factor to use in GLM scheme
  double glm_alpha; // Mignone & Tzeferacos, JCP (2010) 229, 2117, Eq (27).
};

/**
 * Create a new ideal MHD equation object.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param divb Divergence constraint method
 * @return Pointer to mhd equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_mhd_new(const struct gkyl_wv_mhd_inp *inp);

/**
 * Get gas adiabatic constant.
 *
 * @param wv mhd equation object
 * @return Get gas adiabatic constant
 */
double gkyl_wv_mhd_gas_gamma(const struct gkyl_wv_eqn* wv);
