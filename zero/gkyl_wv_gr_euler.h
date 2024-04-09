#pragma once

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_gr_euler_rp {
  WV_GR_EULER_RP_LAX = 0, // Default (Lax fluxes).
};

// Input context, packaged as a struct.
struct gkyl_wv_gr_euler_inp {
  double gas_gamma; // Adiabatic index.
  enum gkyl_wv_gr_euler_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new general relativistic Euler equations object.
*
* @param gas_gamma Adiabatic index.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the general relativistic Euler equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_euler_new(double gas_gamma, bool use_gpu);

/**
* Create a new general relativistic Euler equations object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the general relativistic Euler equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_euler_inew(const struct gkyl_wv_gr_euler_inp* inp);

/**
* Get adiabatic index.
*
* @param wv General relativistic Euler equations object.
* @return Adiabatic index.
*/
double
gkyl_wv_gr_euler_gas_gamma(const struct gkyl_wv_eqn* wv);