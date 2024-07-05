#pragma once

#include <gkyl_wv_eqn.h>
#include <gkyl_gr_spacetime.h>

// Type of Riemann-solver to use:
enum gkyl_wv_gr_euler_tetrad_rp {
  WV_GR_EULER_TETRAD_RP_LAX = 0, // Default (Lax fluxes).
  WV_GR_EULER_TETRAD_RP_ROE,
};

// Input context, packaged as a struct.
struct gkyl_wv_gr_euler_tetrad_inp {
  double gas_gamma; // Adiabatic index.
  struct gkyl_gr_spacetime *spacetime; // Pointer to base spacetime object.
  enum gkyl_wv_gr_euler_tetrad_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new general relativistic Euler equations object using the tetrad formalism.
*
* @param gas_gamma Adiabatic index.
* @param spacetime Pointer to base spacetime object.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the general relativistic Euler equations object using the tetrad formalism.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_euler_tetrad_new(double gas_gamma, struct gkyl_gr_spacetime* spacetime, bool use_gpu);

/**
* Create a new general relativistic Euler equations object using the tetrad formalism, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the general relativistic Euler equations object using the tetrad formalism.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_euler_tetrad_inew(const struct gkyl_wv_gr_euler_tetrad_inp* inp);

/**
* Get adiabatic index.
*
* @param wv General relativistic Euler equations object using the tetrad formalism.
* @return Adiabatic index.
*/
double
gkyl_wv_gr_euler_tetrad_gas_gamma(const struct gkyl_wv_eqn* wv);

/**
* Get base spacetime object.
*
* @param wv General relativistic Euler equations object using the tetrad formalism.
* @return Pointer to the base spacetime object.
*/
struct gkyl_gr_spacetime*
gkyl_wv_gr_euler_tetrad_spacetime(const struct gkyl_wv_eqn* wv);