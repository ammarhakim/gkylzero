#pragma once

#include <gkyl_wv_eqn.h>
#include <gkyl_gr_spacetime.h>

// Type of Riemann-solver to use:
enum gkyl_wv_gr_ultra_rel_euler_tetrad_rp {
  WV_GR_ULTRA_REL_EULER_TETRAD_RP_LAX = 0, // Default (Lax fluxes).
  WV_GR_ULTRA_REL_EULER_TETRAD_RP_ROE,
};

// Input context, packaged as a struct.
struct gkyl_wv_gr_ultra_rel_euler_tetrad_inp {
  double gas_gamma; // Adiabatic index.
  enum gkyl_spacetime_gauge spacetime_gauge; // Spacetime gauge choice.
  struct gkyl_gr_spacetime *spacetime; // Pointer to base spacetime object.

  enum gkyl_wv_gr_ultra_rel_euler_tetrad_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new general relativistic Euler equations object in the tetrad basis with ultra-relativistic equation of state.
*
* @param gas_gamma Adiabatic index.
* @param spacetime_gauge Spacetime gauge choice.
* @param spacetime Pointer to base spacetime object.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the general relativistic Euler equations object in the tetrad basis with ultra-relativistic equation of state.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_ultra_rel_euler_tetrad_new(double gas_gamma, enum gkyl_spacetime_gauge spacetime_gauge, struct gkyl_gr_spacetime* spacetime, bool use_gpu);

/**
* Create a new general relativistic Euler equations object in the tetrad basis with ultra-relativistic equation of state, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the general relativistic Euler equations object in the tetrad basis with ultra-relativistic equation of state.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_ultra_rel_euler_tetrad_inew(const struct gkyl_wv_gr_ultra_rel_euler_tetrad_inp* inp);

/**
* Get adiabatic index.
*
* @param eqn General relativistic Euler equations object in the tetrad basis with ultra-relativistic equation of state.
* @return Adiabatic index.
*/
double
gkyl_wv_gr_ultra_rel_euler_tetrad_gas_gamma(const struct gkyl_wv_eqn* eqn);

/**
* Get spacetime gauge choice.
*
* @param eqn General relativistic Euler equations object in the tetrad basis with ultra-relativistic equation of state.
* @return Spacetime gauge choice.
*/
enum gkyl_spacetime_gauge
gkyl_wv_gr_ultra_rel_euler_tetrad_spacetime_gauge(const struct gkyl_wv_eqn* eqn);

/**
* Get base spacetime object.
*
* @param eqn General relativistic Euler equations object in the tetrad basis with ultra-relativistic equation of state.
* @return Pointer to the base spacetime object.
*/
struct gkyl_gr_spacetime*
gkyl_wv_gr_ultra_rel_euler_tetrad_spacetime(const struct gkyl_wv_eqn* eqn);