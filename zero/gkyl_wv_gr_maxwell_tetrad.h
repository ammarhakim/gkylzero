#pragma once

#include <gkyl_wv_eqn.h>
#include <gkyl_gr_spacetime.h>

// Type of Riemann-solver to use:
enum gkyl_wv_gr_maxwell_tetrad_rp {
  WV_GR_MAXWELL_TETRAD_RP_LAX = 0, // Default (Lax fluxes).
  WV_GR_MAXWELL_TETRAD_RP_ROE,
};

// Input context, packaged as a struct.
struct gkyl_wv_gr_maxwell_tetrad_inp {
  double light_speed; // Speed of light.
  double e_fact; // Factor of speed of light for electric field correction.
  double b_fact; // Factor of speed of light for magnetic field correction.

  struct gkyl_gr_spacetime *spacetime; // Pointer to base spacetime object.
  enum gkyl_wv_gr_maxwell_tetrad_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new general relativistic Maxwell equations object in the tetrad basis.
*
* @param light_speed Speed of light.
* @param e_fact Factor of speed of light for electric field correction.
* @param b_fact Factor of speed of light for magnetic field correction.
* @param spacetime Pointer to base spacetime object.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the general relativistic Maxwell equations object in the tetrad basis.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_maxwell_tetrad_new(double light_speed, double e_fact, double b_fact, struct gkyl_gr_spacetime* spacetime, bool use_gpu);

/**
* Create a new general relativistic Maxwell equations object in the tetrad basis, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the general relativistic Maxwell equations object in the tetrad basis.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_maxwell_tetrad_inew(const struct gkyl_wv_gr_maxwell_tetrad_inp* inp);

/**
* Get speed of light.
*
* @param eqn General relativistic Maxwell equations object in the tetrad basis.
* @return Speed of light.
*/
double
gkyl_wv_gr_maxwell_tetrad_light_speed(const struct gkyl_wv_eqn* eqn);

/**
* Get factor of speed of light for electric field correction.
*
* @param eqn General relativistic Maxwell equations object in the tetrad basis.
* @return Factor of speed of light for electric field correction.
*/
double
gkyl_wv_gr_maxwell_tetrad_e_fact(const struct gkyl_wv_eqn* eqn);

/**
* Get factor of speed of light for magnetic field correction.
*
* @param eqn General relativistic Maxwell equations object in the tetrad basis.
* @return Factor of speed of light for magnetic field correction.
*/
double
gkyl_wv_gr_maxwell_tetrad_b_fact(const struct gkyl_wv_eqn* eqn);

/**
* Get base spacetime object.
*
* @param eqn General relativistic Maxwell equations object in the tetrad basis.
* @return Pointer to the base spacetime object.
*/
struct gkyl_gr_spacetime*
gkyl_wv_gr_maxwell_tetrad_spacetime(const struct gkyl_wv_eqn* eqn);