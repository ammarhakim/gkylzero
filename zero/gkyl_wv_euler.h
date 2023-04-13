#pragma once

#include <gkyl_wv_eqn.h>

// Type of Rieman problem solver to use
enum gkyl_wv_euler_rp {
  WV_EULER_RP_ROE = 0, // default
  WV_EULER_RP_HLLC,
  WV_EULER_RP_LAX,
  WV_EULER_RP_HLL,
  WV_EULER_RP_HLLI,
};

// input packaged as a struct
struct gkyl_wv_euler_inp {
  double gas_gamma; // gas adiabatic constant
  enum gkyl_wv_euler_rp rp_type; // type of RP to use
};

/**
 * Create a new Euler equation object.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @return Pointer to Euler equation object.
 */
struct gkyl_wv_eqn *gkyl_wv_euler_new(double gas_gamma);

/**
 * Create a new Euler equation object.
 * 
 * @param inp Input parameters
 * @return Pointer to Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_euler_inew(const struct gkyl_wv_euler_inp *inp);

/**
 * Get gas adiabatic constant.
 * 
 * @param wv Euler equation object
 * @return Gas adiabatic constant
 */
double gkyl_wv_euler_gas_gamma(const struct gkyl_wv_eqn* wv);
