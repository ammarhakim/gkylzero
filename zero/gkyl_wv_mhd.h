#pragma once

#include <gkyl_wv_eqn.h>

/**
 * Create a new ideal MHD equation object.
 *
 * @param gas_gamma Gas adiabatic constant
 * @return Pointer to mhd equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_mhd_new(
    double gas_gamma, const char *divergence_constraint);

/**
 * Get gas adiabatic constant.
 *
 * @param wv mhd equation object
 * @return Get gas adiabatic constant
 */
double gkyl_wv_mhd_gas_gamma(const struct gkyl_wv_eqn* wv);

/**
 * Free wave mhd eqn object.
 *
 * @param ref Reference counter for mhd eqn
 */
void
gkyl_wv_mhd_free(const struct gkyl_ref_count *ref);

/**
 * Create a new ideal MHD equation object on gpu.
 *
 * @param gas_gamma Gas adiabatic constant
 * @return Pointer to mhd equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_mhd_cu_dev_new(
    double gas_gamma, const char *divergence_constraint);

