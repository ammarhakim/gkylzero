#pragma once

#include <gkyl_wv_eqn.h>

/**
 * Create a new incompressible Euler equation object for use 
 * by canonical Poisson bracket DG updater.
 * 
 * @return Pointer to incompressible Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_can_pb_incompress_euler_new();

/**
 * Create a new Hasegawa-Mima equation object for use 
 * by canonical Poisson bracket DG updater.
 * 
 * @param kappa Constant density gradient scale length (for turbulence drive).
 * @return Pointer to Hasegawa-Mima equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_can_pb_hasegawa_mima_new(double kappa);

/**
 * Get constant density gradient scale length parameter in Hasegawa-Mima system.
 * 
 * @param wv Hasegawa-Mima equation object
 * @return Constant density gradient scale length (for turbulence drive). 
 */
double gkyl_wv_can_pb_hasegawa_mima_kappa(const struct gkyl_wv_eqn* eqn);

/**
 * Create a new Hasegawa-Wakatani equation object for use 
 * by canonical Poisson bracket DG updater.
 * 
 * @param alpha Adiabaticity parameter for adiabatic coupling of vorticity and density.
 * @param kappa Constant density gradient scale length (for turbulence drive). 
 * @param is_modified Boolean parameter for if we are doing the modified Hasegawa-Wakatani
 *                    system and need to flux surface average the vorticity and density 
 *                    to find the non-zonal contribution to adiabatic coupling RHS
 * @return Pointer to Hasegawa Wakatani equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_can_pb_hasegawa_wakatani_new(double alpha, double kappa, bool is_modified);

/**
 * Get adiabatic coupling parameter
 * 
 * @param wv Hasegawa-Wakatani equation object
 * @return Adiabaticity parameter for adiabatic coupling of vorticity and density.
 */
double gkyl_wv_can_pb_hasegawa_wakatani_alpha(const struct gkyl_wv_eqn* eqn);

/**
 * Get constant density gradient scale length parameter in Hasegawa-Wakatani system.
 * 
 * @param wv Hasegawa-Wakatani equation object
 * @return Constant density gradient scale length (for turbulence drive). 
 */
double gkyl_wv_can_pb_hasegawa_wakatani_kappa(const struct gkyl_wv_eqn* eqn);

/**
 * Determine if Hasegawa-Wakatani system is modified or not.
 * 
 * @param wv Hasegawa-Wakatani equation object
 * @return Boolean is_modified for whether the Hasegawa-Wakatani system is modified.
 */
bool gkyl_wv_can_pb_hasegawa_wakatani_is_modified(const struct gkyl_wv_eqn* eqn);
