#pragma once

/**
 * Compute maximum absolute speed.
 * 
 * @param vt Thermal velocity
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
double gkyl_iso_euler_max_abs_speed(double vt, const double q[4]);

/**
 * Compute flux.
 * 
 * @param vt Thermal velocity
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_iso_euler_flux(double vt, const double q[4], double flux[4]);
