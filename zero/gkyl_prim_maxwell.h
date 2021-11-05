#pragma once

/**
 * Compute maximum absolute speed.
 * 
 * @param c Speed of light
 * @param e_fact Correction speed for div(E) correction
 * @param b_fact Correction speed for div(B) correction
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
double gkyl_maxwell_max_abs_speed(double c, double e_fact, double b_fact, const double q[8]);

/**
 * Compute flux.
 * 
 * @param c Speed of light
 * @param e_fact Correction speed for div(E) correction
 * @param b_fact Correction speed for div(B) correction
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_maxwell_flux(double c, double e_fact, double b_fact, const double q[8], double flux[8]);
