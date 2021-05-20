#pragma once

/**
 * Compute maximum absolute speed.
 * 
 * @param dir Direction 
 * @param c Speed of light
 * @param e_fact Correction speed for div(E) correction
 * @param b_fact Correction speed for div(B) correction
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
double gkyl_maxwell_max_abs_speed(int dir, double c, double e_fact, double b_fact, const double q[8]);

/**
 * Compute flux in direction 'dir'.
 * 
 * @param dir Direction to compute flux
 * @param c Speed of light
 * @param e_fact Correction speed for div(E) correction
 * @param b_fact Correction speed for div(B) correction
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_maxwell_flux(int dir, double c, double e_fact, double b_fact, const double q[8], double flux[8]);
