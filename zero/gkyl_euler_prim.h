#pragma once

/**
 * Computes the scalar pressure given the conserved variables.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 */
static inline double gkyl_euler_pressure(double gas_gamma, const double *q)
{
  return (gas_gamma-1)*(q[4]-0.5*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3])/q[0]);
}

/**
 * Compute maximum absolute speed.
 * 
 * @param dir Direction 
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
double gkyl_euler_max_abs_speed(int dir, double gas_gamma, const double *q);

/**
 * Compute flux in direction 'dir'.
 * 
 * @param dir Direction to compute flux
 * @param gas_gamma Gas adiabatic constant
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_euler_flux(int dir, double gas_gamma, const double *q, double *flux);
