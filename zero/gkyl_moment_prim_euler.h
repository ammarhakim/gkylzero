#pragma once

/**
 * Computes the scalar pressure given the conserved variables.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 */
static inline double
gkyl_euler_pressure(double gas_gamma, const double q[5])
{
  return (gas_gamma-1)*(q[4]-0.5*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3])/q[0]);
}

/**
 * Compute primitive variables given conserved variables.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @param v Primitive variables (output)
 */
static inline void
gkyl_euler_prim_vars(double gas_gamma, const double q[5], double v[5])
{
  v[0] = q[0];
  v[1] = q[1]/q[0];
  v[2] = q[2]/q[0];
  v[3] = q[3]/q[0];
  v[4] = gkyl_euler_pressure(gas_gamma, q);
}

/**
 * Compute maximum absolute speed.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
double gkyl_euler_max_abs_speed(double gas_gamma, const double q[5]);

/**
 * Compute flux given conserved variables
 * 
 * @param gas_gamma Gas adiabatic constant
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_euler_flux(double gas_gamma, const double q[5], double flux[5]);
