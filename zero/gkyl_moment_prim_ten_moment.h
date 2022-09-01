#pragma once

/**
 * Computes the primitive variables given the conserved variables.
 * 
 * @param q Conserved variables
 * @param out Primitive variables
 */
static inline void
gkyl_ten_moment_primitive(const double q[10], double out[10])
{
  out[0] = q[0]; 
  out[1] = q[1]/q[0]; 
  out[2] = q[2]/q[0]; 
  out[3] = q[3]/q[0]; 
  out[4] = q[4]-(q[1]*q[1])/q[0]; 
  out[5] = q[5]-(q[1]*q[2])/q[0]; 
  out[6] = q[6]-(q[1]*q[3])/q[0]; 
  out[7] = q[7]-(q[2]*q[2])/q[0]; 
  out[8] = q[8]-(q[2]*q[3])/q[0]; 
  out[9] = q[9]-(q[3]*q[3])/q[0];
}

/**
 * Computes the diagonal components of the pressure tensor
 * 
 * @param q Conserved variables
 * @param out [Pxx, Pyy, Pzz]
 */
static inline void
gkyl_ten_moment_diag_pressure(const double q[10], double out[3])
{
  out[0] = q[4]-(q[1]*q[1])/q[0]; // pxx
  out[1] = q[7]-(q[2]*q[2])/q[0]; // pyy
  out[2] = q[9]-(q[3]*q[3])/q[0]; // pzz
}

/**
 * Compute maximum absolute speed.
 * 
 * @param dir Direction 
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
double gkyl_ten_moment_max_abs_speed(const double q[10]);

/**
 * Compute flux in direction 'dir'.
 * 
 * @param dir Direction to compute flux
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_ten_moment_flux(const double q[10], double flux[10]);
