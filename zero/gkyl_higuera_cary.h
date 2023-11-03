#pragma once

#include <math.h>

/**
 * Computes the update of du/dt = q/m (E + v x B) using the Higuera and Cary
 * algorithm
 * 
 * @param u inital momentum
 * @param q species charge
 * @param m species mass
 * @param dt timestep size
 * @param c speed of light
 * @param E electric field (required at the half-timestep)
 * @param B mangetic field (required at the half-timestep)
 */
static inline void
higuera_cary_push(double u[3], const double q, const double m, const double dt,
const double c, const double E[3], const double B[3])
{
  const double qmdt = q*0.5*dt/m;
  const double E_0 = qmdt*E[0];
  const double E_1 = qmdt*E[1];
  const double E_2 = qmdt*E[2];
  const double B_0 = qmdt*B[0];
  const double B_1 = qmdt*B[1];
  const double B_2 = qmdt*B[2];

  const double u_0_minus = u[0] + E_0;
  const double u_1_minus = u[1] + E_1;
  const double u_2_minus = u[2] + E_2;

  const double u_star = u_0_minus*(B_0/c) + u_1_minus*(B_1/c) + u_2_minus*(B_2/c); 
  const double gamma_minus = sqrt(1 + (u_0_minus*u_0_minus + u_1_minus*u_1_minus + u_2_minus*u_2_minus)/(c*c)); 
  const double dot_tau_tau = (B_0*B_0 + B_1*B_1 + B_2*B_2); 
  const double sigma = gamma_minus*gamma_minus - dot_tau_tau; 
  const double gamma_new = sqrt(  0.5*(   sigma + sqrt( sigma*sigma + 4*(dot_tau_tau + u_star*u_star ) ) )  ); 

  const double t_0 = B_0/gamma_new; 
  const double t_1 = B_1/gamma_new; 
  const double t_2 = B_2/gamma_new; 
  const double s = 1/(1+(t_0*t_0 + t_1*t_1 + t_2*t_2)); 

  const double umt = u_0_minus*t_0 + u_1_minus*t_1 + u_2_minus*t_2;
  const double u_0_plus = s*( u_0_minus + umt*t_0 + (u_1_minus*t_2 - u_2_minus*t_1));
  const double u_1_plus = s*( u_1_minus + umt*t_1 + (u_2_minus*t_0 - u_0_minus*t_2));
  const double u_2_plus = s*( u_2_minus + umt*t_2 + (u_0_minus*t_1 - u_1_minus*t_0));
  u[0] = u_0_plus + E_0 + (u_1_plus*t_2 - u_2_plus*t_1);
  u[1] = u_1_plus + E_1 + (u_2_plus*t_0 - u_0_plus*t_2);
  u[2] = u_2_plus + E_2 + (u_0_plus*t_1 - u_1_plus*t_0);
}