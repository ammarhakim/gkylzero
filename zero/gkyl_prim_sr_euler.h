#pragma once

#include <math.h>

/**
 * Compute primitive variables given conserved variables. 
 * Method from Marti and Muller Living Reviews Comp Astro 2015, section 4.6
 *
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @param v Primitive variables (output)
 */
static inline void gkyl_sr_euler_prim_vars(double gas_gamma, const double q[5], double v[5])
{
  double us=0., vs=0., ws=0., q2s = 0., cs2 = 0.;
  double gammas=0., rhos=0., rhoEpss = 0., fs = 0., dfs = 0., fac0 = 1.;
  double g1 = gas_gamma - 1;
  double ps = 0., ps2 = 1.; 
  double tol = 1.e-6;
  // int iter = 0;
  
  while (fabs(ps2 - ps) > tol) {
    // iter += 1;
    // printf("Iteration %i \n", iter);
    ps = ps2;
    fac0 = q[1] + ps;
    us = q[2] / fac0;  
    vs = q[3] / fac0;
    ws = q[4] / fac0;
    q2s = us*us + vs*vs + ws*ws;
    gammas = 1. / sqrt(1. - q2s);
    
    rhos = q[0] / gammas;
    
    rhoEpss = (q[1]  - gammas*q[0] + ps*(1 - gammas*gammas)) / (gammas*gammas);

    fs = g1*rhoEpss - ps; //(gas_gamma - 1)*rhos*eps - ps = p(rhos,epss) - ps, eqn 55

    cs2 = gas_gamma*gammas*gammas*ps / fac0;
    dfs = q2s*cs2 - 1; //eqn 60 for df / dp

    ps2 = ps - fs / dfs;
   

  }
  
  fac0 = q[1] + ps2;
  us = q[2] / fac0;  
  vs = q[3] / fac0;
  ws = q[4] / fac0;
  q2s = us*us + vs*vs + ws*ws;
  gammas = 1 / sqrt(1 - q2s);
						     
  
  v[0] = q[0] / gammas; // rho
  v[1] = ps2; //p
  v[2] = us; 
  v[3] = vs;
  v[4] = ws;
}

/**
 * Compute maximum absolute speed.
 *  
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
double gkyl_sr_euler_max_abs_speed(double gas_gamma, const double q[5]);

/**
 * Compute flux. Assumes rotation to local coordinate system.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_sr_euler_flux(double gas_gamma, const double q[5], double flux[5]);
