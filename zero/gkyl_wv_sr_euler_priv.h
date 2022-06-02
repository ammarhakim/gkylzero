#pragma once

// Private header, not for direct use in user code

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct wv_sr_euler {
  struct gkyl_wv_eqn eqn; // base object
  double gas_gamma; // gas adiabatic constant
};

/**
 * Free SR Euler eqn object.
 *
 * @param ref Reference counter for SR Euler eqn
 */
void gkyl_sr_euler_free(const struct gkyl_ref_count *ref);

/**
 * Compute primitive variables given conserved variables.
 * 
 * See "Grid-based Methods in Relativistic Hydrodynamics and
 * Magnetohydrodynamics", Marti and Muller, Living Reviews in
 * Comp. Astro. vol 1 (3), 2015. URL:
 * https://link.springer.com/article/10.1007/lrca-2015-3
 *
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @param v Primitive variables (output)
 */
GKYL_CU_D
static inline void
gkyl_sr_euler_prim_vars(double gas_gamma, const double q[5], double v[5])
{
  double us=0., vs=0., ws=0., q2s = 0., cs2 = 0.;
  double gammas=0., rhos=0., rhoEpss = 0., fs = 0., dfs = 0., fac0 = 1.;
  double g1 = gas_gamma - 1;
  double ps = 0., ps2 = 1.; 
  double tol = 1.e-6;
  size_t iter = 0;
  
  while (fabs(ps2 - ps) > tol) {
    iter += 1;
    ps = ps2;
    fac0 = q[1] + ps;
    us = q[2] / fac0;
    vs = q[3] / fac0;
    ws = q[4] / fac0;
    q2s = us*us + vs*vs + ws*ws;
    //gammas = 1. / sqrt(1. - q2s);
    gammas = pow(fabs(1. - q2s), -0.5);
    
    rhos = q[0] / gammas;
    rhoEpss = (q[1]  - gammas*q[0] + ps*(1 - gammas*gammas)) / (gammas*gammas);
    fs = g1*rhoEpss - ps; // (gas_gamma - 1)*rhos*eps - ps = p(rhos,epss) - ps, eqn 55
    cs2 = gas_gamma*gammas*gammas*ps / fac0;
    dfs = q2s*cs2 - 1; // eqn 60 for df / dp
    ps2 = ps - fs / dfs;

    //printf("---> Iteration %ld (%g) \n", iter, ps2);
    //printf(" %lg %lg %lg %lg %lg\n", q[0], q[1], q[2], q[3], q[4]);
  }
  //printf("Iterations %ld. Error %lg \n", iter, fabs(ps2-ps));
  
  fac0 = q[1] + ps2;
  us = q[2] / fac0;  
  vs = q[3] / fac0;
  ws = q[4] / fac0;
  q2s = us*us + vs*vs + ws*ws;
  //gammas = 1 / sqrt(1 - q2s);
  gammas = pow(fabs(1. - q2s), -0.5);
  
  v[0] = q[0] / gammas; // rho
  v[1] = ps2; // p
  v[2] = us; 
  v[3] = vs;
  v[4] = ws;
}

/**
 * Compute flux. Assumes rotation to local coordinate system.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
GKYL_CU_D
static void
gkyl_sr_euler_flux(double gas_gamma, const double q[5], double flux[5])
{
  double v[5];
  gkyl_sr_euler_prim_vars(gas_gamma, q, v);
  double pr = v[1];

  double fac0 = q[1] + pr;
  flux[0] = q[0]*q[2]/fac0; // gamma*rho*u
  flux[1] = q[2]; //gamma^2*rho*h*u
  flux[2] = q[2]*q[2]/fac0 + pr; // gamma^2*rho*h*u*u + pr
  flux[3] = q[2]*q[3]/fac0; // gamma^2*rho*h*u*v
  flux[4] = q[2]*q[4]/fac0; // gamma^2*rho*h*u*w
}

GKYL_CU_D
static inline void
rot_to_local(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  // Mass density and energy are scalars
  qlocal[0] = qglobal[0];
  qlocal[1] = qglobal[1];
  // Rotate momentum to local coordinates
  qlocal[2] = qglobal[2]*norm[0] + qglobal[3]*norm[1] + qglobal[4]*norm[2];
  qlocal[3] = qglobal[2]*tau1[0] + qglobal[3]*tau1[1] + qglobal[4]*tau1[2];
  qlocal[4] = qglobal[2]*tau2[0] + qglobal[3]*tau2[1] + qglobal[4]*tau2[2];
}

GKYL_CU_D
static inline void
rot_to_global(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  // Mass density and energy are scalars
  qglobal[0] = qlocal[0];
  qglobal[1] = qlocal[1];
  // Rotate momentum back to global coordinates
  qglobal[2] = qlocal[2]*norm[0] + qlocal[3]*tau1[0] + qlocal[4]*tau2[0];
  qglobal[3] = qlocal[2]*norm[1] + qlocal[3]*tau1[1] + qlocal[4]*tau2[1];
  qglobal[4] = qlocal[2]*norm[2] + qlocal[3]*tau1[2] + qlocal[4]*tau2[2];
  
}

// Waves and speeds using Roe averaging
GKYL_CU_D
static double
wave_roe(const struct gkyl_wv_eqn *eqn, 
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  const struct wv_sr_euler *sr_euler = container_of(eqn, struct wv_sr_euler, eqn);
  double vl[5], vr[5];
  double gas_gamma = sr_euler->gas_gamma;
  double g1 = gas_gamma - 1.;
  double gFrac = gas_gamma/g1;

  // Get prim variables rho, p, u, v, w. 
  gkyl_sr_euler_prim_vars(gas_gamma, ql, vl);
  gkyl_sr_euler_prim_vars(gas_gamma, qr, vr);
  double gammal = 1. / sqrt(1. - (vl[2]*vl[2] + vl[3]*vl[3] + vl[4]*vl[4]));
  double pl = vl[1];
  double gammar = 1. / sqrt(1. - (vr[2]*vr[2] + vr[3]*vr[3] + vr[4]*vr[4]));
  double pr = vr[1];
        
  //Equation numbers in all of the following follows Eulderink AASS 110, 587 (1995)
  //Roe Averages
  //double Kl = sqrt(rhol + gFrac*pl),  Kr = sqrt(rhor + gFrac*pr); //sqrt(10.3)
   
  double Kl = sqrt(ql[1] + pl) / gammal,  Kr = sqrt(qr[1] + pr) / gammar;  //sqrt(10.3)
  double ravgK = 1./(Kl + Kr);
  double v0 = (Kl*gammal + Kr*gammar)*ravgK; //10.7
  double v1 = (Kl*gammal*vl[2] + Kr*gammar*vr[2])*ravgK;
  double v2 = (Kl*gammal*vl[3] + Kr*gammar*vr[3])*ravgK;  
  double v3 = (Kl*gammal*vl[4] + Kr*gammar*vr[4])*ravgK;
  double v4 = (pl/Kl + pr/Kr)*ravgK;
  double cm = 1. - gFrac*v4, cp = 1. + gFrac*v4;
  
  double vava = -v0*v0 + v1*v1 + v2*v2 + v3*v3; //v_alpha v^alpha

  double s2 = 0.5*gas_gamma*v4*(1-vava) - 0.5*g1*(1+vava); //10.13
  double e = v0*v0 - v1*v1; //10.14
  double y = sqrt((1-gas_gamma*v4)*e + s2); //10.14

  
  // Compute projections of jump, Eq 10.16
  double k = v0*delta[1] - v1*delta[2];
  double vada = -v0*delta[1] + v1*delta[2] + v2*delta[3] + v3*delta[4]; // v_alpha Delta^alpha
  double a1 = -(s2*k + sqrt(s2)*y*(v0*delta[2] - v1*delta[1]) + g1*e*(delta[0]+cp*vada))/(2.*e*s2);
  double a2 = -(s2*k - sqrt(s2)*y*(v0*delta[2] - v1*delta[1]) + g1*e*(delta[0]+cp*vada))/(2.*e*s2);
  double a3 = (2.*s2*k + g1*e*(delta[0]+cp*vada))/(e*s2);
  double a4 = delta[3] - k*v2 / e;
  double a5 = delta[4] - k*v3 / e;

  double *wv;
  // Wave 1: eigenvalue is lambda- 
  wv = &waves[0];
  wv[0] = a1*cm;
  wv[1] = a1*(v0-sqrt(s2)*v1/y);
  wv[2] = a1*(v1-sqrt(s2)*v0/y);
  wv[3] = a1*v2;
  wv[4] = a1*v3;
  s[0] = ((1.-gas_gamma*v4)*v0*v1 - sqrt(s2)*y)/((1.-gas_gamma*v4)*v0*v0+s2);//10.12

  // Wave 2: eigenvalue is u, u, u three waves are lumped into one
  wv = &waves[5];
  wv[0] = a3*(cm+s2/g1) - a4*cp*v2 - a5*cp*v3;
  wv[1] = a3*v0;
  wv[2] = a3*v1;
  wv[3] = a3*v2 + a4;
  wv[4] = a3*v3 + a5;
  s[1] = v1 / v0;

  // Wave 3: eigenvalue is lambda+
  wv = &waves[10];
  wv[0] = a2*cm;
  wv[1] = a2*(v0+sqrt(s2)*v1/y);
  wv[2] = a2*(v1+sqrt(s2)*v0/y);
  wv[3] = a2*v2;
  wv[4] = a2*v3;;
  s[2] = ((1.-gas_gamma*v4)*v0*v1 + sqrt(s2)*y)/((1.-gas_gamma*v4)*v0*v0+s2);//10.12
  
  return ((1.-gas_gamma*v4)*v0*fabs(v1) + sqrt(s2)*y)/((1.-gas_gamma*v4)*v0*v0+s2);
}

GKYL_CU_D
static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, 
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[5], *w2 = &waves[10];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i=0; i<5; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i] + s2m*w2[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i] + s2p*w2[i];
  }
}

GKYL_CU_D
static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_sr_euler *sr_euler = container_of(eqn, struct wv_sr_euler, eqn);
  double v[5] = {0.0};
  double gas_gamma = sr_euler->gas_gamma;
  gkyl_sr_euler_prim_vars(gas_gamma, q, v);

  double pr = v[1];
  double gamma = 1. / sqrt(1. - (v[2]*v[2] + v[3]*v[3] + v[4]*v[4]));
  double fac0 = q[1] + pr;
  double v4 = gamma*gamma*pr/fac0;
  double fac1 = 1 - gas_gamma*v4;
  double fac2 = gas_gamma*pr*fac1*(fac0 - q[2]*q[2] / fac0) + gas_gamma*gas_gamma*pr*pr;
  
  return (fac1*q[2] + sqrt(fac2)) / (fac1*fac0 + gas_gamma*pr);
}
