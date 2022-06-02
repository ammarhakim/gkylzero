#pragma once

// Private header, not for direct use in user code

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct wv_euler {
  struct gkyl_wv_eqn eqn; // base object
  double gas_gamma; // gas adiabatic constant
};

/**
 * Free Euler eqn object.
 *
 * @param ref Reference counter for Euler eqn
 */
void gkyl_euler_free(const struct gkyl_ref_count *ref);

/**
 * Computes the scalar pressure given the conserved variables.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 */
GKYL_CU_D
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
GKYL_CU_D
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
 * Compute flux. Assumes rotation to local coordinate system.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
GKYL_CU_D
static void
gkyl_euler_flux(double gas_gamma, const double q[5], double flux[5])
{
  double pr = gkyl_euler_pressure(gas_gamma, q), u = q[1]/q[0];
  flux[0] = q[1]; // rho*u
  flux[1] = q[1]*u + pr; // rho*u*u + pr
  flux[2] = q[2]*u; // rho*v*u
  flux[3] = q[3]*u; // rho*w*u
  flux[4] = (q[4]+pr)*u; // (E+p)*u
}

// Euler perfectly reflecting wall
GKYL_CU_D
static void
euler_wall(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  // copy density and pressure
  ghost[0] = skin[0];
  ghost[4] = skin[4];

  // zero-normal for momentum
  ghost[1] = -skin[1];
  ghost[2] = skin[2];
  ghost[3] = skin[3];
}

GKYL_CU_D
static inline void
rot_to_local(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = qglobal[1]*norm[0] + qglobal[2]*norm[1] + qglobal[3]*norm[2];
  qlocal[2] = qglobal[1]*tau1[0] + qglobal[2]*tau1[1] + qglobal[3]*tau1[2];
  qlocal[3] = qglobal[1]*tau2[0] + qglobal[2]*tau2[1] + qglobal[3]*tau2[2];
  qlocal[4] = qglobal[4];
}

GKYL_CU_D
static inline void
rot_to_global(const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = qlocal[1]*norm[0] + qlocal[2]*tau1[0] + qlocal[3]*tau2[0];
  qglobal[2] = qlocal[1]*norm[1] + qlocal[2]*tau1[1] + qlocal[3]*tau2[1];
  qglobal[3] = qlocal[1]*norm[2] + qlocal[2]*tau1[2] + qlocal[3]*tau2[2];
  qglobal[4] = qlocal[4];
}

// Waves and speeds using Roe averaging
GKYL_CU_D
static double
wave_roe(const struct gkyl_wv_eqn *eqn, 
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double gas_gamma = euler->gas_gamma;
  double g1 = gas_gamma - 1;

  double rhol = ql[0], rhor = qr[0];
  double pl = gkyl_euler_pressure(gas_gamma, ql), pr = gkyl_euler_pressure(gas_gamma, qr);

  // Roe averages: see Roe's original 1981 paper or LeVeque book
  double srrhol = sqrt(rhol), srrhor = sqrt(rhor);
  double ravgl1 = 1/srrhol, ravgr1 = 1/srrhor;
  double ravg2 = 1/(srrhol+srrhor);
  double u = (ql[1]*ravgl1 + qr[1]*ravgr1)*ravg2;
  double v = (ql[2]*ravgl1 + qr[2]*ravgr1)*ravg2;
  double w = (ql[3]*ravgl1 + qr[3]*ravgr1)*ravg2;
  double enth = ((ql[4]+pl)*ravgl1 + (qr[4]+pr)*ravgr1)*ravg2;

  // See http://ammar-hakim.org/sj/euler-eigensystem.html for notation
  // and meaning of these terms
  double q2 = u*u+v*v+w*w;
  double aa2 = g1*(enth-0.5*q2);
  double a = sqrt(aa2);
  double g1a2 = g1/aa2, euv = enth-q2;

  // Compute projections of jump
  double a4 = g1a2*(euv*delta[0] + u*delta[1] + v*delta[2] + w*delta[3] - delta[4]);
  double a2 = delta[2] - v*delta[0];
  double a3 = delta[3] - w*delta[0];
  double a5 = 0.5*(delta[1] + (a-u)*delta[0] - a*a4)/a;
  double a1 = delta[0] - a4 - a5;

  double *wv;
  // Wave 1: eigenvalue is u-c
  wv = &waves[0];
  wv[0] = a1;
  wv[1] = a1*(u-a);
  wv[2] = a1*v;
  wv[3] = a1*w;
  wv[4] = a1*(enth-u*a);
  s[0] = u-a;

  // Wave 2: eigenvalue is u, u, u three waves are lumped into one
  wv = &waves[5];
  wv[0] = a4;
  wv[1] = a4*u;
  wv[2] = a4*v + a2;
  wv[3] = a4*w + a3;
  wv[4] = a4*0.5*q2 + a2*v + a3*w;
  s[1] = u;

  // Wave 3: eigenvalue is u+c
  wv = &waves[10];
  wv[0] = a5;
  wv[1] = a5*(u+a);
  wv[2] = a5*v;
  wv[3] = a5*w;
  wv[4] = a5*(enth+u*a);
  s[2] = u+a;
  
  return fabs(u)+a;
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
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double gas_gamma = euler->gas_gamma;
  double pr = gkyl_euler_pressure(gas_gamma, q), u = q[1]/q[0];
  return fabs(u) + sqrt(gas_gamma*pr/q[0]);
}
