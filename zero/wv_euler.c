#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_euler_prim.h>
#include <gkyl_wv_euler.h>

static const int dir_shuffle[][4] = {
  { 0, 1, 2, 3},
  { 0, 2, 3, 1},
  { 0, 3, 1, 2}
};

struct wv_euler {
    struct gkyl_wv_eqn eqn; // base object
    double gas_gamma; // gas adiabatic constant
};

static void
euler_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_euler *euler = container_of(base, struct wv_euler, eqn);
  gkyl_free(euler);
}

// Waves and speeds using Roe averaging
static double
wave_roe(const struct gkyl_wv_eqn *eqn, 
  int dir, const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  const int *d = dir_shuffle[dir];
  double gas_gamma = euler->gas_gamma;
  double g1 = gas_gamma - 1;

  double rhol = ql[0], rhor = qr[0];
  double pl = gkyl_euler_pressure(gas_gamma, ql), pr = gkyl_euler_pressure(gas_gamma, qr);

  // Roe averages: see Roe's original 1986 paper or LeVeque book
  double srrhol = sqrt(rhol), srrhor = sqrt(rhor);
  double ravgl1 = 1/srrhol, ravgr1 = 1/srrhor;
  double ravg2 = 1/(srrhol+srrhor);
  double u = (ql[d[1]]*ravgl1 + qr[d[1]]*ravgr1)*ravg2;
  double v = (ql[d[2]]*ravgl1 + qr[d[2]]*ravgr1)*ravg2;
  double w = (ql[d[3]]*ravgl1 + qr[d[3]]*ravgr1)*ravg2;
  double enth = ((ql[4]+pl)*ravgl1 + (qr[4]+pr)*ravgr1)*ravg2;

  // See http://ammar-hakim.org/sj/euler-eigensystem.html for notation
  // and meaning of these terms
  double q2 = u*u+v*v+w*w;
  double aa2 = g1*(enth-0.5*q2);
  double a = sqrt(aa2);
  double g1a2 = g1/aa2, euv = enth-q2;

  // Compute projections of jump
  double a4 = g1a2*(euv*delta[0] + u*delta[d[1]] + v*delta[d[2]] + w*delta[d[3]] - delta[4]);
  double a2 = delta[d[2]] - v*delta[0];
  double a3 = delta[d[3]] - w*delta[0];
  double a5 = 0.5*(delta[d[1]] + (a-u)*delta[0] - a*a4)/a;
  double a1 = delta[0] - a4 - a5;

  double *wv;
  // Wave 1: eigenvalue is u-c
  wv = &waves[0];
  wv[0] = a1;
  wv[d[1]] = a1*(u-a);
  wv[d[2]] = a1*v;
  wv[d[3]] = a1*w;
  wv[4] = a1*(enth-u*a);
  s[0] = u-a;

  // Wave 2: eigenvalue is u, u, u three waves are lumped into one
  wv = &waves[5];
  wv[0] = a4;
  wv[d[1]] = a4*u;
  wv[d[2]] = a4*v + a2;
  wv[d[3]] = a4*w + a3;
  wv[4] = a4*0.5*q2 + a2*v + a3*w;
  s[1] = u;

  // Wave 3: eigenvalue is u+c
  wv = &waves[10];
  wv[0] = a5;
  wv[d[1]] = a5*(u+a);
  wv[d[2]] = a5*v;
  wv[d[3]] = a5*w;
  wv[4] = a5*(enth+u*a);
  s[2] = u+a;
  
  return fabs(u)+a;
}

static void
qfluct(const struct gkyl_wv_eqn *eqn, 
  int dir, const double *ql, const double *qr, const double *waves, const double *s,
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

struct gkyl_wv_eqn*
gkyl_wv_euler_new(double gas_gamma)
{
  struct wv_euler *euler = gkyl_malloc(sizeof(struct wv_euler));

  euler->eqn.num_equations = 5;
  euler->eqn.num_waves = 3;
  euler->gas_gamma = gas_gamma;
  euler->eqn.waves_func = wave_roe;
  euler->eqn.qfluct_func = qfluct;

  euler->eqn.ref_count = (struct gkyl_ref_count) { euler_free, 1 };

  return &euler->eqn;
}
