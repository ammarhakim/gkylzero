#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_prim_euler.h>
#include <gkyl_wv_euler.h>

static const int dir_shuffle[][3] = {
  {1, 2, 3},
  {2, 3, 1},
  {3, 1, 2}
};

// Make indexing cleaner with the dir_shuffle
#define RHOU d[0]
#define RHOV d[1]
#define RHOW d[2]

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

  // Roe averages: see Roe's original 1981 paper or LeVeque book
  double srrhol = sqrt(rhol), srrhor = sqrt(rhor);
  double ravgl1 = 1/srrhol, ravgr1 = 1/srrhor;
  double ravg2 = 1/(srrhol+srrhor);
  double u = (ql[RHOU]*ravgl1 + qr[RHOU]*ravgr1)*ravg2;
  double v = (ql[RHOV]*ravgl1 + qr[RHOV]*ravgr1)*ravg2;
  double w = (ql[RHOW]*ravgl1 + qr[RHOW]*ravgr1)*ravg2;
  double enth = ((ql[4]+pl)*ravgl1 + (qr[4]+pr)*ravgr1)*ravg2;

  // See http://ammar-hakim.org/sj/euler-eigensystem.html for notation
  // and meaning of these terms
  double q2 = u*u+v*v+w*w;
  double aa2 = g1*(enth-0.5*q2);
  double a = sqrt(aa2);
  double g1a2 = g1/aa2, euv = enth-q2;

  // Compute projections of jump
  double a4 = g1a2*(euv*delta[0] + u*delta[RHOU] + v*delta[RHOV] + w*delta[RHOW] - delta[4]);
  double a2 = delta[RHOV] - v*delta[0];
  double a3 = delta[RHOW] - w*delta[0];
  double a5 = 0.5*(delta[RHOU] + (a-u)*delta[0] - a*a4)/a;
  double a1 = delta[0] - a4 - a5;

  double *wv;
  // Wave 1: eigenvalue is u-c
  wv = &waves[0];
  wv[0] = a1;
  wv[RHOU] = a1*(u-a);
  wv[RHOV] = a1*v;
  wv[RHOW] = a1*w;
  wv[4] = a1*(enth-u*a);
  s[0] = u-a;

  // Wave 2: eigenvalue is u, u, u three waves are lumped into one
  wv = &waves[5];
  wv[0] = a4;
  wv[RHOU] = a4*u;
  wv[RHOV] = a4*v + a2;
  wv[RHOW] = a4*w + a3;
  wv[4] = a4*0.5*q2 + a2*v + a3*w;
  s[1] = u;

  // Wave 3: eigenvalue is u+c
  wv = &waves[10];
  wv[0] = a5;
  wv[RHOU] = a5*(u+a);
  wv[RHOV] = a5*v;
  wv[RHOW] = a5*w;
  wv[4] = a5*(enth+u*a);
  s[2] = u+a;
  
  return fabs(u)+a;
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, 
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

static double
max_speed(const struct gkyl_wv_eqn *eqn, int dir, const double *q)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  return gkyl_euler_max_abs_speed(dir, euler->gas_gamma, q);
}

struct gkyl_wv_eqn*
gkyl_wv_euler_new(double gas_gamma)
{
  struct wv_euler *euler = gkyl_malloc(sizeof(struct wv_euler));

  euler->eqn.type = GKYL_EQN_EULER;
  euler->eqn.num_equations = 5;
  euler->eqn.num_waves = 3;
  euler->gas_gamma = gas_gamma;
  euler->eqn.waves_func = wave_roe;
  euler->eqn.qfluct_func = qfluct_roe;
  euler->eqn.max_speed_func = max_speed;

  euler->eqn.ref_count = (struct gkyl_ref_count) { euler_free, 1 };

  return &euler->eqn;
}

double
gkyl_wv_euler_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  return euler->gas_gamma;
}
