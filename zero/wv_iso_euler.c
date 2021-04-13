#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_wv_iso_euler.h>

static const int dir_shuffle[][4] = {
  { 0, 1, 2, 3},
  { 0, 2, 3, 1},
  { 0, 3, 1, 2}
};

struct wv_iso_euler {
    struct gkyl_wv_eqn eqn; // base object
    double vt; // thermal velocity
};

static void
iso_euler_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_iso_euler *iso_euler = container_of(base, struct wv_iso_euler, eqn);
  gkyl_free(iso_euler);
}

// Waves and speeds using Roe averaging
static double
wave_roe(const struct gkyl_wv_eqn *eqn, 
  int dir, const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  const struct wv_iso_euler *iso_euler = container_of(eqn, struct wv_iso_euler, eqn);
  const int *d = dir_shuffle[dir];
  double vt = iso_euler->vt;

  double rhol = ql[0], rhor = qr[0];

  // Roe averages: see Roe's original 1981 paper or LeVeque book
  double srrhol = sqrt(rhol), srrhor = sqrt(rhor);
  double ravgl1 = 1/srrhol, ravgr1 = 1/srrhor;
  double ravg2 = 1/(srrhol+srrhor);
  double u = (ql[d[1]]*ravgl1 + qr[d[1]]*ravgr1)*ravg2;
  double v = (ql[d[2]]*ravgl1 + qr[d[2]]*ravgr1)*ravg2;
  double w = (ql[d[3]]*ravgl1 + qr[d[3]]*ravgr1)*ravg2;

  // Compute projections of jump
  double a0 = delta[0]*(vt+u)/vt/2.0-delta[d[1]]/vt/2.0;
  double a1 = delta[d[2]]-delta[0]*v;
  double a2 = delta[d[3]]-delta[0]*w;
  double a3 = delta[0]*(vt-u)/vt/2.0+delta[d[1]]/vt/2.0;

  double *wv;
  // Wave 1: eigenvalue is u-vt
  wv = &waves[0];
  wv[0] = a0;
  wv[d[1]] = a0*(u-vt);
  wv[d[2]] = a0*v;
  wv[d[3]] = a0*w;
  s[0] = u-vt;

  // Wave 2: eigenvalue is u & u, two waves are lumped into one
  wv = &waves[4];
  wv[0] = 0.0;
  wv[d[1]] = 0.0;
  wv[d[2]] = a1;
  wv[d[3]] = a2;
  s[1] = u;

  // Wave 3: eigenvalue is u+vt
  wv = &waves[8];
  wv[0] = a3;
  wv[d[1]] = a3*(u+vt);
  wv[d[2]] = a3*v;
  wv[d[3]] = a3*w;
  s[2] = u+vt;
  
  return fabs(u)+vt;
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn, 
  int dir, const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[4], *w2 = &waves[8];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i=0; i<4; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i] + s2m*w2[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i] + s2p*w2[i];
  }
}

struct gkyl_wv_eqn*
gkyl_wv_iso_euler_new(double vt)
{
  struct wv_iso_euler *iso_euler = gkyl_malloc(sizeof(struct wv_iso_euler));

  iso_euler->eqn.num_equations = 4;
  iso_euler->eqn.num_waves = 3;
  iso_euler->vt = vt;
  iso_euler->eqn.waves_func = wave_roe;
  iso_euler->eqn.qfluct_func = qfluct_roe;

  iso_euler->eqn.ref_count = (struct gkyl_ref_count) { iso_euler_free, 1 };

  return &iso_euler->eqn;
}
