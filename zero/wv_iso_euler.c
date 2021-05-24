#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_prim_iso_euler.h>
#include <gkyl_wv_iso_euler.h>

static const int dir_shuffle[][3] = {
  {1, 2, 3},
  {2, 3, 1},
  {3, 1, 2}
};

// Make indexing cleaner with the dir_shuffle
#define RHOU d[0]
#define RHOV d[1]
#define RHOW d[2]

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
  double u = (ql[RHOU]*ravgl1 + qr[RHOU]*ravgr1)*ravg2;
  double v = (ql[RHOV]*ravgl1 + qr[RHOV]*ravgr1)*ravg2;
  double w = (ql[RHOW]*ravgl1 + qr[RHOW]*ravgr1)*ravg2;

  // Compute projections of jump
  double a0 = delta[0]*(vt+u)/vt/2.0-delta[RHOU]/vt/2.0;
  double a1 = delta[RHOV]-delta[0]*v;
  double a2 = delta[RHOW]-delta[0]*w;
  double a3 = delta[0]*(vt-u)/vt/2.0+delta[RHOU]/vt/2.0;

  double *wv;
  // Wave 1: eigenvalue is u-vt
  wv = &waves[0];
  wv[0] = a0;
  wv[RHOU] = a0*(u-vt);
  wv[RHOV] = a0*v;
  wv[RHOW] = a0*w;
  s[0] = u-vt;

  // Wave 2: eigenvalue is u & u, two waves are lumped into one
  wv = &waves[4];
  wv[0] = 0.0;
  wv[RHOU] = 0.0;
  wv[RHOV] = a1;
  wv[RHOW] = a2;
  s[1] = u;

  // Wave 3: eigenvalue is u+vt
  wv = &waves[8];
  wv[0] = a3;
  wv[RHOU] = a3*(u+vt);
  wv[RHOV] = a3*v;
  wv[RHOW] = a3*w;
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

static double
max_speed(const struct gkyl_wv_eqn *eqn, int dir, const double *q)
{
  const struct wv_iso_euler *iso_euler = container_of(eqn, struct wv_iso_euler, eqn);
  return gkyl_iso_euler_max_abs_speed(dir, iso_euler->vt, q);
}

struct gkyl_wv_eqn*
gkyl_wv_iso_euler_new(double vt)
{
  struct wv_iso_euler *iso_euler = gkyl_malloc(sizeof(struct wv_iso_euler));

  iso_euler->eqn.type = GKYL_ISO_EULER;
  iso_euler->eqn.num_equations = 4;
  iso_euler->eqn.num_waves = 3;
  iso_euler->vt = vt;
  iso_euler->eqn.waves_func = wave_roe;
  iso_euler->eqn.qfluct_func = qfluct_roe;
  iso_euler->eqn.max_speed_func = max_speed;

  iso_euler->eqn.ref_count = (struct gkyl_ref_count) { iso_euler_free, 1 };

  return &iso_euler->eqn;
}
