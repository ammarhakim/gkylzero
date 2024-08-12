#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_moment_prim_iso_euler.h>
#include <gkyl_wv_iso_euler.h>

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

static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<4; ++i)
    wout[i] = qin[i];
}
static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<4; ++i)
    qout[i] = win[i];
}

// Isothermal Euler perfectly reflecting wall
static void
iso_euler_wall(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  // copy density
  ghost[0] = skin[0];

  // zero-normal for momentum
  ghost[1] = -skin[1];
  ghost[2] = skin[2];
  ghost[3] = skin[3];
}

// Isothermal Euler no-slip wall
static void
iso_euler_no_slip(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  // copy density and pressure
  ghost[0] = skin[0];

  // zero-normal for momentum
  ghost[1] = -skin[1];
  ghost[2] = -skin[2];
  ghost[3] = -skin[3];
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qglobal, double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = qglobal[1]*norm[0] + qglobal[2]*norm[1] + qglobal[3]*norm[2];
  qlocal[2] = qglobal[1]*tau1[0] + qglobal[2]*tau1[1] + qglobal[3]*tau1[2];
  qlocal[3] = qglobal[1]*tau2[0] + qglobal[2]*tau2[1] + qglobal[3]*tau2[2];
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* GKYL_RESTRICT qlocal, double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = qlocal[1]*norm[0] + qlocal[2]*tau1[0] + qlocal[3]*tau2[0];
  qglobal[2] = qlocal[1]*norm[1] + qlocal[2]*tau1[1] + qlocal[3]*tau2[1];
  qglobal[3] = qlocal[1]*norm[2] + qlocal[2]*tau1[2] + qlocal[3]*tau2[2];
}

// Waves and speeds using Roe averaging
static double
wave_roe(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *ql, const double *qr, 
  double *waves, double *s)
{
  const struct wv_iso_euler *iso_euler = container_of(eqn, struct wv_iso_euler, eqn);
  double vt = iso_euler->vt;

  double rhol = ql[0], rhor = qr[0];

  // Roe averages: see Roe's original 1981 paper or LeVeque book
  double srrhol = sqrt(rhol), srrhor = sqrt(rhor);
  double ravgl1 = 1/srrhol, ravgr1 = 1/srrhor;
  double ravg2 = 1/(srrhol+srrhor);
  double u = (ql[1]*ravgl1 + qr[1]*ravgr1)*ravg2;
  double v = (ql[2]*ravgl1 + qr[2]*ravgr1)*ravg2;
  double w = (ql[3]*ravgl1 + qr[3]*ravgr1)*ravg2;

  // Compute projections of jump
  double a0 = delta[0]*(vt+u)/vt/2.0-delta[1]/vt/2.0;
  double a1 = delta[2]-delta[0]*v;
  double a2 = delta[3]-delta[0]*w;
  double a3 = delta[0]*(vt-u)/vt/2.0+delta[1]/vt/2.0;

  double *wv;
  // Wave 1: eigenvalue is u-vt
  wv = &waves[0];
  wv[0] = a0;
  wv[1] = a0*(u-vt);
  wv[2] = a0*v;
  wv[3] = a0*w;
  s[0] = u-vt;

  // Wave 2: eigenvalue is u & u, two waves are lumped into one
  wv = &waves[4];
  wv[0] = 0.0;
  wv[1] = 0.0;
  wv[2] = a1;
  wv[3] = a2;
  s[1] = u;

  // Wave 3: eigenvalue is u+vt
  wv = &waves[8];
  wv[0] = a3;
  wv[1] = a3*(u+vt);
  wv[2] = a3*v;
  wv[3] = a3*w;
  s[2] = u+vt;
  
  return fabs(u)+vt;
}

static void
qfluct_roe(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, const double *waves, const double *s,
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

// Waves and speeds using Lax fluxes
static double
wave_lax(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *ql, const double *qr, 
  double *waves, double *s)
{
  const struct wv_iso_euler *iso_euler = container_of(eqn, struct wv_iso_euler, eqn);
  double vt = iso_euler->vt;
  
  double sl = gkyl_iso_euler_max_abs_speed(vt, ql);
  double sr = gkyl_iso_euler_max_abs_speed(vt, qr);
  double amax = fmax(sl, sr);

  double fl[4], fr[4];
  gkyl_iso_euler_flux(vt, ql, fl);
  gkyl_iso_euler_flux(vt, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[4];
  for (int i=0; i<4; ++i) {
    w0[i] = 0.5*((qr[i]-ql[i]) - (fr[i]-fl[i])/amax);
    w1[i] = 0.5*((qr[i]-ql[i]) + (fr[i]-fl[i])/amax);
  }

  s[0] = -amax;
  s[1] = amax;
  
  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[4];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i=0; i<4; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i];
  }
}

static double
wave(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, 
  double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return wave_roe(eqn, delta, ql, qr, waves, s);
  else
    return wave_lax(eqn, delta, ql, qr, waves, s);

  return 0.0; // can't happen
}

static void
qfluct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return qfluct_roe(eqn, ql, qr, waves, s, amdq, apdq);
  else
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return q[0] > 0.0; // density should be positive
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_iso_euler *iso_euler = container_of(eqn, struct wv_iso_euler, eqn);
  return gkyl_iso_euler_max_abs_speed(iso_euler->vt, q);
}

static inline void
iso_euler_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 4; i++) {
    sout[i] = 0.0;
  }
}

struct gkyl_wv_eqn*
gkyl_wv_iso_euler_new(double vt)
{
  struct wv_iso_euler *iso_euler = gkyl_malloc(sizeof(struct wv_iso_euler));

  iso_euler->eqn.type = GKYL_EQN_ISO_EULER;
  iso_euler->eqn.num_equations = 4;
  iso_euler->eqn.num_waves = 3;
  iso_euler->eqn.num_diag = 4;
  
  iso_euler->vt = vt;
  iso_euler->eqn.waves_func = wave;
  iso_euler->eqn.qfluct_func = qfluct;

  iso_euler->eqn.check_inv_func = check_inv;
  iso_euler->eqn.max_speed_func = max_speed;
  
  iso_euler->eqn.rotate_to_local_func = rot_to_local;
  iso_euler->eqn.rotate_to_global_func = rot_to_global;

  iso_euler->eqn.cons_to_riem = cons_to_riem;
  iso_euler->eqn.riem_to_cons = riem_to_cons;

  iso_euler->eqn.wall_bc_func = iso_euler_wall;
  iso_euler->eqn.no_slip_bc_func = iso_euler_no_slip;

  iso_euler->eqn.cons_to_diag = gkyl_default_cons_to_diag;

  iso_euler->eqn.source_func = iso_euler_source;

  iso_euler->eqn.ref_count = gkyl_ref_count_init(iso_euler_free);

  return &iso_euler->eqn;
}

double
gkyl_wv_iso_euler_vt(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_iso_euler *iso_euler = container_of(eqn, struct wv_iso_euler, eqn);
  return iso_euler->vt;
}
