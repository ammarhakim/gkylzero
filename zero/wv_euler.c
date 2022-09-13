#include <math.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_moment_prim_euler.h>
#include <gkyl_wv_euler.h>

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

static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<5; ++i)
    wout[i] = qin[i];
}

// Euler perfectly reflecting wall
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

// Euler no-slip wall
static void
euler_no_slip(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  // copy density and pressure
  ghost[0] = skin[0];
  ghost[4] = skin[4];

  // zero-normal for momentum
  ghost[1] = -skin[1];
  ghost[2] = -skin[2];
  ghost[3] = -skin[3];
}

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

// Waves and speeds using Lax fluxes
static double
wave_lax(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double gas_gamma = euler->gas_gamma;

  double rhol = ql[0], rhor = qr[0];
  double ul = ql[1]/ql[0], ur = qr[1]/qr[0];
  double pl = gkyl_euler_pressure(gas_gamma, ql), pr = gkyl_euler_pressure(gas_gamma, qr);
  double sl = fabs(ul) + sqrt(gas_gamma*pl/rhol), sr = fabs(ur) + sqrt(gas_gamma*pr/rhor);

  double *wv = &waves[0]; // single wave
  for (int i=0; i<5; ++i)  wv[i] = delta[i];

  s[0] = 0.5*(sl+sr);
  
  return s[0];
}

static void
qfluct_lax(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{

  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double gas_gamma = euler->gas_gamma;

  double rhol = ql[0], rhor = qr[0];
  double ul = ql[1]/ql[0], ur = qr[1]/qr[0];
  double pl = gkyl_euler_pressure(gas_gamma, ql), pr = gkyl_euler_pressure(gas_gamma, qr);
  double sl = fabs(ul) + sqrt(gas_gamma*pl/rhol), sr = fabs(ur) + sqrt(gas_gamma*pr/rhor);
  double amax = fmax(sl, sr);

  double fl[5], fr[5];
  gkyl_euler_flux(gas_gamma, ql, fl);
  gkyl_euler_flux(gas_gamma, qr, fr);

  for (int i=0; i<5; ++i) {
    amdq[i] = 0.5*(fr[i]-fl[i] - amax*(qr[i]-ql[i]));
    apdq[i] = 0.5*(fr[i]-fl[i] + amax*(qr[i]-ql[i]));
  }
}

static double
wave_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
    return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

// Waves and speeds using Roe averaging
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

static double
wave_roe_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return wave_roe(eqn, delta, ql, qr, waves, s);
  else
    return wave_lax(eqn, delta, ql, qr, waves, s);

  return 0.0; // can't happen
}

static void
qfluct_roe_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return qfluct_roe(eqn, ql, qr, waves, s, amdq, apdq);
  else
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

// HLLC
static double
wave_hllc(const struct gkyl_wv_eqn *eqn, const double *dQ, const double *ql,
  const double *qr, double *waves, double *speeds)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double g = euler->gas_gamma;

  // r->rho, u->u_x, p->pressure, c->sound speed; l->left, r->right, m->middle
  double rl = ql[0];
  double ul = ql[1] / rl;
  double pl = gkyl_euler_pressure(g, ql);
  double cl = sqrt(g*pl/rl);

  double rr = qr[0];
  double ur = qr[1] / rr;
  double pr = gkyl_euler_pressure(g, qr);
  double cr = sqrt(g*pr/rr);

  // STEP 1. compute min and max wave speeds; here, Toro (10.59) is chosen
  //   first, estimate middle pressure, Toro (10.61)
  double ra = 0.5 * (rl + rr);
  double ca = 0.5 * (cl + cr);
  double pm = 0.5 * (pl + pr) - 0.5 * (ur - ul) * ra * ca;
  //   second, compue the q coefficients in Toro (10.60)
  double coeffl = pm <= pl ? 1 : sqrt(1 + 0.5 * (1 + 1/g) * (pm/pl - 1));
  double coeffr = pm <= pr ? 1 : sqrt(1 + 0.5 * (1 + 1/g) * (pm/pr - 1));
  //   finally, compute Toro (10.59)
  double sl = ul - cl * coeffl;
  double sr = ur + cr * coeffr;

  // STEP 2. compute middle wave speed, Toro (10.37)
  double sm = (pr-pl+rl*ul*(sl-ul)-rr*ur*(sr-ur)) / (rl*(sl-ul)-rr*(sr-ur));

  // STEP 3. compute left and right intermediate states, Toro (10.39)
  double qml[5], qmr[5];

  qml[0] = rl * (sl-ul) / (sl-sm);
  qml[1] = qml[0] * sm;
  qml[2] = qml[0] * ql[2] / ql[0];
  qml[3] = qml[0] * ql[3] / ql[0];
  qml[4] = qml[0] * (ql[4]/rl + (sm-ul) * (sm + pl / rl / (sl-ul)));

  qmr[0] = rr * (sr-ur) / (sr-sm);
  qmr[1] = qmr[0] * sm;
  qmr[2] = qmr[0] * qr[2] / qr[0];
  qmr[3] = qmr[0] * qr[3] / qr[0];
  qmr[4] = qmr[0] * (qr[4]/rr + (sm-ur) * (sm + pr / rr / (sr-ur)));

  // STEP 4. collect all waves and speeds
  double *wv;

  wv = waves;
  speeds[0] = sl;
  for (int i=0; i<5; ++i)  wv[i] = qml[i] - ql[i];

  wv += 5;
  speeds[1] = sm;
  for (int i=0; i<5; ++i)  wv[i] = qmr[i] - qml[i];

  wv += 5;
  speeds[2] = sr;
  for (int i=0; i<5; ++i)  wv[i] = qr[i] - qmr[i];

  return sr;
}

static void
qfluct_hllc(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  const double *waves, const double *s, double *amdq, double *apdq)
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
wave_hllc_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return wave_hllc(eqn, delta, ql, qr, waves, s);
  else
    return wave_lax(eqn, delta, ql, qr, waves, s);

  return 0.0; // can't happen
}

static void
qfluct_hllc_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return qfluct_hllc(eqn, ql, qr, waves, s, amdq, apdq);
  else
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}


static bool
check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  
  if (q[0] < 0.0)
    return false;

  double pr = gkyl_euler_pressure(euler->gas_gamma, q);
  if (pr < 0.0)
    return false;

  return true;
}

static double
max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  return gkyl_euler_max_abs_speed(euler->gas_gamma, q);
}

struct gkyl_wv_eqn*
gkyl_wv_euler_inew(const struct gkyl_wv_euler_inp *inp)
{
  struct wv_euler *euler = gkyl_malloc(sizeof(struct wv_euler));

  euler->eqn.type = GKYL_EQN_EULER;
  euler->eqn.num_equations = 5;
  euler->gas_gamma = inp->gas_gamma;

  switch (inp->rp_type) {
    case WV_EULER_RP_ROE:
      euler->eqn.num_waves = 3;  
      euler->eqn.waves_func = wave_roe_l;
      euler->eqn.qfluct_func = qfluct_roe_l;
      break;

    case WV_EULER_RP_HLLC:
      euler->eqn.num_waves = 3;  
      euler->eqn.waves_func = wave_hllc_l;
      euler->eqn.qfluct_func = qfluct_hllc_l;
      break;
      
    case WV_EULER_RP_LAX:
      euler->eqn.num_waves = 1;
      euler->eqn.waves_func = wave_lax_l;
      euler->eqn.qfluct_func = qfluct_lax_l;
      break;      
  }
      
  
  euler->eqn.max_speed_func = max_speed;

  euler->eqn.rotate_to_local_func = rot_to_local;
  euler->eqn.rotate_to_global_func = rot_to_global;

  euler->eqn.check_inv_func = check_inv;

  euler->eqn.wall_bc_func = euler_wall;
  euler->eqn.no_slip_bc_func = euler_no_slip;

  euler->eqn.cons_to_riem = cons_to_riem;

  euler->eqn.ref_count = gkyl_ref_count_init(euler_free);

  return &euler->eqn;  
}

struct gkyl_wv_eqn*
gkyl_wv_euler_new(double gas_gamma)
{
  return gkyl_wv_euler_inew( &(struct gkyl_wv_euler_inp) {
      .gas_gamma = gas_gamma,
      .rp_type = WV_EULER_RP_ROE
    }
  );
}

double
gkyl_wv_euler_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  return euler->gas_gamma;
}
