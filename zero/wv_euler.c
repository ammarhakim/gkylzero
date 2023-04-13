#include <math.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_moment_prim_euler.h>
#include <gkyl_wv_euler.h>

enum hll_speed_type {
  HLL_SPEED_PVRS = 0, // default
  HLL_SPEED_MINMAX = 1,
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

static inline void
cons_to_riem(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<5; ++i)
    wout[i] = qin[i];
}
static inline void
riem_to_cons(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout)
{
  // TODO: this should use proper L matrix
  for (int i=0; i<5; ++i)
    qout[i] = win[i];
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
  double amax = fmax(sl, sr);

  double fl[5], fr[5];
  gkyl_euler_flux(gas_gamma, ql, fl);
  gkyl_euler_flux(gas_gamma, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[5];
  for (int i=0; i<5; ++i) {
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
  const double *w0 = &waves[0], *w1 = &waves[5];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i=0; i<5; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i];
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

// project column vector delta onto the right eigenvectors of the flux Jacobian
// evaluated at avg = {u, v, w, enth}
static double
proj_onto_euler_eigvect(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *avg, double *waves, double *s)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double gas_gamma = euler->gas_gamma;
  double g1 = gas_gamma - 1;

  double u = avg[0], v = avg[1], w = avg[2], enth = avg[3];

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

inline static void
roe_avg(const struct gkyl_wv_eqn *eqn, const double *ql,
    const double *qr, double *avg)
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

  avg[0] = u; avg[1] = v; avg[2] = w; avg[3] = enth;
}

// Waves and speeds using Roe averaging
static double
wave_roe(const struct gkyl_wv_eqn *eqn,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  double avg[4];
  roe_avg(eqn, ql, qr, avg);

  return proj_onto_euler_eigvect(eqn, delta, avg, waves, s);
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

// HLL
static void
states_hll_common(const struct gkyl_wv_eqn *eqn, const double *ql,
  const double *qr, double state[8], const enum hll_speed_type speed_type)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double g = euler->gas_gamma;

  // r->rho, u->u_x, p->pressure, c->sound speed; l->left, r->right, m->middle
  double rl = ql[0];
  double ul = ql[1] / rl;
  double pl = gkyl_euler_pressure(g, ql);

  double rr = qr[0];
  double ur = qr[1] / rr;
  double pr = gkyl_euler_pressure(g, qr);

  // STEP 1. compute min and max wave speeds
  double sl, sr;
  if (speed_type==HLL_SPEED_PVRS) {

    // using the pressure-based estimates of Toro (10.59)
    //   first, estimate middle pressure; using PVRS Toro (10.61) but Toro Fig 9.4
    //   presents an adaptive method to switch between PVRS and TRRS (10.63) and
    //   TSRS (10.65); for ideal gases, max(0, p_PVRS) might also be good
    double cl = sqrt(g*pl/rl);
    double cr = sqrt(g*pr/rr);
    double ra = 0.5 * (rl + rr);
    double ca = 0.5 * (cl + cr);
    double pm = 0.5 * (pl + pr) - 0.5 * (ur - ul) * ra * ca;
    //   second, compue the q coefficients in Toro (10.60)
    double coeffl = pm <= pl ? 1 : sqrt(1 + 0.5 * (1 + 1/g) * (pm/pl - 1));
    double coeffr = pm <= pr ? 1 : sqrt(1 + 0.5 * (1 + 1/g) * (pm/pr - 1));
    //   finally, compute Toro (10.59)
    sl = ul - cl * coeffl;
    sr = ur + cr * coeffr;

  } else if (speed_type==HLL_SPEED_MINMAX) {

    double cl = sqrt(g*pl/rl);
    double cr = sqrt(g*pr/rr);

    double qm[5];
    for (int c=0; c<5; ++c)
      qm[c] = 0.5 * (ql[c] + qr[c]);

    double rm = qm[0];
    double um = qm[1] / rm;
    double pm = gkyl_euler_pressure(g, qm);
    double cm = sqrt(g*pm/rm);

    /* double sl = fmin(ul-cl, fmin(ur-cr, um-cm)); */
    /* double sr = fmax(ul+cl, fmax(ur+cr, um+cm)); */
    sl = fmin(ul-cl, um-cm);
    sr = fmax(ur+cr, um+cm);

    double avg[4];
    roe_avg(eqn, ql, qr, avg);
    double u_roe=avg[0], v_roe=avg[1], w_roe=avg[2], enth_roe=avg[3];
    double q2_roe = u_roe*u_roe + v_roe*v_roe + w_roe*w_roe;
    double c_roe = sqrt((g-1) * (enth_roe - 0.5*q2_roe));

    sl = fmin(sl, u_roe-c_roe);
    sr = fmax(sr, u_roe+c_roe);
  }

  state[0] = rl;
  state[1] = ul;
  state[2] = pl;
  state[3] = rr;
  state[4] = ur;
  state[5] = pr;
  state[6] = sl;
  state[7] = sr;
}

// HLL
static void
states_hll(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  double *speeds, double *qm, const enum hll_speed_type speed_type)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);

  // STEP 1. compute min and max wave speeds
  double state[8];
  states_hll_common(eqn, ql, qr, state, speed_type);
  double rl = state[0], ul = state[1], pl = state[2];
  double rr = state[3], ur = state[4], pr = state[5];
  double sl = state[6], sr = state[7];

  // STEP 2. compute left and right fluxes
  double fl[5];
  fl[0] = ql[0] * ul;
  fl[1] = ql[1] * ul + pl;
  fl[2] = ql[2] * ul;
  fl[3] = ql[3] * ul;
  fl[4] = ql[4] * ul + pl * ul;

  double fr[5];
  fr[0] = qr[0] * ur;
  fr[1] = qr[1] * ur + pr;
  fr[2] = qr[2] * ur;
  fr[3] = qr[3] * ur;
  fr[4] = qr[4] * ur + pr * ur;

  // STEP 3. compute standard HLL intermediate states
  qm[0] = (sr*qr[0]-sl*ql[0]+fl[0]-fr[0]) / (sr-sl);
  qm[1] = (sr*qr[1]-sl*ql[1]+fl[1]-fr[1]) / (sr-sl);
  qm[2] = (sr*qr[2]-sl*ql[2]+fl[2]-fr[2]) / (sr-sl);
  qm[3] = (sr*qr[3]-sl*ql[3]+fl[3]-fr[3]) / (sr-sl);
  qm[4] = (sr*qr[4]-sl*ql[4]+fl[4]-fr[4]) / (sr-sl);

  // STEP 4. collect all speeds
  speeds[0] = sl;
  speeds[1] = sr;
}
  
static double
wave_hll(const struct gkyl_wv_eqn *eqn, const double *dQ, const double *ql,
  const double *qr, double *waves, double *speeds)
{
  double qm[5];
  states_hll(eqn, ql, qr, speeds, qm, HLL_SPEED_PVRS);

  double *wv;

  wv = waves;
  for (int i=0; i<5; ++i)  wv[i] = qm[i] - ql[i];

  wv += 5;
  for (int i=0; i<5; ++i)  wv[i] = qr[i] - qm[i];

  return fmax(fabs(speeds[0]), fabs(speeds[1]));
}

static void
qfluct_hll(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  const double *waves, const double *s, double *amdq, double *apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[5];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i=0; i<5; ++i) {
    amdq[i] = s0m*w0[i] + s1m*w1[i];
    apdq[i] = s0p*w0[i] + s1p*w1[i];
  }
}

static double
wave_hll_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return wave_hll(eqn, delta, ql, qr, waves, s);
  else // FIXME perhaps not needed
    return wave_lax(eqn, delta, ql, qr, waves, s);

  return 0.0; // can't happen
}

static void
qfluct_hll_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return qfluct_hll(eqn, ql, qr, waves, s, amdq, apdq);
  else
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

// HLLI

// XXX: Following implementation is ad hoc. Original HLLI algorithm corrected
// the qfluct instead of waves, as the wave corrections are dependent on the
// self-similarity variable, i.e., not a constant. However, in our tests
// directly correcting qfluct causes inconsistency between the HLLI qfluct and
// the HLL waves used for limiting. Following algorithm applies a plausible
// correction to the waves instead, so that the waves and the qfluct are
// consistent. However, the results with 2nd-order limiting/correction is much
// more diffusive than HLLC/Roe (though better than HLL).
static double
wave_hlli(const struct gkyl_wv_eqn *eqn, const double *dQ, const double *ql,
  const double *qr, double *waves, double *speeds)
{
  double qm[5];
  states_hll(eqn, ql, qr, speeds, qm, HLL_SPEED_MINMAX);

  double *wv;

  wv = waves;
  for (int i=0; i<5; ++i)  wv[i] = qm[i] - ql[i];

  wv += 5;
  for (int i=0; i<5; ++i)  wv[i] = qr[i] - qm[i];

  double sl = speeds[0], sr = speeds[1];
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double my_jump[5], my_waves[3*5], my_s[3];

  // 1. Compute jumps
  for (int c=0; c<5; ++c)
    my_jump[c] = qr[c] - ql[c];

  // 2. Project jumps onto eigenvectors to get waves R*(L*(QR-QL)) in eq (28)
  int avg_type = 0;
  if (avg_type==0) { // using eigenstructure due to arithmetic averages

    double qm[5];
    for (int c=0; c<5; ++c) {
      qm[c] = 0.5 * (ql[c] + qr[c]);
    }

    double u=qm[1]/qm[0], v=qm[2]/qm[0], w=qm[3]/qm[0];
    double pr = gkyl_euler_pressure(euler->gas_gamma, qm);
    double enth = (qm[4] + pr) / qm[0];
    double avg[4] = {u, v, w, enth};

    proj_onto_euler_eigvect(eqn, my_jump, avg, my_waves, my_s);

  } else if (avg_type==1) { // using eigenstructure due to Roe averages
    wave_roe(eqn, my_jump, ql, qr, my_waves, my_s);
  }

  // 3. Compute the delta^p weights in eq (29)
  double del[3];
  for (int p=0; p<3; ++p) { // TODO zero sl/sr; upper/lower limits
    del[p] = 1 - fmin(my_s[p],0)/sl - fmax(my_s[p],0)/sr;
    del[p] = fmax(fmin(del[p], 1), 0);
  }

  // 4. Accumulate HLLI corrections to the flucuations as in eq (28)
  double corr[5] = {0};
  for (int p=0; p<3; ++p) { // FIXME apply on intermediate waves only?
    double *wv = my_waves + p*5;
    for (int c=0; c<5; ++c) {
      corr[c] += (1 / (sr - sl)) * del[p] * wv[c];
    }
  }

  // 4. APPLY HLLI corrections to the fluctuations as in eq (28)
  double *wl = waves;
  double *wr = waves + 5;
  for (int c=0; c<5; ++c) {
    wl[c] -= corr[c] * sr;
    wr[c] += corr[c] * sl;
  }

  return fmax(fabs(speeds[0]), fabs(speeds[1]));
}

static double
wave_hlli_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double *waves, double *s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX)
    return wave_hlli(eqn, delta, ql, qr, waves, s);
  else // FIXME perhaps not needed
    return wave_lax(eqn, delta, ql, qr, waves, s);

  return 0.0; // can't happen
}

static void
qfluct_hlli_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *s,
  double *amdq, double *apdq)
{
  qfluct_hll_l(eqn, type, ql, qr, waves, s, amdq, apdq);
}

// HLLC
static void
states_hllc(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  double *speeds, double *qml, double *qmr)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);

  // STEP 1. compute min and max wave speeds
  double state[8];
  states_hll_common(eqn, ql, qr, state, HLL_SPEED_PVRS);
  double rl = state[0], ul = state[1], pl = state[2];
  double rr = state[3], ur = state[4], pr = state[5];
  double sl = state[6], sr = state[7];

  // STEP 2. compute middle wave speed, Toro (10.37)
  double sm = (pr-pl+rl*ul*(sl-ul)-rr*ur*(sr-ur)) / (rl*(sl-ul)-rr*(sr-ur));

  // STEP 3. compute left and right intermediate states, Toro (10.39)
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

  // STEP 4. collect all speeds
  speeds[0] = sl;
  speeds[1] = sm;
  speeds[2] = sr;
}
  
static double
wave_hllc(const struct gkyl_wv_eqn *eqn, const double *dQ, const double *ql,
  const double *qr, double *waves, double *speeds)
{
  double qml[5], qmr[5];
  states_hllc(eqn, ql, qr, speeds, qml, qmr);

  double *wv;

  wv = waves;
  for (int i=0; i<5; ++i)  wv[i] = qml[i] - ql[i];

  wv += 5;
  for (int i=0; i<5; ++i)  wv[i] = qmr[i] - qml[i];

  wv += 5;
  for (int i=0; i<5; ++i)  wv[i] = qr[i] - qmr[i];

  return fmax(fabs(speeds[0]), fabs(speeds[2]));
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

static void
qfluct_hllc_direct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *speeds,
  double *amdq, double *apdq)
{
  double s[3], qml[5], qmr[5];
  states_hllc(eqn, ql, qr, s, qml, qmr);

  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i=0; i<5; ++i) {
    double w0 = qml[i] - ql[i];
    double w1 = qmr[i] - qml[i];
    double w2 = qr[i] - qmr[i];
    amdq[i] = s0m*w0 + s1m*w1 + s2m*w2;
    apdq[i] = s0p*w0 + s1p*w1 + s2p*w2;
  }
}

static double
flux_jump(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);

  double fr[5], fl[5];
  gkyl_euler_flux(euler->gas_gamma, ql, fl);
  gkyl_euler_flux(euler->gas_gamma, qr, fr);

  for (int m=0; m<5; ++m) flux_jump[m] = fr[m]-fl[m];

  double amaxl = gkyl_euler_max_abs_speed(euler->gas_gamma, ql);
  double amaxr = gkyl_euler_max_abs_speed(euler->gas_gamma, qr);  

  return fmax(amaxl, amaxr);
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

static inline void
euler_cons_to_diag(const struct gkyl_wv_eqn *eqn,
  const double *qin, double *diag)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  // density and moment as copied as-is
  for (int i=0; i<4; ++i) diag[i] = qin[i];
  double ke = 0.5*(qin[1]*qin[1] + qin[2]*qin[2] + qin[3]*qin[3])/qin[0];
  diag[4] = ke; 
  diag[5] = qin[4]-ke;
}

struct gkyl_wv_eqn*
gkyl_wv_euler_inew(const struct gkyl_wv_euler_inp *inp)
{
  struct wv_euler *euler = gkyl_malloc(sizeof(struct wv_euler));

  euler->eqn.type = GKYL_EQN_EULER;
  euler->eqn.num_equations = 5;
  euler->eqn.num_diag = 6; // KE and PE stored separate
  
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
      euler->eqn.num_waves = 2;
      euler->eqn.waves_func = wave_lax_l;
      euler->eqn.qfluct_func = qfluct_lax_l;
      break;   

    case WV_EULER_RP_HLL:
      euler->eqn.num_waves = 2;  
      euler->eqn.waves_func = wave_hll_l;
      euler->eqn.qfluct_func = qfluct_hll_l;
      break;

    case WV_EULER_RP_HLLI:
      euler->eqn.num_waves = 2;
      euler->eqn.waves_func = wave_hlli_l;
      euler->eqn.qfluct_func = qfluct_hlli_l;
      break;

  }
      

  euler->eqn.flux_jump = flux_jump;
  euler->eqn.max_speed_func = max_speed;

  euler->eqn.rotate_to_local_func = rot_to_local;
  euler->eqn.rotate_to_global_func = rot_to_global;

  euler->eqn.check_inv_func = check_inv;

  euler->eqn.wall_bc_func = euler_wall;
  euler->eqn.no_slip_bc_func = euler_no_slip;

  euler->eqn.cons_to_riem = cons_to_riem;
  euler->eqn.riem_to_cons = riem_to_cons;

  euler->eqn.cons_to_diag = euler_cons_to_diag;

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
