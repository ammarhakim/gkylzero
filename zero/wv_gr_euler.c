#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h>

static inline void
gkyl_gr_euler_prim_vars(double gas_gamma, const double q[5], double v[5])
{
  double D = q[0];
  double momx = q[1];
  double momy = q[2];
  double momz = q[3];
  double Etot = q[4];

  double C = D / sqrt(((Etot + D) * (Etot + D)) - ((momx * momx) + (momy * momy) + (momz * momz)));
  double C0 = (D + Etot) / sqrt(((Etot + D) * (Etot + D)) - ((momx * momx) + (momy * momy) + (momz * momz)));
  if (((Etot + D) * (Etot + D)) - ((momx * momx) + (momy * momy) + (momz * momz)) < pow(10.0, -8.0)) {
    C = D / sqrt(pow(10.0, -8.0));
    C0 = (D + Etot) / sqrt(pow(10.0, -8.0));
  }

  double alpha0 = -1.0 / (gas_gamma * gas_gamma);
  double alpha1 = -2.0 * C * ((gas_gamma - 1.0) / (gas_gamma * gas_gamma));
  double alpha2 = ((gas_gamma - 2.0) / gas_gamma) * ((C0 * C0) - 1.0) + 1.0 - (C * C) * ((gas_gamma - 1.0) / gas_gamma) * ((gas_gamma - 1.0) / gas_gamma);
  double alpha4 = (C0 * C0) - 1.0;
  double eta = 2.0 * C *((gas_gamma - 1.0) / gas_gamma);

  double guess = 1.0;
  int iter = 0;

  while (iter < 1000) {
    double poly = (alpha4 * (guess * guess * guess) * (guess - eta)) + (alpha2 * (guess * guess)) + (alpha1 * guess) + alpha0;
    double poly_der = alpha1 + (2.0 * alpha2 * guess) + (4.0 * alpha4 * (guess * guess * guess)) - (3.0 * eta * alpha4 * (guess * guess));

    double guess_new = guess - (poly / poly_der);

    if (fabs(guess - guess_new) < pow(10.0, -8.0)) {
      iter = 1000;
    }
    else {
      iter += 1;
      guess = guess_new;
    }
  }

  double W = 0.5 * C0 * guess * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma - 1.0) / gas_gamma) * ((1.0 - (C * guess)) / ((C0 * C0) * (guess * guess))))));
  double h = 1.0 / (C * guess);

  v[0] = D / W;
  v[1] = momx / (v[0] * h * (W * W));
  v[2] = momy / (v[0] * h * (W * W));
  v[3] = momz / (v[0] * h * (W * W));
  v[4] = (v[0] * h * (W * W)) - D - Etot;

  if (v[0] < pow(10.0, -8.0)) {
    v[0] = pow(10.0, -8.0);
  }
  if (v[4] < pow(10.0, -8.0)) {
    v[4] = pow(10.0, -8.0);
  }
}

static inline double
gkyl_gr_euler_max_abs_speed(double gas_gamma, const double q[5])
{
  double v[5] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);
  double rho = v[0];
  double p = v[4];

  double num = (gas_gamma * p) / rho;
  double den = 1.0 * ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));
  double cs = sqrt(num / den);

  double v_sq = sqrt((v[1] * v[1]) + (v[2] * v[2]) + (v[3] * v[3]));
  return fabs(v_sq) + cs;
}

static void
gkyl_gr_euler_flux(double gas_gamma, const double q[5], double flux[5])
{
  double v[5] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);
  double rho =  v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double W = 1.0 / (sqrt(1.0 - ((vx * vx) + (vy * vy) + (vz * vz))));
  double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

  flux[0] = rho * W * vx;
  flux[1] = rho * h * (W * W) * (vx * vx) + p;
  flux[2] = rho * h * (W * W) * (vx * vy);
  flux[3] = rho * h * (W * W) * (vx * vz);
  flux[4] = ((rho * h * (W * W)) - p - (rho * W)) * vx + (p * vx);
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 5; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 5; i++) {
    qout[i] = win[i];
  }
}

static void
gr_euler_wall(double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 5; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

static void
gr_euler_no_slip(double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 1; i < 4; i++) {
    ghost[i] = -skin[i];
  }

  ghost[0] = skin[0];
  ghost[4] = skin[4];
}

static inline void
rot_to_local(const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal, double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);
  qlocal[4] = qglobal[4];
}

static inline void
rot_to_global(const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal, double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);
  qglobal[4] = qlocal[4];
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  double sl = gkyl_gr_euler_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_gr_euler_max_abs_speed(gas_gamma, qr);
  double amax = fmax(sl, sr);

  double fl[5], fr[5];
  gkyl_gr_euler_flux(gas_gamma, ql, fl);
  gkyl_gr_euler_flux(gas_gamma, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[5];
  for (int i = 0; i < 5; i++) {
    w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
    w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
  }

  s[0] = -amax;
  s[1] = amax;

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[5];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 5; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);

  double fr[5], fl[5];
  gkyl_gr_euler_flux(gr_euler->gas_gamma, ql, fl);
  gkyl_gr_euler_flux(gr_euler->gas_gamma, qr, fr);

  for (int m = 0; m < 5; m++) {
    flux_jump[m] = fr[m] - fl[m];
  }

  double amaxl = gkyl_gr_euler_max_abs_speed(gr_euler->gas_gamma, ql);
  double amaxr = gkyl_gr_euler_max_abs_speed(gr_euler->gas_gamma, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  double v[5] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);

  if (v[0] < 0.0 || v[4] < 0.0) {
    return false;
  }
  else {
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  return gkyl_gr_euler_max_abs_speed(gas_gamma, q);
}

static inline void
gr_euler_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 5; i++) {
    diag[i] = qin[i];
  }
}

void
gkyl_gr_euler_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_euler *gr_euler = container_of(base->on_dev, struct wv_gr_euler, eqn);
    gkyl_cu_free(gr_euler);
  }

  struct wv_gr_euler *gr_euler = container_of(base, struct wv_gr_euler, eqn);
  gkyl_free(gr_euler);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_new(double gas_gamma, bool use_gpu)
{
  return gkyl_wv_gr_euler_inew(&(struct gkyl_wv_gr_euler_inp) {
      .gas_gamma = gas_gamma,
      .rp_type = WV_GR_EULER_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_inew(const struct gkyl_wv_gr_euler_inp* inp)
{
  struct wv_gr_euler *gr_euler = gkyl_malloc(sizeof(struct wv_gr_euler));

  gr_euler->eqn.type = GKYL_EQN_GR_EULER;
  gr_euler->eqn.num_equations = 5;
  gr_euler->eqn.num_diag = 5;

  gr_euler->gas_gamma = inp->gas_gamma;

  if (inp->rp_type == WV_GR_EULER_RP_LAX) {
    gr_euler->eqn.num_waves = 2;
    gr_euler->eqn.waves_func = wave_lax_l;
    gr_euler->eqn.qfluct_func = qfluct_lax_l;
  }

  gr_euler->eqn.flux_jump = flux_jump;
  gr_euler->eqn.check_inv_func = check_inv;
  gr_euler->eqn.max_speed_func = max_speed;
  gr_euler->eqn.rotate_to_local_func = rot_to_local;
  gr_euler->eqn.rotate_to_global_func = rot_to_global;

  gr_euler->eqn.wall_bc_func = gr_euler_wall;
  gr_euler->eqn.no_slip_bc_func = gr_euler_no_slip;

  gr_euler->eqn.cons_to_riem = cons_to_riem;
  gr_euler->eqn.riem_to_cons = riem_to_cons;

  gr_euler->eqn.cons_to_diag = gr_euler_cons_to_diag;

  gr_euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_euler->eqn.flags);
  gr_euler->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_euler_free);
  gr_euler->eqn.on_dev = &gr_euler->eqn; // On the CPU, the equation object points to itself.

  return &gr_euler->eqn;
}

double
gkyl_wv_gr_euler_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  return gas_gamma;
}