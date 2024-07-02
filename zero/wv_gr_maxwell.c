#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_maxwell.h>
#include <gkyl_wv_gr_maxwell_priv.h>

static inline double
gkyl_gr_maxwell_max_abs_speed(double light_speed, const double q[11])
{
  double spatial_det = q[8];
  double lapse = q[9];

  bool in_excision_region = false;
  if (q[10] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    return light_speed / (lapse * sqrt(spatial_det));
  }
  else {
    return light_speed;
  }
}

static void
gkyl_gr_maxwell_flux(double light_speed, double e_fact, double b_fact, const double q[11], double flux[11])
{
  double Ex = q[0], Ey = q[1], Ez = q[2];
  double Bx = q[3], By = q[4], Bz = q[5];
  double phi = q[6];
  double psi = q[7];
  double spatial_det = q[8];
  double lapse = q[9];

  bool in_excision_region = false;
  if (q[10] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    flux[0] = (1.0 / ((lapse * lapse) * spatial_det)) * e_fact * (light_speed * light_speed) * phi;
    flux[1] = (1.0 / ((lapse * lapse) * spatial_det)) * (light_speed * light_speed) * Bz;
    flux[2] = -(1.0 / ((lapse * lapse) * spatial_det)) * (light_speed * light_speed) * By;
    flux[3] = b_fact * psi;
    flux[4] = -Ez;
    flux[5] = Ey;
    flux[6] = e_fact * Ex;
    flux[7] = (1.0 / ((lapse * lapse) * spatial_det)) * b_fact * (light_speed * light_speed) * Bx;

    for (int i = 8; i < 11; i++) {
      flux[i] = 0.0;
    }
  }
  else {
    for (int i = 0; i < 11; i++) {
      flux[i] = 0.0;
    }
  }
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 11; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 11; i++) {
    qout[i] = win[i];
  }
}

static void
gr_maxwell_wall(double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 11; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
  ghost[2] = -ghost[2];
  ghost[3] = -ghost[3];
  ghost[6] = -ghost[6];
}

static inline void
rot_to_local(const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal, double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = (qglobal[0] * norm[0]) + (qglobal[1] * norm[1]) + (qglobal[2] * norm[2]);
  qlocal[1] = (qglobal[0] * tau1[0]) + (qglobal[1] * tau1[1]) + (qglobal[2] * tau1[2]);
  qlocal[2] = (qglobal[0] * tau2[0]) + (qglobal[1] * tau2[1]) + (qglobal[2] * tau2[2]);

  qlocal[3] = (qglobal[3] * norm[0]) + (qglobal[4] * norm[1]) + (qglobal[5] * norm[2]);
  qlocal[4] = (qglobal[3] * tau1[0]) + (qglobal[4] * tau1[1]) + (qglobal[5] * tau1[2]);
  qlocal[5] = (qglobal[3] * tau2[0]) + (qglobal[4] * tau2[1]) + (qglobal[5] * tau2[2]);

  qlocal[6] = qglobal[6];
  qlocal[7] = qglobal[7];
  qlocal[8] = qglobal[8];
  qlocal[9] = qglobal[9];
  qlocal[10] = qglobal[10];
}

static inline void
rot_to_global(const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal, double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = (qlocal[0] * norm[0]) + (qlocal[1] * tau1[0]) + (qlocal[2] * tau2[0]);
  qglobal[1] = (qlocal[0] * norm[1]) + (qlocal[1] * tau1[1]) + (qlocal[2] * tau2[1]);
  qglobal[2] = (qlocal[0] * norm[2]) + (qlocal[1] * tau1[2]) + (qlocal[2] * tau2[2]);

  qglobal[3] = (qlocal[3] * norm[0]) + (qlocal[4] * tau1[0]) + (qlocal[5] * tau2[0]);
  qglobal[4] = (qlocal[3] * norm[1]) + (qlocal[4] * tau1[1]) + (qlocal[5] * tau2[1]);
  qglobal[5] = (qlocal[3] * norm[2]) + (qlocal[4] * tau1[2]) + (qlocal[5] * tau2[2]);

  qglobal[6] = qlocal[6];
  qglobal[7] = qlocal[7];
  qglobal[8] = qlocal[8];
  qglobal[9] = qlocal[9];
  qglobal[10] = qlocal[10];
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
  double light_speed = gr_maxwell->light_speed;
  double e_fact = gr_maxwell->e_fact;
  double b_fact = gr_maxwell->b_fact;

  double sl = gkyl_gr_maxwell_max_abs_speed(light_speed, ql);
  double sr = gkyl_gr_maxwell_max_abs_speed(light_speed, qr);
  double amax = fmax(sl, sr);

  double fl[11], fr[11];
  gkyl_gr_maxwell_flux(light_speed, e_fact, b_fact, ql, fl);
  gkyl_gr_maxwell_flux(light_speed, e_fact, b_fact, qr, fr);

  bool in_excision_region_l = false;
  if (ql[10] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[10] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  double *w0 = &waves[0], *w1 = &waves[11];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 11; i++) {
      w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
      w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
    }
  }
  else {
    for (int i = 0; i < 11; i++) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }

  s[0] = -amax;
  s[1] = amax;

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double* w0 = &waves[0], *w1 = &waves[11];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 11; i++) {
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
qfluct_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
wave_roe(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
  double light_speed = gr_maxwell->light_speed;
  double e_fact = gr_maxwell->e_fact;
  double b_fact = gr_maxwell->b_fact;

  double spatial_det_l = ql[8];
  double lapse_l = ql[9];

  bool in_excision_region_l = false;
  if (ql[10] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  double spatial_det_r = qr[8];
  double lapse_r = qr[9];

  bool in_excision_region_r = false;
  if (qr[10] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  double spatial_det = (0.5 * spatial_det_l) + (0.5 * spatial_det_r);
  double lapse = (0.5 * lapse_l) + (0.5 * lapse_r);

  double a1 = 0.5 * (delta[3] - ((lapse * sqrt(spatial_det)) / light_speed) * delta[7]);
  double a2 = 0.5 * (delta[3] + ((lapse * sqrt(spatial_det)) / light_speed) * delta[7]);
  double a3 = 0.5 * (delta[0] - (light_speed / (lapse * sqrt(spatial_det))) * delta[6]);
  double a4 = 0.5 * (delta[0] + (light_speed / (lapse * sqrt(spatial_det))) * delta[6]);
  double a5 = 0.5 * (delta[1] - (light_speed / (lapse * sqrt(spatial_det))) * delta[5]);
  double a6 = 0.5 * ((light_speed / (lapse * sqrt(spatial_det))) * delta[4] + delta[2]);
  double a7 = 0.5 * ((light_speed / (lapse * sqrt(spatial_det))) * delta[5] + delta[1]);
  double a8 = 0.5 * (delta[2] - (light_speed / (lapse * sqrt(spatial_det))) * delta[4]);

  for (int i = 0; i < (11 * 6); i++) {
    waves[i] = 0.0;
  }

  double *wv = 0;

  if (!in_excision_region_l && !in_excision_region_r) {
    wv = &waves[0 * 11];
    wv[3] = a1;
    wv[7] = -a1 * (light_speed / (lapse * sqrt(spatial_det)));
    s[0] = -(light_speed / (lapse * sqrt(spatial_det))) * b_fact;

    wv = &waves[1 * 11];
    wv[3] = a2;
    wv[7] = a2 * (light_speed / (lapse * sqrt(spatial_det)));
    s[1] = (light_speed / (lapse * sqrt(spatial_det))) * b_fact;

    wv = &waves[2 * 11];
    wv[0] = a3;
    wv[6] = -a3 * ((lapse * sqrt(spatial_det)) / light_speed);
    s[2] = -(light_speed / (lapse * sqrt(spatial_det))) * e_fact;

    wv = &waves[3 * 11];
    wv[0] = a4;
    wv[6] = a4 * ((lapse * sqrt(spatial_det)) / light_speed);
    s[3] = (light_speed / (lapse * sqrt(spatial_det))) * e_fact;

    wv = &waves[4 * 11];
    wv[1] = a5;
    wv[2] = a6;
    wv[4] = a6 * ((lapse * sqrt(spatial_det)) / light_speed);
    wv[5] = -a5 * ((lapse * sqrt(spatial_det)) / light_speed);
    s[4] = -(light_speed / (lapse * sqrt(spatial_det)));

    wv = &waves[5 * 11];
    wv[1] = a7;
    wv[2] = a8;
    wv[4] = -a8 * ((lapse * sqrt(spatial_det)) / light_speed);
    wv[5] = a7 * ((lapse * sqrt(spatial_det)) / light_speed);
    s[5] = light_speed / (lapse * sqrt(spatial_det));

    return light_speed / (lapse * sqrt(spatial_det));
  }
  else {
    for (int i = 0; i < 5; i++) {
      s[0] = pow(10.0, -8.0);
    }

    return pow(10.0, -8.0);
  }
}

static void
qfluct_roe(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0 * 11], *w1 = &waves[1 * 11], *w2 = &waves[2 * 11];
  const double *w3 = &waves[3 * 11], *w4 = &waves[4 * 11], *w5 = &waves[5 * 11];

  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s3m = fmin(0.0, s[3]), s4m = fmin(0.0, s[4]), s5m = fmin(0.0, s[5]);

  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);
  double s3p = fmax(0.0, s[3]), s4p = fmax(0.0, s[4]), s5p = fmax(0.0, s[5]);

  for (int i = 0; i < 8; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]) + (s2m * w2[i]) + (s3m * w3[i]) + (s4m * w4[i]) + (s5m * w5[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]) + (s2p * w2[i]) + (s3p * w3[i]) + (s4p * w4[i]) + (s5p * w5[i]);
  }

  for (int i = 8; i < 11; i++) {
    amdq[i] = 0.0;
    apdq[i] = 0.0;
  }
}

static double
wave_roe_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return wave_roe(eqn, delta, ql, qr, waves, s);
  }
  else {
    return wave_lax(eqn, delta, ql, qr, waves, s);
  }
}

static void
qfluct_roe_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return qfluct_roe(eqn, ql, qr, waves, s, amdq, apdq);
  }
  else {
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
  }
}

static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);

  double fr[11], fl[11];
  gkyl_gr_maxwell_flux(gr_maxwell->light_speed, gr_maxwell->e_fact, gr_maxwell->b_fact, ql, fl);
  gkyl_gr_maxwell_flux(gr_maxwell->light_speed, gr_maxwell->e_fact, gr_maxwell->b_fact, qr, fr);

  bool in_excision_region_l = false;
  if (ql[10] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[10] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  if (!in_excision_region_l && !in_excision_region_r) {
    for (int m = 0; m < 11; m++) {
      flux_jump[m] = fr[m] - fl[m];
    }
  }
  else {
    for (int m = 0; m < 11; m++) {
      flux_jump[m] = 0.0;
    }
  }

  double amaxl = gkyl_gr_maxwell_max_abs_speed(gr_maxwell->light_speed, ql);
  double amaxr = gkyl_gr_maxwell_max_abs_speed(gr_maxwell->light_speed, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  return true;
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
  double light_speed = gr_maxwell -> light_speed;

  return gkyl_gr_maxwell_max_abs_speed(light_speed, q);
}

static inline void
gr_maxwell_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 6; i++) {
    diag[i] = qin[i] * qin[i];
  }
}

void
gkyl_gr_maxwell_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object,
    struct wv_gr_maxwell *gr_maxwell = container_of(base->on_dev, struct wv_gr_maxwell, eqn);
    gkyl_cu_free(gr_maxwell);
  }

  struct wv_gr_maxwell *gr_maxwell = container_of(base, struct wv_gr_maxwell, eqn);
  gkyl_free(gr_maxwell);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_maxwell_new(double light_speed, double e_fact, double b_fact, struct gkyl_gr_spacetime* spacetime, bool use_gpu)
{
  return gkyl_wv_gr_maxwell_inew(&(struct gkyl_wv_gr_maxwell_inp) {
      .light_speed = light_speed,
      .e_fact = e_fact,
      .b_fact = b_fact,
      .spacetime = spacetime,
      .rp_type = WV_GR_MAXWELL_RP_ROE,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_maxwell_inew(const struct gkyl_wv_gr_maxwell_inp* inp)
{
  struct wv_gr_maxwell *gr_maxwell = gkyl_malloc(sizeof(struct wv_gr_maxwell));

  gr_maxwell->eqn.type = GKYL_EQN_GR_MAXWELL;
  gr_maxwell->eqn.num_equations = 11;
  gr_maxwell->eqn.num_diag = 6;

  gr_maxwell->light_speed = inp->light_speed;
  gr_maxwell->e_fact = inp->e_fact;
  gr_maxwell->b_fact = inp->b_fact;
  gr_maxwell->spacetime = inp->spacetime;

  if (inp->rp_type == WV_GR_MAXWELL_RP_LAX) {
    gr_maxwell->eqn.num_waves = 2;
    gr_maxwell->eqn.waves_func = wave_lax_l;
    gr_maxwell->eqn.qfluct_func = qfluct_lax_l;
  }

  if (inp->rp_type == WV_GR_MAXWELL_RP_ROE) {
    gr_maxwell->eqn.num_waves = 6;
    gr_maxwell->eqn.waves_func = wave_roe_l;
    gr_maxwell->eqn.qfluct_func = qfluct_roe_l;
  }

  gr_maxwell->eqn.flux_jump = flux_jump;
  gr_maxwell->eqn.check_inv_func = check_inv;
  gr_maxwell->eqn.max_speed_func = max_speed;
  gr_maxwell->eqn.rotate_to_local_func = rot_to_local;
  gr_maxwell->eqn.rotate_to_global_func = rot_to_global;

  gr_maxwell->eqn.wall_bc_func = gr_maxwell_wall;

  gr_maxwell->eqn.cons_to_riem = cons_to_riem;
  gr_maxwell->eqn.riem_to_cons = riem_to_cons;

  gr_maxwell->eqn.cons_to_diag = gr_maxwell_cons_to_diag;

  gr_maxwell->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_maxwell->eqn.flags);
  gr_maxwell->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_maxwell_free);
  gr_maxwell->eqn.on_dev = &gr_maxwell->eqn; // On the CPU, the equation object points to itself.

  return &gr_maxwell->eqn;
}

double
gkyl_wv_gr_maxwell_light_speed(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
  double light_speed = gr_maxwell->light_speed;
  
  return light_speed;
}

double
gkyl_wv_gr_maxwell_e_fact(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
  double e_fact = gr_maxwell->e_fact;

  return e_fact;
}

double
gkyl_wv_gr_maxwell_b_fact(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
  double b_fact = gr_maxwell->b_fact;

  return b_fact;
}

struct gkyl_gr_spacetime*
gkyl_wv_gr_maxwell_spacetime(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
  struct gkyl_gr_spacetime* spacetime = gr_maxwell->spacetime;

  return spacetime;
}