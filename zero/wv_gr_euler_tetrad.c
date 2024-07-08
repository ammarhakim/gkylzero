#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h>
#include <gkyl_wv_gr_euler_tetrad.h>
#include <gkyl_wv_gr_euler_tetrad_priv.h>

static void
gkyl_gr_euler_tetrad_flux(double gas_gamma, const double q[29], double flux[29])
{
  double v[29] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);
  double rho =  v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  bool in_excision_region = false;
  if (v[28] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double W = 1.0 / (sqrt(1.0 - (vx * vx) - (vy * vy) - (vz * vz)));
    if ((vx * vx) + (vy * vy) + (vz * vz) > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

    flux[0] = rho * W * vx;
    flux[1] = (rho * h * (W * W) * (vx * vx)) + p;
    flux[2] = rho * h * (W * W) * (vy * vx);
    flux[3] = rho * h * (W * W) * (vz * vx);
    flux[4] = ((rho * h * (W * W)) - (rho * W)) * vx;

    for (int i = 5; i < 29; i++) {
      flux[i] = 0.0;
    }
  }
  else {
    for (int i = 0; i < 29; i++) {
      flux[i] = 0.0;
    }
  }
}

static void
gkyl_gr_euler_tetrad_extrapolate_flux(double gas_gamma, const double q[29], const double flux[29], double flux_extrap[29])
{
  double v[29];
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);
  double rho = v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double spatial_det = v[5];
  double lapse = v[6];
  double shift_x = v[7];
  double shift_y = v[8];
  double shift_z = v[9];

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_metric[0][0] = v[10]; spatial_metric[0][1] = v[11]; spatial_metric[0][2] = v[12];
  spatial_metric[1][0] = v[13]; spatial_metric[1][1] = v[14]; spatial_metric[1][2] = v[15];
  spatial_metric[2][0] = v[16]; spatial_metric[2][1] = v[17]; spatial_metric[2][2] = v[18];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  inv_spatial_metric[0][0] = v[19]; inv_spatial_metric[0][1] = v[20]; inv_spatial_metric[0][2] = v[21];
  inv_spatial_metric[1][0] = v[22]; inv_spatial_metric[1][1] = v[23]; inv_spatial_metric[1][2] = v[24];
  inv_spatial_metric[2][0] = v[25]; inv_spatial_metric[2][1] = v[26]; inv_spatial_metric[2][2] = v[27];

  bool in_excision_region = false;
  if (v[28] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }
  
  if (!in_excision_region) {
    double W_flat = 1.0 / (sqrt(1.0 - (vx * vx) - (vy * vy) - (vz * vz)));
    if ((vx * vx) + (vy * vy) + (vz * vz) > 1.0 - pow(1.0, -8.0)) {
      W_flat = 1.0 / sqrt(pow(1.0, -8.0));
    }

    double *vel = gkyl_malloc(sizeof(double[3]));
    double v_sq = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W_curved = 1.0 / sqrt(1.0 - v_sq);
    if(v_sq > 1.0 - pow(10.0, -8.0)) {
      W_curved = 1.0 / sqrt(pow(1.0, -8.0));
    }

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

    if (W_curved > 1.5) {
      printf("%f = %f\n", W_flat, W_curved);
    }
    
    if (fabs(vx) < pow(10.0, -16.0)) {
      flux_extrap[0] = flux[0] + (lapse * sqrt(spatial_det)) * rho * W_curved * (- (shift_x / lapse));
      flux_extrap[1] = (flux[1] - p) + (lapse * sqrt(spatial_det)) * p;
      flux_extrap[2] = flux[2] + (lapse * sqrt(spatial_det)) * (rho * h * (W_curved * W_curved) * (vy * (- (shift_x / lapse))));
      flux_extrap[3] = flux[3] + (lapse * sqrt(spatial_det)) * (rho * h * (W_curved * W_curved) * (vz * (- (shift_x / lapse))));
      flux_extrap[4] = flux[4] + (lapse * sqrt(spatial_det)) * (((rho * h * (W_curved * W_curved)) - p - (rho * W_curved)) * (- (shift_x / lapse)));
    }
    else {
      flux_extrap[0] = (flux[0] / W_flat) * (lapse * sqrt(spatial_det)) * W_curved * (1.0 - (shift_x / (lapse * vx)));
      flux_extrap[1] = ((flux[1] - p) / (W_flat * W_flat)) * (lapse * sqrt(spatial_det)) * (W_curved * W_curved) * (1.0 - (shift_x / (lapse * vx)))
        + (lapse * sqrt(spatial_det)) * p;
      flux_extrap[2] = (flux[2] / (W_flat * W_flat)) * (lapse * sqrt(spatial_det)) * (W_curved * W_curved) * (1.0 - (shift_x / (lapse * vx)));
      flux_extrap[3] = (flux[3] / (W_flat * W_flat)) * (lapse * sqrt(spatial_det)) * (W_curved * W_curved) * (1.0 - (shift_x / (lapse * vx)));
      // The general extrapolation for energy is horrifying, so instead we extract pressure, transform that, and then transform back.
      flux_extrap[4] = (flux[4] / (((rho * h * (W_flat * W_flat)) - (rho * W_flat)) * vx)) * (lapse * sqrt(spatial_det)) *
        (((rho * h * (W_curved * W_curved)) - p - (rho * W_curved)) * (vx - (shift_x / lapse)) + (p * vx));
    }

    for (int i = 5; i < 29; i++) {
      flux_extrap[i] = 0.0;
    }

    gkyl_free(vel);
  }
  else {
    for (int i = 0; i < 29; i++) {
      flux_extrap[i] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(inv_spatial_metric);
}

static double
wave_lax_tetrad(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_euler_tetrad *gr_euler_tetrad = container_of(eqn, struct wv_gr_euler_tetrad, eqn);
  double gas_gamma = gr_euler_tetrad->gas_gamma;

  double sl = gkyl_gr_euler_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_gr_euler_max_abs_speed(gas_gamma, qr);
  double amax = fmax(sl, sr);

  double fl[29], fr[29];
  gkyl_gr_euler_tetrad_flux(gas_gamma, ql, fl);
  gkyl_gr_euler_tetrad_flux(gas_gamma, qr, fr);

  double fl_extrap[29], fr_extrap[29];
  gkyl_gr_euler_tetrad_extrapolate_flux(gas_gamma, ql, fl, fl_extrap);
  gkyl_gr_euler_tetrad_extrapolate_flux(gas_gamma, qr, fr, fr_extrap);

  bool in_excision_region_l = false;
  if (ql[28] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[28] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  double *w0 = &waves[0], *w1 = &waves[29];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 29; i++) {
      w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr_extrap[i] - fl_extrap[i]) / amax);
      w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr_extrap[i] - fl_extrap[i]) / amax);
    }
  }
  else {
    for (int i = 0; i < 29; i++) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }

  s[0] = -amax;
  s[1] = amax;

  return s[1];
}

static void
qfluct_lax_tetrad(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[29];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 29; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_lax_l_tetrad(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  return wave_lax_tetrad(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l_tetrad(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  return qfluct_lax_tetrad(eqn, ql, qr, waves, s, amdq, apdq);
}

void
gkyl_gr_euler_tetrad_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_euler_tetrad *gr_euler_tetrad = container_of(base->on_dev, struct wv_gr_euler_tetrad, eqn);
    gkyl_cu_free(gr_euler_tetrad);
  }

  struct wv_gr_euler_tetrad *gr_euler_tetrad = container_of(base, struct wv_gr_euler_tetrad, eqn);
  gkyl_free(gr_euler_tetrad);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_tetrad_new(double gas_gamma, struct gkyl_gr_spacetime* spacetime, bool use_gpu)
{
  return gkyl_wv_gr_euler_tetrad_inew(&(struct gkyl_wv_gr_euler_tetrad_inp) {
      .gas_gamma = gas_gamma,
      .spacetime = spacetime,
      .rp_type = WV_GR_EULER_TETRAD_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_euler_tetrad_inew(const struct gkyl_wv_gr_euler_tetrad_inp* inp)
{
  struct wv_gr_euler_tetrad *gr_euler_tetrad = gkyl_malloc(sizeof(struct wv_gr_euler_tetrad));

  gr_euler_tetrad->eqn.type = GKYL_EQN_GR_EULER_TETRAD;
  gr_euler_tetrad->eqn.num_equations = 29;
  gr_euler_tetrad->eqn.num_diag = 5;

  gr_euler_tetrad->gas_gamma = inp->gas_gamma;
  gr_euler_tetrad->spacetime = inp->spacetime;

  if (inp->rp_type == WV_GR_EULER_RP_LAX) {
    gr_euler_tetrad->eqn.num_waves = 2;
    gr_euler_tetrad->eqn.waves_func = wave_lax_l_tetrad;
    gr_euler_tetrad->eqn.qfluct_func = qfluct_lax_l_tetrad;
  }

  gr_euler_tetrad->eqn.flux_jump = flux_jump;
  gr_euler_tetrad->eqn.check_inv_func = check_inv;
  gr_euler_tetrad->eqn.max_speed_func = max_speed;
  gr_euler_tetrad->eqn.rotate_to_local_func = rot_to_local;
  gr_euler_tetrad->eqn.rotate_to_global_func = rot_to_global;

  gr_euler_tetrad->eqn.wall_bc_func = gr_euler_wall;
  gr_euler_tetrad->eqn.no_slip_bc_func = gr_euler_no_slip;

  gr_euler_tetrad->eqn.cons_to_riem = cons_to_riem;
  gr_euler_tetrad->eqn.riem_to_cons = riem_to_cons;

  gr_euler_tetrad->eqn.cons_to_diag = gr_euler_cons_to_diag;

  gr_euler_tetrad->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_euler_tetrad->eqn.flags);
  gr_euler_tetrad->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_euler_free);
  gr_euler_tetrad->eqn.on_dev = &gr_euler_tetrad->eqn; // On the CPU, the equation object points to itself.

  return &gr_euler_tetrad->eqn;
}

double
gkyl_wv_gr_euler_tetrad_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_euler_tetrad *gr_euler_tetrad = container_of(eqn, struct wv_gr_euler_tetrad, eqn);
  double gas_gamma = gr_euler_tetrad->gas_gamma;

  return gas_gamma;
}

struct gkyl_gr_spacetime*
gkyl_wv_gr_euler_tetrad_spacetime(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_euler_tetrad *gr_euler_tetrad = container_of(eqn, struct wv_gr_euler_tetrad, eqn);
  struct gkyl_gr_spacetime* spacetime = gr_euler_tetrad->spacetime;

  return spacetime;
}