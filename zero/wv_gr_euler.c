#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h>

static inline void
gkyl_gr_euler_prim_vars(double gas_gamma, const double q[29], double v[29])
{
  double spatial_det = q[5];
  double lapse = q[6];
  double shift_x = q[7];
  double shift_y = q[8];
  double shift_z = q[9];

  double **spatial_metric = malloc(sizeof(double*) * 3);
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = malloc(sizeof(double) * 3);
  }

  spatial_metric[0][0] = q[10]; spatial_metric[0][1] = q[11]; spatial_metric[0][2] = q[12];
  spatial_metric[1][0] = q[13]; spatial_metric[1][1] = q[14]; spatial_metric[1][2] = q[15];
  spatial_metric[2][0] = q[16]; spatial_metric[2][1] = q[17]; spatial_metric[2][2] = q[18];
  
  double **inv_spatial_metric = malloc(sizeof(double*) * 3);
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = malloc(sizeof(double) * 3);
  }

  inv_spatial_metric[0][0] = q[19]; inv_spatial_metric[0][1] = q[20]; inv_spatial_metric[0][2] = q[21];
  inv_spatial_metric[1][0] = q[22]; inv_spatial_metric[1][1] = q[23]; inv_spatial_metric[1][2] = q[24];
  inv_spatial_metric[2][0] = q[25]; inv_spatial_metric[2][1] = q[26]; inv_spatial_metric[2][2] = q[27];

  bool in_excision_region = false;
  if (q[28] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double D = q[0] / sqrt(spatial_det);
    double momx = q[1] / sqrt(spatial_det);
    double momy = q[2] / sqrt(spatial_det);
    double momz = q[3] / sqrt(spatial_det);
    double Etot = q[4] / sqrt(spatial_det);

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

    v[5] = spatial_det;
    v[6] = lapse;
    v[7] = shift_x;
    v[8] = shift_y;
    v[9] = shift_z;

    v[10] = spatial_metric[0][0]; v[11] = spatial_metric[0][1]; v[12] = spatial_metric[0][2];
    v[13] = spatial_metric[1][0]; v[14] = spatial_metric[1][1]; v[15] = spatial_metric[1][2];
    v[16] = spatial_metric[2][0]; v[17] = spatial_metric[2][1]; v[18] = spatial_metric[2][2];

    v[19] = inv_spatial_metric[0][0]; v[20] = inv_spatial_metric[0][1]; v[21] = inv_spatial_metric[0][2];
    v[22] = inv_spatial_metric[1][0]; v[23] = inv_spatial_metric[1][1]; v[24] = inv_spatial_metric[1][2];
    v[25] = inv_spatial_metric[2][0]; v[26] = inv_spatial_metric[2][1]; v[27] = inv_spatial_metric[2][2];

    v[28] = 1.0;
  }
  else {
    for (int i = 0; i < 28; i++) {
      v[i] = 0.0;
    }
    
    v[28] = -1.0;
  }
}

static inline double
gkyl_gr_euler_max_abs_speed(double gas_gamma, const double q[29])
{
  double v[29] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);

  bool in_excision_region = false;
  if (v[28] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double rho = v[0];
    double p = v[4];

    double num = (gas_gamma * p) / rho;
    double den = 1.0 * ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));
    double cs = sqrt(num / den);

    double v_sq = sqrt((v[1] * v[1]) + (v[2] * v[2]) + (v[3] * v[3]));
    return fabs(v_sq) + cs;
  }
  else {
    return pow(10.0, -8.0);
  }
}

static void
gkyl_gr_euler_flux(double gas_gamma, const double q[29], double flux[29])
{
  double v[29] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, q, v);
  double rho =  v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = v[4];

  double spatial_det = v[5];
  double lapse = v[6];
  double shift_x = v[7];
  double shift_y = v[8];
  double shift_z = v[9];

  double **spatial_metric = malloc(sizeof(double*) * 3);
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = malloc(sizeof(double) * 3);
  }

  spatial_metric[0][0] = v[10]; spatial_metric[0][1] = v[11]; spatial_metric[0][2] = v[12];
  spatial_metric[1][0] = v[13]; spatial_metric[1][1] = v[14]; spatial_metric[1][2] = v[15];
  spatial_metric[2][0] = v[16]; spatial_metric[2][1] = v[17]; spatial_metric[2][2] = v[18];

  double **inv_spatial_metric = malloc(sizeof(double*) * 3);
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = malloc(sizeof(double) * 3);
  }

  inv_spatial_metric[0][0] = v[19]; inv_spatial_metric[0][1] = v[20]; inv_spatial_metric[0][2] = v[21];
  inv_spatial_metric[1][0] = v[22]; inv_spatial_metric[1][1] = v[23]; inv_spatial_metric[1][2] = v[24];
  inv_spatial_metric[2][0] = v[25]; inv_spatial_metric[2][1] = v[26]; inv_spatial_metric[2][2] = v[27];

  bool in_excision_region = false;
  if (v[28] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double *vel = malloc(sizeof(double) * 3);
    double v_sq = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W = 1.0 / (sqrt(1.0 - v_sq));
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(1.0 - pow(10.0, -8.0));
    }

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

    flux[0] = (lapse * sqrt(spatial_det)) * (rho * W * (vx - (shift_x / lapse)));
    flux[1] = (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vx * (vx - (shift_x / lapse))) + p);
    flux[2] = (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vy * (vx - (shift_x / lapse))));
    flux[3] = (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vz * (vx - (shift_x / lapse))));
    flux[4] = (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vx - (shift_x / lapse)) + (p * vx));

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

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 29; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 29; i++) {
    qout[i] = win[i];
  }
}

static void
gr_euler_wall(double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 29; i++) {
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

  for (int i = 5; i < 29; i++) {
    ghost[i] = skin[i];
  }
}

static inline void
rot_to_local(const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal, double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);
  qlocal[4] = qglobal[4];

  qlocal[5] = qglobal[5];
  qlocal[6] = qglobal[6];
  qlocal[7] = (qglobal[7] * norm[0]) + (qglobal[8] * norm[1]) + (qglobal[9] * norm[2]);
  qlocal[8] = (qglobal[7] * tau1[0]) + (qglobal[8] * tau1[1]) + (qglobal[9] * tau1[2]);
  qlocal[9] = (qglobal[7] * tau2[0]) + (qglobal[8] * tau2[1]) + (qglobal[9] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qglobal[10] * norm[0]) + (qglobal[11] * norm[1]) + (qglobal[12] * norm[2]);
  r1[1] = (qglobal[10] * tau1[0]) + (qglobal[11] * tau1[1]) + (qglobal[12] * tau1[2]);
  r1[2] = (qglobal[10] * tau2[0]) + (qglobal[11] * tau2[1]) + (qglobal[12] * tau2[2]);

  r2[0] = (qglobal[13] * norm[0]) + (qglobal[14] * norm[1]) + (qglobal[15] * norm[2]);
  r2[1] = (qglobal[13] * tau1[0]) + (qglobal[14] * tau1[1]) + (qglobal[15] * tau1[2]);
  r2[2] = (qglobal[13] * tau2[0]) + (qglobal[14] * tau2[1]) + (qglobal[15] * tau2[2]);

  r3[0] = (qglobal[16] * norm[0]) + (qglobal[17] * norm[1]) + (qglobal[18] * norm[2]);
  r3[1] = (qglobal[16] * tau1[0]) + (qglobal[17] * tau1[1]) + (qglobal[18] * tau1[2]);
  r3[2] = (qglobal[16] * tau2[0]) + (qglobal[17] * tau2[1]) + (qglobal[18] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double v1[3], v2[3], v3[3];
  v1[0] = (r1[0] * norm[0]) + (r2[0] * norm[1]) + (r3[0] * norm[2]);
  v1[1] = (r1[0] * tau1[0]) + (r2[0] * tau1[1]) + (r3[0] * tau1[2]);
  v1[2] = (r1[0] * tau2[0]) + (r2[0] * tau2[1]) + (r3[0] * tau2[2]);

  v2[0] = (r1[1] * norm[0]) + (r2[1] * norm[1]) + (r3[1] * norm[2]);
  v2[1] = (r1[1] * tau1[0]) + (r2[1] * tau1[1]) + (r3[1] * tau1[2]);
  v2[2] = (r1[1] * tau2[0]) + (r2[1] * tau2[1]) + (r3[1] * tau2[2]);

  v3[0] = (r1[2] * norm[0]) + (r2[2] * norm[1]) + (r3[2] * norm[2]);
  v3[1] = (r1[2] * tau1[0]) + (r2[2] * tau1[1]) + (r3[2] * tau1[2]);
  v3[2] = (r1[2] * tau2[0]) + (r2[2] * tau2[1]) + (r3[2] * tau2[2]);

  // Rotate spatial metric tensor to global coordinate frame.
  qlocal[10] = v1[0]; qlocal[11] = v1[1]; qlocal[12] = v1[2];
  qlocal[13] = v2[0]; qlocal[14] = v2[1]; qlocal[15] = v2[2];
  qlocal[16] = v3[0]; qlocal[17] = v3[1]; qlocal[18] = v3[2];

  // Temporary arrays to store rotated inverse column vectors.
  double inv_r1[3], inv_r2[3], inv_r3[3];
  inv_r1[0] = (qglobal[19] * norm[0]) + (qglobal[20] * norm[1]) + (qglobal[21] * norm[2]);
  inv_r1[1] = (qglobal[19] * tau1[0]) + (qglobal[20] * tau1[1]) + (qglobal[21] * tau1[2]);
  inv_r1[2] = (qglobal[19] * tau2[0]) + (qglobal[20] * tau2[1]) + (qglobal[21] * tau2[2]);

  inv_r2[0] = (qglobal[22] * norm[0]) + (qglobal[23] * norm[1]) + (qglobal[24] * norm[2]);
  inv_r2[1] = (qglobal[22] * tau1[0]) + (qglobal[23] * tau1[1]) + (qglobal[24] * tau1[2]);
  inv_r2[2] = (qglobal[22] * tau2[0]) + (qglobal[23] * tau2[1]) + (qglobal[24] * tau2[2]);

  inv_r3[0] = (qglobal[25] * norm[0]) + (qglobal[26] * norm[1]) + (qglobal[27] * norm[2]);
  inv_r3[1] = (qglobal[25] * tau1[0]) + (qglobal[26] * tau1[1]) + (qglobal[27] * tau1[2]);
  inv_r3[2] = (qglobal[25] * tau2[0]) + (qglobal[26] * tau2[1]) + (qglobal[27] * tau2[2]);

  // Temporary arrays to store rotated inverse row vectors.
  double inv_v1[3], inv_v2[3], inv_v3[3];
  inv_v1[0] = (inv_r1[0] * norm[0]) + (inv_r2[0] * norm[1]) + (inv_r3[0] * norm[2]);
  inv_v1[1] = (inv_r1[0] * tau1[0]) + (inv_r2[0] * tau1[1]) + (inv_r3[0] * tau1[2]);
  inv_v1[2] = (inv_r1[0] * tau2[0]) + (inv_r2[0] * tau2[1]) + (inv_r3[0] * tau2[2]);

  inv_v2[0] = (inv_r1[1] * norm[0]) + (inv_r2[1] * norm[1]) + (inv_r3[1] * norm[2]);
  inv_v2[1] = (inv_r1[1] * tau1[0]) + (inv_r2[1] * tau1[1]) + (inv_r3[1] * tau1[2]);
  inv_v2[2] = (inv_r1[1] * tau2[0]) + (inv_r2[1] * tau2[1]) + (inv_r3[1] * tau2[2]);

  inv_v3[0] = (inv_r1[2] * norm[0]) + (inv_r2[2] * norm[1]) + (inv_r3[2] * norm[2]);
  inv_v3[1] = (inv_r1[2] * tau1[0]) + (inv_r2[2] * tau1[1]) + (inv_r3[2] * tau1[2]);
  inv_v3[2] = (inv_r1[2] * tau2[0]) + (inv_r2[2] * tau2[1]) + (inv_r3[2] * tau2[2]);

  // Rotate inverse spatial metric tensor to global coordinate frame.
  qlocal[19] = inv_v1[0]; qlocal[20] = inv_v1[1]; qlocal[21] = inv_v1[2];
  qlocal[22] = inv_v2[0]; qlocal[23] = inv_v2[1]; qlocal[24] = inv_v2[2];
  qlocal[25] = inv_v3[0]; qlocal[26] = inv_v3[1]; qlocal[27] = inv_v3[2];

  qlocal[28] = qglobal[28];
}

static inline void
rot_to_global(const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal, double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);
  qglobal[4] = qlocal[4];

  qglobal[5] = qlocal[5];
  qglobal[6] = qlocal[6];
  qglobal[7] = (qlocal[7] * norm[0]) + (qlocal[8] * tau1[0]) + (qlocal[9] * tau2[0]);
  qglobal[8] = (qlocal[7] * norm[1]) + (qlocal[8] * tau1[1]) + (qlocal[9] * tau2[1]);
  qglobal[9] = (qlocal[7] * norm[2]) + (qlocal[8] * tau1[2]) + (qlocal[9] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qlocal[10] * norm[0]) + (qlocal[11] * tau1[0]) + (qlocal[12] * tau2[0]);
  r1[1] = (qlocal[10] * norm[1]) + (qlocal[11] * tau1[1]) + (qlocal[12] * tau2[1]);
  r1[2] = (qlocal[10] * norm[2]) + (qlocal[11] * tau1[2]) + (qlocal[12] * tau2[2]);

  r2[0] = (qlocal[13] * norm[0]) + (qlocal[14] * tau1[0]) + (qlocal[15] * tau2[0]);
  r2[1] = (qlocal[13] * norm[1]) + (qlocal[14] * tau1[1]) + (qlocal[15] * tau2[1]);
  r2[2] = (qlocal[13] * norm[2]) + (qlocal[14] * tau1[2]) + (qlocal[15] * tau2[2]);

  r3[0] = (qlocal[16] * norm[0]) + (qlocal[17] * tau1[0]) + (qlocal[18] * tau2[0]);
  r3[1] = (qlocal[16] * norm[1]) + (qlocal[17] * tau1[1]) + (qlocal[18] * tau2[1]);
  r3[2] = (qlocal[16] * norm[2]) + (qlocal[17] * tau1[2]) + (qlocal[18] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double v1[3], v2[3], v3[3];
  v1[0] = (r1[0] * norm[0]) + (r2[0] * tau1[0]) + (r3[0] * tau2[0]);
  v1[1] = (r1[0] * norm[1]) + (r2[0] * tau1[1]) + (r3[0] * tau2[1]);
  v1[2] = (r1[0] * norm[2]) + (r2[0] * tau1[2]) + (r3[0] * tau2[2]);

  v2[0] = (r1[1] * norm[0]) + (r2[1] * tau1[0]) + (r3[1] * tau2[0]);
  v2[1] = (r1[1] * norm[1]) + (r2[1] * tau1[1]) + (r3[1] * tau2[1]);
  v2[2] = (r1[1] * norm[2]) + (r2[1] * tau1[2]) + (r3[1] * tau2[2]);

  v3[0] = (r1[2] * norm[0]) + (r2[2] * tau1[0]) + (r3[2] * tau2[0]);
  v3[1] = (r1[2] * norm[1]) + (r2[2] * tau1[1]) + (r3[2] * tau2[1]);
  v3[2] = (r1[2] * norm[2]) + (r2[2] * tau1[2]) + (r3[2] * tau2[2]);

  // Rotate spatial metric tensor back to local coordinate frame.
  qglobal[10] = v1[0]; qglobal[11] = v1[1]; qglobal[12] = v1[2];
  qglobal[13] = v2[0]; qglobal[14] = v2[1]; qglobal[15] = v2[2];
  qglobal[16] = v3[0]; qglobal[17] = v3[1]; qglobal[18] = v3[2];

  // Temporary arrays to store rotated inverse column vectors.
  double inv_r1[3], inv_r2[3], inv_r3[3];
  inv_r1[0] = (qlocal[19] * norm[0]) + (qlocal[20] * tau1[0]) + (qlocal[21] * tau2[0]);
  inv_r1[1] = (qlocal[19] * norm[1]) + (qlocal[20] * tau1[1]) + (qlocal[21] * tau2[1]);
  inv_r1[2] = (qlocal[19] * norm[2]) + (qlocal[20] * tau1[2]) + (qlocal[21] * tau2[2]);

  inv_r2[0] = (qlocal[22] * norm[0]) + (qlocal[23] * tau1[0]) + (qlocal[24] * tau2[0]);
  inv_r2[1] = (qlocal[22] * norm[1]) + (qlocal[23] * tau1[1]) + (qlocal[24] * tau2[1]);
  inv_r2[2] = (qlocal[22] * norm[2]) + (qlocal[23] * tau1[2]) + (qlocal[24] * tau2[2]);

  inv_r3[0] = (qlocal[25] * norm[0]) + (qlocal[26] * tau1[0]) + (qlocal[27] * tau2[0]);
  inv_r3[1] = (qlocal[25] * norm[1]) + (qlocal[26] * tau1[1]) + (qlocal[27] * tau2[1]);
  inv_r3[2] = (qlocal[25] * norm[2]) + (qlocal[26] * tau1[2]) + (qlocal[27] * tau2[2]);

  // Temporary arrays to store rotated inverse row vectors.
  double inv_v1[3], inv_v2[3], inv_v3[3];
  inv_v1[0] = (inv_r1[0] * norm[0]) + (inv_r2[0] * tau1[0]) + (inv_r3[0] * tau2[0]);
  inv_v1[1] = (inv_r1[0] * norm[1]) + (inv_r2[0] * tau1[1]) + (inv_r3[0] * tau2[1]);
  inv_v1[2] = (inv_r1[0] * norm[2]) + (inv_r2[0] * tau1[2]) + (inv_r3[0] * tau2[2]);

  inv_v2[0] = (inv_r1[1] * norm[0]) + (inv_r2[1] * tau1[0]) + (inv_r3[1] * tau2[0]);
  inv_v2[1] = (inv_r1[1] * norm[1]) + (inv_r2[1] * tau1[1]) + (inv_r3[1] * tau2[1]);
  inv_v2[2] = (inv_r1[1] * norm[2]) + (inv_r2[1] * tau1[2]) + (inv_r3[1] * tau2[2]);

  inv_v3[0] = (inv_r1[2] * norm[0]) + (inv_r2[2] * tau1[0]) + (inv_r3[2] * tau2[0]);
  inv_v3[1] = (inv_r1[2] * norm[1]) + (inv_r2[2] * tau1[1]) + (inv_r3[2] * tau2[1]);
  inv_v3[2] = (inv_r1[2] * norm[2]) + (inv_r2[2] * tau1[2]) + (inv_r3[2] * tau2[2]);

  // Rotate inverse spatial metric tensor back to local coordinate frame.
  qglobal[19] = inv_v1[0]; qglobal[20] = inv_v1[1]; qglobal[21] = inv_v1[2];
  qglobal[22] = inv_v2[0]; qglobal[23] = inv_v2[1]; qglobal[24] = inv_v2[2];
  qglobal[25] = inv_v3[0]; qglobal[26] = inv_v3[1]; qglobal[27] = inv_v3[2];

  qglobal[28] = qlocal[28];
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  double sl = gkyl_gr_euler_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_gr_euler_max_abs_speed(gas_gamma, qr);
  double amax = fmax(sl, sr);

  double fl[29], fr[29];
  gkyl_gr_euler_flux(gas_gamma, ql, fl);
  gkyl_gr_euler_flux(gas_gamma, qr, fr);

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
      w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
      w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
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
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
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
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);

  double fr[29], fl[29];
  gkyl_gr_euler_flux(gr_euler->gas_gamma, ql, fl);
  gkyl_gr_euler_flux(gr_euler->gas_gamma, qr, fr);

  bool in_excision_region_l = false;
  if (ql[28] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[28] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  if (!in_excision_region_l && !in_excision_region_r) {
    for (int m = 0; m < 29; m++) {
      flux_jump[m] = fr[m] - fl[m];
    }
  }
  else {
    for (int m = 0; m < 29; m++) {
      flux_jump[m] = 0.0;
    }
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

  double v[29] = { 0.0 };
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
gkyl_wv_gr_euler_new(double gas_gamma, struct gkyl_gr_spacetime* spacetime, bool use_gpu)
{
  return gkyl_wv_gr_euler_inew(&(struct gkyl_wv_gr_euler_inp) {
      .gas_gamma = gas_gamma,
      .spacetime = spacetime,
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
  gr_euler->eqn.num_equations = 29;
  gr_euler->eqn.num_diag = 5;

  gr_euler->gas_gamma = inp->gas_gamma;
  gr_euler->spacetime = inp->spacetime;

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