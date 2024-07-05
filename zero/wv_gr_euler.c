#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h>

void
gkyl_gr_euler_prim_vars(double gas_gamma, const double q[29], double v[29])
{
  double spatial_det = q[5];
  double lapse = q[6];
  double shift_x = q[7];
  double shift_y = q[8];
  double shift_z = q[9];

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_metric[0][0] = q[10]; spatial_metric[0][1] = q[11]; spatial_metric[0][2] = q[12];
  spatial_metric[1][0] = q[13]; spatial_metric[1][1] = q[14]; spatial_metric[1][2] = q[15];
  spatial_metric[2][0] = q[16]; spatial_metric[2][1] = q[17]; spatial_metric[2][2] = q[18];
  
  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
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

    while (iter < 100) {
      double poly = (alpha4 * (guess * guess * guess) * (guess - eta)) + (alpha2 * (guess * guess)) + (alpha1 * guess) + alpha0;
      double poly_der = alpha1 + (2.0 * alpha2 * guess) + (4.0 * alpha4 * (guess * guess * guess)) - (3.0 * eta * alpha4 * (guess * guess));

      double guess_new = guess - (poly / poly_der);

      if (fabs(guess - guess_new) < pow(10.0, -8.0)) {
        iter = 100;
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

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(inv_spatial_metric);
}

double
gkyl_gr_euler_max_abs_speed(double gas_gamma, const double q[29])
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

  double num = (gas_gamma * p) / rho;
  double den = 1.0 + ((p / rho) * (gas_gamma) / (gas_gamma - 1.0));
  double c_s = sqrt(num / den);

  bool in_excision_region = false;
  if (v[28] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  bool curved_spacetime = false;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        if (fabs(spatial_metric[i][j] - 1.0) > pow(10.0, -8.0)) {
          curved_spacetime = true;
        }
      }
      else {
        if (fabs(spatial_metric[i][j]) > pow(10.0, -8.0)) {
          curved_spacetime = true;
        }
      }
    }
  }
  if (fabs(spatial_det - 1.0) > pow(10.0, -8.0) || fabs(lapse - 1.0) > pow(10.0, -8.0) || fabs(shift_x) > pow(10.0, -8.0) ||
    fabs(shift_y) > pow(10.0, -8.0) || fabs(shift_z) > pow(10.0, -8.0)) {
    curved_spacetime = true;
  }

  if (!in_excision_region) {
    if (curved_spacetime) {
      double *vel = gkyl_malloc(sizeof(double[3]));
      double v_sq = 0.0;
      vel[0] = vx; vel[1] = vy; vel[2] = vz;

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double *shift = gkyl_malloc(sizeof(double[3]));
      shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

      double *material_eigs = gkyl_malloc(sizeof(double[3]));
      double *fast_acoustic_eigs = gkyl_malloc(sizeof(double[3]));
      double *slow_acoustic_eigs = gkyl_malloc(sizeof(double[3]));

      for (int i = 0; i < 3; i++) {
        material_eigs[i] = (lapse * vel[i]) - shift[i];

        fast_acoustic_eigs[i] = (lapse / (1.0 - (v_sq * (c_s * c_s)))) * ((vel[i] * (1.0 - (c_s * c_s))) +
          (c_s * sqrt((1.0 - v_sq) * (inv_spatial_metric[i][i] * (1.0 - (v_sq * (c_s * c_s))) - (vel[i] * vel[i]) * (1.0 - (c_s * c_s)))))) - shift[i];
        
        slow_acoustic_eigs[i] = (lapse / (1.0 - (v_sq * (c_s * c_s)))) * ((vel[i] * (1.0 - (c_s * c_s))) -
          (c_s * sqrt((1.0 - v_sq) * (inv_spatial_metric[i][i] * (1.0 - (v_sq * (c_s * c_s))) - (vel[i] * vel[i]) * (1.0 - (c_s * c_s)))))) - shift[i];
      }

      double max_eig = 0.0;
      for (int i = 0; i < 3; i++) {
        if (fabs(material_eigs[i]) > max_eig) {
          max_eig = fabs(material_eigs[i]);
        }
        if (fabs(fast_acoustic_eigs[i]) > max_eig) {
          max_eig = fabs(fast_acoustic_eigs[i]);
        }
        if (fabs(slow_acoustic_eigs[i]) > max_eig) {
          max_eig = fabs(slow_acoustic_eigs[i]);
        }
      }

      gkyl_free(vel);
      gkyl_free(shift);
      gkyl_free(material_eigs);
      gkyl_free(fast_acoustic_eigs);
      gkyl_free(slow_acoustic_eigs);

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(inv_spatial_metric[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(inv_spatial_metric);

      return fabs(v_sq) + max_eig;
    }
    else {
      double v_sq = sqrt((vx * vx) + (vy * vy) + (vz * vz));

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(inv_spatial_metric[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(inv_spatial_metric);

      return fabs(v_sq) + c_s;
    }
  }
  else {
    for (int i = 0; i < 3; i++) {
      gkyl_free(spatial_metric[i]);
      gkyl_free(inv_spatial_metric[i]);
    }
    gkyl_free(spatial_metric);
    gkyl_free(inv_spatial_metric);

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
    double *vel = gkyl_malloc(sizeof(double[3]));
    double v_sq = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W = 1.0 / (sqrt(1.0 - v_sq));
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(pow(10.0, -8.0));
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

    gkyl_free(vel);
  }
  else {
    for (int i = 0; i < 29; i++) {
      flux[i] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(inv_spatial_metric);
}

void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 29; i++) {
    wout[i] = qin[i];
  }
}

void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 29; i++) {
    qout[i] = win[i];
  }
}

void
gr_euler_wall(double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 29; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

void
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

void
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

void
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
wave_roe(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  double* delta_mod = gkyl_malloc(sizeof(double[5]));
  delta_mod[0] = delta[0];
  delta_mod[1] = delta[4];
  delta_mod[2] = delta[1];
  delta_mod[3] = delta[2];
  delta_mod[4] = delta[3];

  double vl[29] = { 0.0 };
  double vr[29] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, ql, vl);
  gkyl_gr_euler_prim_vars(gas_gamma, qr, vr);

  double rho_l = vl[0];
  double vx_l = vl[1];
  double vy_l = vl[2];
  double vz_l = vl[3];
  double p_l = vl[4];

  double spatial_det_l = vl[5];
  double lapse_l = vl[6];
  double shift_x_l = vl[7];
  double shift_y_l = vl[8];
  double shift_z_l = vl[9];

  double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_metric_l[0][0] = vl[10]; spatial_metric_l[0][1] = vl[11]; spatial_metric_l[0][2] = vl[12];
  spatial_metric_l[1][0] = vl[13]; spatial_metric_l[1][1] = vl[14]; spatial_metric_l[1][2] = vl[15];
  spatial_metric_l[2][0] = vl[16]; spatial_metric_l[2][1] = vl[17]; spatial_metric_l[2][2] = vl[18];

  double **inv_spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
  }

  inv_spatial_metric_l[0][0] = vl[19]; inv_spatial_metric_l[0][1] = vl[20]; inv_spatial_metric_l[0][2] = vl[21];
  inv_spatial_metric_l[1][0] = vl[22]; inv_spatial_metric_l[1][1] = vl[23]; inv_spatial_metric_l[1][2] = vl[24];
  inv_spatial_metric_l[2][0] = vl[25]; inv_spatial_metric_l[2][1] = vl[26]; inv_spatial_metric_l[2][2] = vl[27];

  bool in_excision_region_l = false;
  if (ql[28] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  double rho_r = vr[0];
  double vx_r = vr[1];
  double vy_r = vr[2];
  double vz_r = vr[3];
  double p_r = vr[4];

  double spatial_det_r = vr[5];
  double lapse_r = vr[6];
  double shift_x_r = vr[7];
  double shift_y_r = vr[8];
  double shift_z_r = vr[9];

  double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_metric_r[0][0] = vr[10]; spatial_metric_r[0][1] = vr[11]; spatial_metric_r[0][2] = vr[12];
  spatial_metric_r[1][0] = vr[13]; spatial_metric_r[1][1] = vr[14]; spatial_metric_r[1][2] = vr[15];
  spatial_metric_r[2][0] = vr[16]; spatial_metric_r[2][1] = vr[17]; spatial_metric_r[2][2] = vr[18];

  double **inv_spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
  }

  inv_spatial_metric_r[0][0] = vr[19]; inv_spatial_metric_r[0][1] = vr[20]; inv_spatial_metric_r[0][2] = vr[21];
  inv_spatial_metric_r[1][0] = vr[22]; inv_spatial_metric_r[1][1] = vr[23]; inv_spatial_metric_r[1][2] = vr[24];
  inv_spatial_metric_r[2][0] = vr[25]; inv_spatial_metric_r[2][1] = vr[26]; inv_spatial_metric_r[2][2] = vr[27];

  bool in_excision_region_r = false;
  if (qr[28] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  bool curved_spacetime = false;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        if (fabs(spatial_metric_l[i][j] - 1.0) > pow(10.0, -8.0) || fabs(spatial_metric_r[i][j] - 1.0) > pow(10.0, -8.0)) {
          curved_spacetime = true;
        }
      }
      else {
        if (fabs(spatial_metric_l[i][j]) > pow(10.0, -8.0) || fabs(spatial_metric_r[i][j]) > pow(10.0, -8.0)) {
          curved_spacetime = true;
        }
      }
    }
  }
  if (fabs(spatial_det_l - 1.0) > pow(10.0, -8.0) || fabs(spatial_det_r - 1.0) > pow(10.0, -8.0) || fabs(lapse_l - 1.0) > pow(10.0, -8.0) ||
    fabs(lapse_r - 1.0) > pow(10.0, -8.0) || fabs(shift_x_l) > pow(10.0, -8.0) || fabs(shift_x_r) > pow(10.0, -8.0) ||
    fabs(shift_y_l) > pow(10.0, -8.0) || fabs(shift_y_r) > pow(10.0, -8.0) || fabs(shift_z_l) > pow(10.0, -8.0) ||
    fabs(shift_z_r) > pow(10.0, -8.0)) {
    curved_spacetime = true;
  }

  if (curved_spacetime) {
    double *vel_l = gkyl_malloc(sizeof(double[3]));
    double v_sq_l = 0.0;
    vel_l[0] = vx_l; vel_l[1] = vy_l; vel_l[2] = vz_l;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_l += spatial_metric_l[i][j] * vel_l[i] * vel_l[j];
      }
    }

    double W_l = 1.0 / sqrt(1.0 - v_sq_l);
    if (v_sq_l > 1.0 - pow(10.0, -8.0)) {
      W_l = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double *vel_r = gkyl_malloc(sizeof(double[3]));
    double v_sq_r = 0.0;
    vel_r[0] = vx_r; vel_r[1] = vy_r; vel_r[2] = vz_r;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_r += spatial_metric_r[i][j] * vel_r[i] * vel_r[j];
      }
    }

    double W_r = 1.0 / sqrt(1.0 - v_sq_r);
    if (v_sq_r > 1.0 - pow(10.0, -8.0)) {
      W_r = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double *shift_l = gkyl_malloc(sizeof(double[3]));
    double *coshift_l = gkyl_malloc(sizeof(double[3]));
    shift_l[0] = shift_x_l; shift_l[1] = shift_y_l; shift_l[2] = shift_z_l;
    coshift_l[0] = 0.0; coshift_l[1] = 0.0; coshift_l[2] = 0.0;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        coshift_l[i] += spatial_metric_l[i][j] * shift_l[j];
      }
    }

    double **spacetime_metric_l = gkyl_malloc(sizeof(double*[4]));
    for (int i = 0; i < 4; i++) {
      spacetime_metric_l[i] = gkyl_malloc(sizeof(double[4]));

      for (int j = 0; j < 4; j++) {
        spacetime_metric_l[i][j] = 0.0;
      }
    }

    for (int i = 0; i < 3; i++) {
      spacetime_metric_l[0][0] += shift_l[i] * coshift_l[i];
    }
    spacetime_metric_l[0][0] -= (lapse_l * lapse_l);

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        spacetime_metric_l[0][i + 1] += spatial_metric_l[i][j] * shift_l[j];
        spacetime_metric_l[i + 1][0] += spatial_metric_l[i][j] * shift_l[j];

        spacetime_metric_l[i + 1][j + 1] = spatial_metric_l[i][j];
      }
    }

    double **inv_spacetime_metric_l = gkyl_malloc(sizeof(double*[4]));
    for (int i = 0; i < 4; i++) {
      inv_spacetime_metric_l[i] = gkyl_malloc(sizeof(double[4]));

      for (int j = 0; j < 4; j++) {
        inv_spacetime_metric_l[i][j] = 0.0;
      }
    }

    inv_spacetime_metric_l[0][0] = -1.0 / (lapse_l * lapse_l);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spacetime_metric_l[0][i + 1] += (spatial_metric_l[i][j] * shift_l[j]) / (lapse_l * lapse_l);
        inv_spacetime_metric_l[i + 1][0] += (spatial_metric_l[i][j] * shift_l[j]) / (lapse_l * lapse_l);

        inv_spacetime_metric_l[i + 1][j + 1] = inv_spatial_metric_l[i][j] - (shift_l[i] * shift_l[j]) / (lapse_l * lapse_l);
      }
    }

    double *shift_r = gkyl_malloc(sizeof(double[3]));
    double *coshift_r = gkyl_malloc(sizeof(double[3]));
    shift_r[0] = shift_x_r; shift_r[1] = shift_y_r; shift_r[2] = shift_z_r;
    coshift_r[0] = 0.0; coshift_r[1] = 0.0; coshift_r[2] = 0.0;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        coshift_r[i] += spatial_metric_r[i][j] * shift_r[j];
      }
    }

    double **spacetime_metric_r = gkyl_malloc(sizeof(double*[4]));
    for (int i = 0; i < 4; i++) {
      spacetime_metric_r[i] = gkyl_malloc(sizeof(double[4]));

      for (int j = 0; j < 4; j++) {
        spacetime_metric_r[i][j] = 0.0;
      }
    }

    for (int i = 0; i < 3; i++) {
      spacetime_metric_r[0][0] += shift_r[i] * coshift_r[i];
    }
    spacetime_metric_r[0][0] -= (lapse_r * lapse_r);

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        spacetime_metric_r[0][i + 1] += spatial_metric_r[i][j] * shift_r[j];
        spacetime_metric_r[i + 1][0] += spatial_metric_r[i][j] * shift_r[j];

        spacetime_metric_r[i + 1][j + 1] = spatial_metric_r[i][j];
      }
    }

    double **inv_spacetime_metric_r = gkyl_malloc(sizeof(double*[4]));
    for (int i = 0; i < 4; i++) {
      inv_spacetime_metric_r[i] = gkyl_malloc(sizeof(double[4]));

      for (int j = 0; j < 4; j++) {
        inv_spacetime_metric_r[i][j] = 0.0;
      }
    }

    inv_spacetime_metric_r[0][0] = -1.0 / (lapse_r * lapse_r);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spacetime_metric_r[0][i + 1] += (spatial_metric_r[i][j] * shift_r[j]) / (lapse_r * lapse_r);
        inv_spacetime_metric_r[i + 1][0] += (spatial_metric_r[i][j] * shift_r[j]) / (lapse_r * lapse_r);

        inv_spacetime_metric_r[i + 1][j + 1] = inv_spatial_metric_r[i][j] - (shift_r[i] * shift_r[j]) / (lapse_r * lapse_r);
      }
    }

    double h_l = 1.0 + ((p_l / rho_l) * (gas_gamma / (gas_gamma - 1.0)));
    double h_r = 1.0 + ((p_r / rho_r) * (gas_gamma / (gas_gamma - 1.0)));

    double K_l = lapse_l * sqrt(spatial_det_l) * rho_l * h_l;
    double K_r = lapse_r * sqrt(spatial_det_r) * rho_r * h_r;
    //double K_l = sqrt(ql[4] + p_l) / W_l;
    //double K_r = sqrt(qr[4] + p_r) / W_r;
    double inv_K_avg = 1.0 / (K_l + K_r);

    double v0 = ((K_l * sqrt(spatial_det_l) * W_l) + (K_r * sqrt(spatial_det_r) * W_r)) * inv_K_avg;
    double v1 = ((K_l * sqrt(spatial_det_l) * W_l * vx_l) + (K_r * sqrt(spatial_det_r) * W_r * vx_r)) * inv_K_avg;
    double v2 = ((K_l * sqrt(spatial_det_l) * W_l * vy_l) + (K_r * sqrt(spatial_det_r) * W_r * vy_r)) * inv_K_avg;
    double v3 = ((K_l * sqrt(spatial_det_l) * W_l * vz_l) + (K_r * sqrt(spatial_det_r) * W_r * vz_r)) * inv_K_avg;
    double v4 = (((K_l * p_l) / (rho_l * h_l)) + ((K_r * p_r) / (rho_r * h_r))) * inv_K_avg;
    //double v4 = (p_l / K_l + p_r / K_r) * inv_K_avg;
    double c_minus = 1.0 - (gas_gamma / (gas_gamma - 1.0)) * v4;
    double c_plus = 1.0 + (gas_gamma / (gas_gamma - 1.0)) * v4;

    double **spacetime_metric_avg = gkyl_malloc(sizeof(double*[4]));
    for (int i = 0; i < 4; i++) {
      spacetime_metric_avg[i] = gkyl_malloc(sizeof(double[4]));
    }

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        spacetime_metric_avg[i][j] = (0.5 * spacetime_metric_l[i][j]) + (0.5 * spacetime_metric_r[i][j]);
      }
    }

    double *v_vect = gkyl_malloc(sizeof(double[4]));
    double *v_covect = gkyl_malloc(sizeof(double[4]));
    v_vect[0] = v0; v_vect[1] = v1; v_vect[2] = v2; v_vect[3] = v3;
    v_covect[0] = 0.0; v_covect[1] = 0.0; v_covect[2] = 0.0; v_covect[3] = 0.0;

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        v_covect[i] += spacetime_metric_l[i][j] * v_vect[j];
      }
    }

    double v_contr = 0.0;
    for (int i = 0; i < 4; i++) {
      v_contr += v_vect[i] * v_covect[i];
    }

    double **inv_spacetime_metric_avg = gkyl_malloc(sizeof(double*[4]));
    for (int i = 0; i < 4; i++) {
      inv_spacetime_metric_avg[i] = gkyl_malloc(sizeof(double[4]));
    }

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        inv_spacetime_metric_avg[i][j] = (0.5 * inv_spacetime_metric_l[i][j]) + (0.5 * inv_spacetime_metric_r[i][j]);
      }
    }

    double s_sq = (0.5 * gas_gamma * v4 * (1.0 - v_contr)) - (0.5 * (gas_gamma - 1.0) * (1.0 + v_contr));
    double e = (inv_spacetime_metric_avg[0][0] * (v1 * v1)) - (2.0 * inv_spacetime_metric_avg[0][1] * (v0 * v1)) +
      (inv_spacetime_metric_avg[1][1] * (v0 * v0));
    double y = sqrt(((1.0 - gas_gamma * v4) * e) + (s_sq * ((inv_spacetime_metric_avg[0][1] * inv_spacetime_metric_avg[0][1]) -
      (inv_spacetime_metric_avg[0][0] - inv_spacetime_metric_avg[1][1]))));
    
    double k = (inv_spacetime_metric_avg[0][0] * v1 * delta_mod[2]) - (inv_spacetime_metric_avg[0][1] * ((v0 * delta_mod[2]) +
      (v1 * delta_mod[1]))) + (inv_spacetime_metric_avg[1][1] * v0 * delta_mod[1]);
    double v_delta_contr = 0.0;
    for (int i = 0; i < 4; i++) {
      v_delta_contr += v_covect[i] * delta_mod[i + 1];
    }

    double a1 = (-1.0 / (2.0 * e * s_sq)) * ((s_sq * k) + (sqrt(s_sq) * y * ((v0 * delta_mod[2]) - (v1 * delta_mod[1]))) +
      ((gas_gamma - 1.0) * e * (delta_mod[0] + (c_plus * v_delta_contr))));
    double a2 = (-1.0 / (2.0 * e * s_sq)) * ((s_sq * k) - (sqrt(s_sq) * y * ((v0 * delta_mod[2]) - (v1 * delta_mod[1]))) +
      ((gas_gamma - 1.0) * e * (delta_mod[0] + (c_plus * v_delta_contr))));
    double a3 = (1.0 / (e * s_sq)) * ((2.0 * s_sq * k) + ((gas_gamma - 1.0) * e * (delta_mod[0] + (c_plus * v_delta_contr))));
    double a4 = delta_mod[3] + ((1.0 / e) * (((inv_spacetime_metric_avg[0][2] * v1) - (inv_spacetime_metric_avg[1][2] * v0)) *
      ((v0 * delta_mod[2]) - (v1 * delta_mod[1])) - (k * v2)));
    double a5 = delta_mod[4] + ((1.0 / e) * (((inv_spacetime_metric_avg[0][3] * v1) - (inv_spacetime_metric_avg[1][3] * v0)) *
      ((v0 * delta_mod[2]) - (v1 * delta_mod[1])) - (k * v3)));

    double *wv;
    wv = &waves[0];

    if (!in_excision_region_l && !in_excision_region_r) {
      wv[0] = a1 * c_minus;
      wv[4] = a1 * (v0 - (sqrt(s_sq) / y) * ((inv_spacetime_metric_avg[1][0] * v0) - (inv_spacetime_metric_avg[0][0] * v1)));
      wv[1] = a1 * (v1 - (sqrt(s_sq) / y) * ((inv_spacetime_metric_avg[1][1] * v0) - (inv_spacetime_metric_avg[0][1] * v1)));
      wv[2] = a1 * (v2 - (sqrt(s_sq) / y) * ((inv_spacetime_metric_avg[1][2] * v0) - (inv_spacetime_metric_avg[0][2] * v1)));
      wv[3] = a1 * (v3 - (sqrt(s_sq) / y) * ((inv_spacetime_metric_avg[1][3] * v0) - (inv_spacetime_metric_avg[0][3] * v1)));
      for (int i = 5; i < 29; i++) {
        wv[i] = 0.0;
      }
      s[0] = (((1.0 - (gas_gamma * v4)) * (v0 * v1)) - (s_sq * inv_spacetime_metric_avg[0][1]) - (sqrt(s_sq) * y)) /
        (((1.0 - (gas_gamma * v4)) * (v0 * v0)) - (s_sq * inv_spacetime_metric_avg[0][0]));
    }
    else {
      for (int i = 0; i < 29; i++) {
        wv[i] = 0.0;
      }
      s[0] = pow(10.0, -8.0);
    }

    wv = &waves[29];
    if (!in_excision_region_l && !in_excision_region_r) {
      wv[0] = (a3 * (c_minus + (s_sq / (gas_gamma - 1.0)))) - (a4 * (c_plus * v2)) - (a5 * (c_plus * v3));
      wv[4] = a3 * v0;
      wv[1] = a3 * v1;
      wv[2] = (a3 * v2) + a4;
      wv[3] = (a3 * v3) + a5;
      for (int i = 5; i < 29; i++) {
        wv[i] = 0.0;
      }
      s[1] = v1 / v0;
    }
    else {
      for (int i = 0; i < 29; i++) {
        wv[i] = 0.0;
      }
      s[1] = pow(10.0, -8.0);
    }

    wv = &waves[58];
    if (!in_excision_region_l && !in_excision_region_r) {
      wv[0] = a2 * c_minus;
      wv[4] = a2 * (v0 + (sqrt(s_sq) / y) * ((inv_spacetime_metric_avg[1][0] * v0) - (inv_spacetime_metric_avg[0][0] * v1)));
      wv[1] = a2 * (v1 + (sqrt(s_sq) / y) * ((inv_spacetime_metric_avg[1][1] * v0) - (inv_spacetime_metric_avg[0][1] * v1)));
      wv[2] = a2 * (v2 + (sqrt(s_sq) / y) * ((inv_spacetime_metric_avg[1][2] * v0) - (inv_spacetime_metric_avg[0][2] * v1)));
      wv[3] = a2 * (v3 + (sqrt(s_sq) / y) * ((inv_spacetime_metric_avg[1][3] * v0) - (inv_spacetime_metric_avg[0][3] * v1)));
      for (int i = 5; i < 29; i++) {
        wv[i] = 0.0;
      }
      s[2] = (((1.0 - (gas_gamma * v4)) * (v0 * v1)) - (s_sq * inv_spacetime_metric_avg[0][1]) + (sqrt(s_sq) * y)) /
        (((1.0 - (gas_gamma * v4)) * (v0 * v0)) - (s_sq * inv_spacetime_metric_avg[0][0]));
      }
    else {
      for (int i = 0; i < 29; i++) {
        wv[i] = 0.0;
      }
      s[2] = pow(10.0, -8.0);
    }

    if (!in_excision_region_l && !in_excision_region_r) {
    return (((1.0 - (gas_gamma * v4)) * (v0 * v1)) - (s_sq * inv_spacetime_metric_avg[0][1]) + (sqrt(s_sq) * y)) /
      (((1.0 - (gas_gamma * v4)) * (v0 * v0)) - (s_sq * inv_spacetime_metric_avg[0][0]));
    }
    else {
      return pow(10.0, -8.0);
    }
  }
  else {
    double W_l = 1.0 / sqrt(1.0 - ((vx_l * vx_l) + (vy_l * vy_l) + (vz_l * vz_l)));
    if ((vx_l * vx_l) + (vy_l * vy_l) + (vz_l * vz_l) > 1.0 - pow(10.0, -8.0)) {
      W_l = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double W_r = 1.0 / sqrt(1.0 - ((vx_r * vx_r) + (vy_r * vy_r) + (vz_r * vz_r)));
    if ((vx_r * vx_r) + (vy_r * vy_r) + (vz_r * vz_r) > 1.0 - pow(10.0, -8.0)) {
      W_r = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double K_l = sqrt(ql[4] + p_l) / W_l;
    double K_r = sqrt(qr[4] + p_r) / W_r;
    double inv_K_avg = 1.0 / (K_l + K_r);

    double v0 = ((K_l * W_l) + (K_r * W_r)) * inv_K_avg;
    double v1 = ((K_l * W_l * vx_l) + (K_r * W_r * vx_r)) * inv_K_avg;
    double v2 = ((K_l * W_l * vy_l) + (K_r * W_r * vy_r)) * inv_K_avg;
    double v3 = ((K_l * W_l * vz_l) + (K_r * W_r * vz_r)) * inv_K_avg;
    double v4 = ((p_l / K_l) + (p_r / K_r)) * inv_K_avg;

    double c_minus = 1.0 - (gas_gamma / (gas_gamma - 1.0)) * v4;
    double c_plus = 1.0 + (gas_gamma / (gas_gamma - 1.0)) * v4;

    double v_contr = -(v0 * v0) + (v1 * v1) + (v2 * v2) + (v3 * v3);

    double s_sq = (0.5 * gas_gamma * v4 * (1.0 - v_contr)) - (0.5 * (gas_gamma - 1.0) * (1.0 + v_contr));
    double e = (v0 * v0) - (v1 * v1);
    double y = sqrt(((1.0 - (gas_gamma * v4)) * e) + s_sq);

    double k = (v0 * delta_mod[1]) - (v1 * delta_mod[2]);
    double v_delta_contr = -(v0 * delta_mod[1]) + (v1 * delta_mod[2]) + (v2 * delta_mod[3]) + (v3 * delta_mod[4]);
    
    double a1 = -((s_sq * k) + (sqrt(s_sq) * y * ((v0 * delta_mod[2]) - (v1 * delta_mod[1]))) +
      ((gas_gamma - 1.0) * e * (delta_mod[0] + (c_plus * v_delta_contr)))) / (2.0 * e * s_sq);
    double a2 = -((s_sq * k) - (sqrt(s_sq) * y * ((v0 * delta_mod[2]) - (v1 * delta_mod[1]))) +
      ((gas_gamma - 1.0) * e * (delta_mod[0] + (c_plus * v_delta_contr)))) / (2.0 * e * s_sq);
    double a3 = ((2.0 * s_sq * k) + ((gas_gamma - 1.0) * e * (delta_mod[0] + (c_plus * v_delta_contr)))) / (e * s_sq);
    double a4 = delta_mod[3] - ((k * v2) / e);
    double a5 = delta_mod[4] - ((k * v3) / e);

    double *wv;
    wv = &waves[0];
    wv[0] = a1 * c_minus;
    wv[4] = a1 * (v0 - (sqrt(s_sq) * v1) / y);
    wv[1] = a1 * (v1 - (sqrt(s_sq) * v0) / y);
    wv[2] = a1 * v2;
    wv[3] = a1 * v3;
    for (int i = 5; i < 29; i++) {
      wv[i] = 0.0;
    }
    s[0] = (((1.0 - (gas_gamma * v4)) * v0 * v1) - (sqrt(s_sq) * y)) / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);

    wv = &waves[29];
    wv[0] = (a3 * (c_minus + (s_sq / (gas_gamma - 1.0)))) - (a4 * c_plus * v2) - (a5 * c_plus * v3);
    wv[4] = a3 * v0;
    wv[1] = a3 * v1;
    wv[2] = (a3 * v2) + a4;
    wv[3] = (a3 * v3) + a5;
    for (int i = 5; i < 29; i++) {
      wv[i] = 0.0;
    }
    s[1] = v1 / v0;

    wv = &waves[58];
    wv[0] = a2 * c_minus;
    wv[4] = a2 * (v0 + (sqrt(s_sq) * v1) / y);
    wv[1] = a2 * (v1 + (sqrt(s_sq) * v0) / y);
    wv[2] = a2 * v2;
    wv[3] = a2 * v3;
    for (int i = 5; i < 29; i++) {
      wv[i] = 0.0;
    }
    s[2] = (((1.0 - (gas_gamma * v4)) * v0 * v1) + (sqrt(s_sq) * y)) / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);

    return (((1.0 - (gas_gamma * v4)) * v0 * fabs(v1)) + (sqrt(s_sq) * y)) / (((1.0 - (gas_gamma * v4)) * v0 * v0) + s_sq);
  }
}

static void
qfluct_roe(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double* w0 = &waves[0];
  const double* w1 = &waves[29];
  const double* w2 = &waves[58];

  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]);

  for (int i = 0; i < 5; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]) + (s2m * w2[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]) + (s2p * w2[i]);
  }
  for (int i = 5; i < 29; i++) {
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

double
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

bool
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

double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  double gas_gamma = gr_euler->gas_gamma;

  return gkyl_gr_euler_max_abs_speed(gas_gamma, q);
}

void
gr_euler_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 29; i++) {
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

  if (inp->rp_type == WV_GR_EULER_RP_ROE) {
    gr_euler->eqn.num_waves = 3;
    gr_euler->eqn.waves_func = wave_roe_l;
    gr_euler->eqn.qfluct_func = qfluct_roe_l;
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

struct gkyl_gr_spacetime*
gkyl_wv_gr_euler_spacetime(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  struct gkyl_gr_spacetime* spacetime = gr_euler->spacetime;

  return spacetime;
}