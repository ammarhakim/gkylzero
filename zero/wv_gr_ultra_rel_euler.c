#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_ultra_rel_euler.h>
#include <gkyl_wv_gr_ultra_rel_euler_priv.h>

void gkyl_gr_ultra_rel_euler_prim_vars(double gas_gamma, const double q[27], double v[27])
{
  double lapse = q[4];
  double shift_x = q[5];
  double shift_y = q[6];
  double shift_z = q[7];

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_metric[0][0] = q[8]; spatial_metric[0][1] = q[9]; spatial_metric[0][2] = q[10];
  spatial_metric[1][0] = q[11]; spatial_metric[1][1] = q[12]; spatial_metric[1][2] = q[13];
  spatial_metric[2][0] = q[14]; spatial_metric[2][1] = q[15]; spatial_metric[2][2] = q[16];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_ultra_rel_euler_inv_spatial_metric(q, &inv_spatial_metric);
  
  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
  }

  extrinsic_curvature[0][0] = q[17]; extrinsic_curvature[0][1] = q[18]; extrinsic_curvature[0][2] = q[19];
  extrinsic_curvature[1][0] = q[20]; extrinsic_curvature[1][1] = q[21]; extrinsic_curvature[1][2] = q[22];
  extrinsic_curvature[2][0] = q[23]; extrinsic_curvature[2][1] = q[24]; extrinsic_curvature[2][2] = q[25];

  bool in_excision_region = false;
  if (q[26] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

    double Etot = q[0] / sqrt(spatial_det);
    double momx = q[1] / sqrt(spatial_det);
    double momy = q[2] / sqrt(spatial_det);
    double momz = q[3] / sqrt(spatial_det);

    double *mom = gkyl_malloc(sizeof(double[3]));
    double mom_sq = 0.0;
    mom[0] = momx; mom[1] = momy; mom[2] = momz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        mom_sq += inv_spatial_metric[i][j] * mom[i] * mom[j];
      }
    }

    gkyl_free(mom);

    double beta = 0.25 * (2.0 - gas_gamma);
    double p = -(2.0 * beta * Etot) + sqrt((4.0 * (beta * beta) * (Etot * Etot)) + ((gas_gamma - 1.0) * ((Etot * Etot) - mom_sq)));
    double rho = p / (gas_gamma - 1.0);

    double cov_vx = momx / (Etot + p);
    double cov_vy = momy / (Etot + p);
    double cov_vz = momz / (Etot + p);

    v[0] = rho;
    v[1] = cov_vx;
    v[2] = cov_vy;
    v[3] = cov_vz;

    v[4] = lapse;
    v[5] = shift_x;
    v[6] = shift_y;
    v[7] = shift_z;

    v[8] = spatial_metric[0][0]; v[9] = spatial_metric[0][1]; v[10] = spatial_metric[0][2];
    v[11] = spatial_metric[1][0]; v[12] = spatial_metric[1][1]; v[13] = spatial_metric[1][2];
    v[14] = spatial_metric[2][0]; v[15] = spatial_metric[2][1]; v[16] = spatial_metric[2][2];

    v[17] = extrinsic_curvature[0][0]; v[18] = extrinsic_curvature[0][1]; v[19] = extrinsic_curvature[0][2];
    v[20] = extrinsic_curvature[1][0]; v[21] = extrinsic_curvature[1][1]; v[22] = extrinsic_curvature[1][2];
    v[23] = extrinsic_curvature[2][0]; v[24] = extrinsic_curvature[2][1]; v[25] = extrinsic_curvature[2][2];

    v[26] = 1.0;
  }
  else {
    for (int i = 0; i < 26; i++) {
      v[i] = 0.0;
    }

    v[26] = -1.0;
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(inv_spatial_metric[i]);
    gkyl_free(extrinsic_curvature[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(inv_spatial_metric);
  gkyl_free(extrinsic_curvature);
}

void 
gkyl_gr_ultra_rel_euler_inv_spatial_metric(const double q[27], double ***inv_spatial_metric)
{
  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_metric[0][0] = q[8]; spatial_metric[0][1] = q[9]; spatial_metric[0][2] = q[10];
  spatial_metric[1][0] = q[11]; spatial_metric[1][1] = q[12]; spatial_metric[1][2] = q[13];
  spatial_metric[2][0] = q[14]; spatial_metric[2][1] = q[15]; spatial_metric[2][2] = q[16];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += spatial_metric[i][i];
  }

  double **spatial_metric_sq = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_sq[i] = gkyl_malloc(sizeof(double[3]));
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_metric_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_sq[i][j] += spatial_metric[i][k] * spatial_metric[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += spatial_metric_sq[i][i];
  }

  double **euclidean_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    euclidean_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        euclidean_metric[i][j] = 1.0;
      }
      else {
        euclidean_metric[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*inv_spatial_metric)[i][j] = (1.0 / spatial_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(spatial_metric_sq[i]);
    gkyl_free(euclidean_metric[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(spatial_metric_sq);
  gkyl_free(euclidean_metric);
}

static inline double
gkyl_gr_ultra_rel_euler_max_abs_speed(double gas_gamma, const double q[27])
{
  double v[27] = { 0.0 };
  gkyl_gr_ultra_rel_euler_prim_vars(gas_gamma, q, v);
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];

  double lapse = v[4];
  double shift_x = v[5];
  double shift_y = v[6];
  double shift_z = v[7];

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_metric[0][0] = v[8]; spatial_metric[0][1] = v[9]; spatial_metric[0][2] = v[10];
  spatial_metric[1][0] = v[11]; spatial_metric[1][1] = v[12]; spatial_metric[1][2] = v[13];
  spatial_metric[2][0] = v[14]; spatial_metric[2][1] = v[15]; spatial_metric[2][2] = v[16];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_ultra_rel_euler_inv_spatial_metric(q, &inv_spatial_metric);

  double c_s = gas_gamma - 1.0;

  bool in_excision_region = false;
  if (v[26] < pow(10.0, -8.0)) {
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
  if (fabs(lapse - 1.0) > pow(10.0, -8.0) || fabs(shift_x) > pow(10.0, -8.0) || fabs(shift_y) > pow(10.0, -8.0) ||
    fabs(shift_z) > pow(10.0, -8.0)) {
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

void
gkyl_gr_ultra_rel_euler_flux(double gas_gamma, const double q[27], double flux[27])
{
  double v[27] = { 0.0 };
  gkyl_gr_ultra_rel_euler_prim_vars(gas_gamma, q, v);
  double rho =  v[0];
  double vx = v[1];
  double vy = v[2];
  double vz = v[3];
  double p = (gas_gamma - 1.0) * rho;

  double lapse = v[4];
  double shift_x = v[5];

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_metric[0][0] = v[8]; spatial_metric[0][1] = v[9]; spatial_metric[0][2] = v[10];
  spatial_metric[1][0] = v[11]; spatial_metric[1][1] = v[12]; spatial_metric[1][2] = v[13];
  spatial_metric[2][0] = v[14]; spatial_metric[2][1] = v[15]; spatial_metric[2][2] = v[16];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  bool in_excision_region = false;
  if (v[26] < pow(10.0, -8.0)) {
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

    flux[0] = (lapse * sqrt(spatial_det)) * ((((rho + p) * (W * W)) - p) * (vx - (shift_x / lapse)) + (p * vx));
    flux[1] = (lapse * sqrt(spatial_det)) * ((rho + p) * (W * W) * (vx * (vx - (shift_x / lapse))) + p);
    flux[2] = (lapse * sqrt(spatial_det)) * ((rho + p) * (W * W) * (vy * (vx - (shift_x / lapse))));
    flux[3] = (lapse * sqrt(spatial_det)) * ((rho + p) * (W * W) * (vz * (vx - (shift_x / lapse))));

    for (int i = 4; i < 27; i++) {
      flux[i] = 0.0;
    }

    gkyl_free(vel);
  }
  else {
    for (int i = 0; i < 27; i++) {
      flux[i] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
  }
  gkyl_free(spatial_metric);
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 27; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 27; i++) {
    qout[i] = win[i];
  }
}

static void
gr_ultra_rel_euler_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 27; i++) {
    ghost[i] = skin[i];
  }

  ghost[1] = -ghost[1];
}

static void
gr_ultra_rel_euler_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 1; i < 4; i++) {
    ghost[i] = -skin[i];
  }

  ghost[0] = skin[0];

  for (int i = 4; i < 27; i++) {
    ghost[i] = skin[i];
  }
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  qlocal[0] = qglobal[0];
  qlocal[1] = (qglobal[1] * norm[0]) + (qglobal[2] * norm[1]) + (qglobal[3] * norm[2]);
  qlocal[2] = (qglobal[1] * tau1[0]) + (qglobal[2] * tau1[1]) + (qglobal[3] * tau1[2]);
  qlocal[3] = (qglobal[1] * tau2[0]) + (qglobal[2] * tau2[1]) + (qglobal[3] * tau2[2]);

  qlocal[4] = qglobal[4];
  qlocal[5] = (qglobal[5] * norm[0]) + (qglobal[6] * norm[1]) + (qglobal[7] * norm[2]);
  qlocal[6] = (qglobal[5] * tau1[0]) + (qglobal[6] * tau1[1]) + (qglobal[7] * tau1[2]);
  qlocal[7] = (qglobal[5] * tau2[0]) + (qglobal[6] * tau2[1]) + (qglobal[7] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qglobal[8] * norm[0]) + (qglobal[9] * norm[1]) + (qglobal[10] * norm[2]);
  r1[1] = (qglobal[8] * tau1[0]) + (qglobal[9] * tau1[1]) + (qglobal[10] * tau1[2]);
  r1[2] = (qglobal[8] * tau2[0]) + (qglobal[9] * tau2[1]) + (qglobal[10] * tau2[2]);

  r2[0] = (qglobal[11] * norm[0]) + (qglobal[12] * norm[1]) + (qglobal[13] * norm[2]);
  r2[1] = (qglobal[11] * tau1[0]) + (qglobal[12] * tau1[1]) + (qglobal[13] * tau1[2]);
  r2[2] = (qglobal[11] * tau2[0]) + (qglobal[12] * tau2[1]) + (qglobal[13] * tau2[2]);

  r3[0] = (qglobal[14] * norm[0]) + (qglobal[15] * norm[1]) + (qglobal[16] * norm[2]);
  r3[1] = (qglobal[14] * tau1[0]) + (qglobal[15] * tau1[1]) + (qglobal[16] * tau1[2]);
  r3[2] = (qglobal[14] * tau2[0]) + (qglobal[15] * tau2[1]) + (qglobal[16] * tau2[2]);

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

  // Rotate spatial metric tensor to local coordinate frame.
  qlocal[8] = v1[0]; qlocal[9] = v1[1]; qlocal[10] = v1[2];
  qlocal[11] = v2[0]; qlocal[12] = v2[1]; qlocal[13] = v2[2];
  qlocal[14] = v3[0]; qlocal[15] = v3[1]; qlocal[16] = v3[2];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qglobal[17] * norm[0]) + (qglobal[18] * norm[1]) + (qglobal[19] * norm[2]);
  extr_r1[1] = (qglobal[17] * tau1[0]) + (qglobal[18] * tau1[1]) + (qglobal[19] * tau1[2]);
  extr_r1[2] = (qglobal[17] * tau2[0]) + (qglobal[18] * tau2[1]) + (qglobal[19] * tau2[2]);

  extr_r2[0] = (qglobal[20] * norm[0]) + (qglobal[21] * norm[1]) + (qglobal[22] * norm[2]);
  extr_r2[1] = (qglobal[20] * tau1[0]) + (qglobal[21] * tau1[1]) + (qglobal[22] * tau1[2]);
  extr_r2[2] = (qglobal[20] * tau2[0]) + (qglobal[21] * tau2[1]) + (qglobal[22] * tau2[2]);

  extr_r3[0] = (qglobal[23] * norm[0]) + (qglobal[24] * norm[1]) + (qglobal[25] * norm[2]);
  extr_r3[1] = (qglobal[23] * tau1[0]) + (qglobal[24] * tau1[1]) + (qglobal[25] * tau1[2]);
  extr_r3[2] = (qglobal[23] * tau2[0]) + (qglobal[24] * tau2[1]) + (qglobal[25] * tau2[2]);

  // Temporary arrays to store rotated extrinsic row vectors.
  double inv_v1[3], inv_v2[3], inv_v3[3];
  inv_v1[0] = (extr_r1[0] * norm[0]) + (extr_r2[0] * norm[1]) + (extr_r3[0] * norm[2]);
  inv_v1[1] = (extr_r1[0] * tau1[0]) + (extr_r2[0] * tau1[1]) + (extr_r3[0] * tau1[2]);
  inv_v1[2] = (extr_r1[0] * tau2[0]) + (extr_r2[0] * tau2[1]) + (extr_r3[0] * tau2[2]);

  inv_v2[0] = (extr_r1[1] * norm[0]) + (extr_r2[1] * norm[1]) + (extr_r3[1] * norm[2]);
  inv_v2[1] = (extr_r1[1] * tau1[0]) + (extr_r2[1] * tau1[1]) + (extr_r3[1] * tau1[2]);
  inv_v2[2] = (extr_r1[1] * tau2[0]) + (extr_r2[1] * tau2[1]) + (extr_r3[1] * tau2[2]);

  inv_v3[0] = (extr_r1[2] * norm[0]) + (extr_r2[2] * norm[1]) + (extr_r3[2] * norm[2]);
  inv_v3[1] = (extr_r1[2] * tau1[0]) + (extr_r2[2] * tau1[1]) + (extr_r3[2] * tau1[2]);
  inv_v3[2] = (extr_r1[2] * tau2[0]) + (extr_r2[2] * tau2[1]) + (extr_r3[2] * tau2[2]);

  // Rotate extrinsic curvature tensor to local coordinate frame.
  qlocal[17] = inv_v1[0]; qlocal[18] = inv_v1[1]; qlocal[19] = inv_v1[2];
  qlocal[20] = inv_v2[0]; qlocal[21] = inv_v2[1]; qlocal[22] = inv_v2[2];
  qlocal[23] = inv_v3[0]; qlocal[24] = inv_v3[1]; qlocal[25] = inv_v3[2];

  qlocal[26] = qglobal[26];
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  qglobal[0] = qlocal[0];
  qglobal[1] = (qlocal[1] * norm[0]) + (qlocal[2] * tau1[0]) + (qlocal[3] * tau2[0]);
  qglobal[2] = (qlocal[1] * norm[1]) + (qlocal[2] * tau1[1]) + (qlocal[3] * tau2[1]);
  qglobal[3] = (qlocal[1] * norm[2]) + (qlocal[2] * tau1[2]) + (qlocal[3] * tau2[2]);

  qglobal[4] = qlocal[4];
  qglobal[5] = (qlocal[5] * norm[0]) + (qlocal[6] * tau1[0]) + (qlocal[7] * tau2[0]);
  qglobal[6] = (qlocal[5] * norm[1]) + (qlocal[6] * tau1[1]) + (qlocal[7] * tau2[1]);
  qglobal[7] = (qlocal[5] * norm[2]) + (qlocal[6] * tau1[2]) + (qlocal[7] * tau2[2]);

  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qlocal[8] * norm[0]) + (qlocal[9] * tau1[0]) + (qlocal[10] * tau2[0]);
  r1[1] = (qlocal[8] * norm[1]) + (qlocal[9] * tau1[1]) + (qlocal[10] * tau2[1]);
  r1[2] = (qlocal[8] * norm[2]) + (qlocal[9] * tau1[2]) + (qlocal[10] * tau2[2]);

  r2[0] = (qlocal[11] * norm[0]) + (qlocal[12] * tau1[0]) + (qlocal[13] * tau2[0]);
  r2[1] = (qlocal[11] * norm[1]) + (qlocal[12] * tau1[1]) + (qlocal[13] * tau2[1]);
  r2[2] = (qlocal[11] * norm[2]) + (qlocal[12] * tau1[2]) + (qlocal[13] * tau2[2]);

  r3[0] = (qlocal[14] * norm[0]) + (qlocal[15] * tau1[0]) + (qlocal[16] * tau2[0]);
  r3[1] = (qlocal[14] * norm[1]) + (qlocal[15] * tau1[1]) + (qlocal[16] * tau2[1]);
  r3[2] = (qlocal[14] * norm[2]) + (qlocal[15] * tau1[2]) + (qlocal[16] * tau2[2]);

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

  // Rotate spatial metric tensor back to global coordinate frame.
  qglobal[8] = v1[0]; qglobal[9] = v1[1]; qglobal[10] = v1[2];
  qglobal[11] = v2[0]; qglobal[12] = v2[1]; qglobal[13] = v2[2];
  qglobal[14] = v3[0]; qglobal[15] = v3[1]; qglobal[16] = v3[2];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qlocal[17] * norm[0]) + (qlocal[18] * tau1[0]) + (qlocal[19] * tau2[0]);
  extr_r1[1] = (qlocal[17] * norm[1]) + (qlocal[18] * tau1[1]) + (qlocal[19] * tau2[1]);
  extr_r1[2] = (qlocal[17] * norm[2]) + (qlocal[18] * tau1[2]) + (qlocal[19] * tau2[2]);

  extr_r2[0] = (qlocal[20] * norm[0]) + (qlocal[21] * tau1[0]) + (qlocal[22] * tau2[0]);
  extr_r2[1] = (qlocal[20] * norm[1]) + (qlocal[21] * tau1[1]) + (qlocal[22] * tau2[1]);
  extr_r2[2] = (qlocal[20] * norm[2]) + (qlocal[21] * tau1[2]) + (qlocal[22] * tau2[2]);

  extr_r3[0] = (qlocal[23] * norm[0]) + (qlocal[24] * tau1[0]) + (qlocal[25] * tau2[0]);
  extr_r3[1] = (qlocal[23] * norm[1]) + (qlocal[24] * tau1[1]) + (qlocal[25] * tau2[1]);
  extr_r3[2] = (qlocal[23] * norm[2]) + (qlocal[24] * tau1[2]) + (qlocal[25] * tau2[2]);

  // Temporary arrays to store rotated extrinsic row vectors.
  double inv_v1[3], inv_v2[3], inv_v3[3];
  inv_v1[0] = (extr_r1[0] * norm[0]) + (extr_r2[0] * tau1[0]) + (extr_r3[0] * tau2[0]);
  inv_v1[1] = (extr_r1[0] * norm[1]) + (extr_r2[0] * tau1[1]) + (extr_r3[0] * tau2[1]);
  inv_v1[2] = (extr_r1[0] * norm[2]) + (extr_r2[0] * tau1[2]) + (extr_r3[0] * tau2[2]);

  inv_v2[0] = (extr_r1[1] * norm[0]) + (extr_r2[1] * tau1[0]) + (extr_r3[1] * tau2[0]);
  inv_v2[1] = (extr_r1[1] * norm[1]) + (extr_r2[1] * tau1[1]) + (extr_r3[1] * tau2[1]);
  inv_v2[2] = (extr_r1[1] * norm[2]) + (extr_r2[1] * tau1[2]) + (extr_r3[1] * tau2[2]);

  inv_v3[0] = (extr_r1[2] * norm[0]) + (extr_r2[2] * tau1[0]) + (extr_r3[2] * tau2[0]);
  inv_v3[1] = (extr_r1[2] * norm[1]) + (extr_r2[2] * tau1[1]) + (extr_r3[2] * tau2[1]);
  inv_v3[2] = (extr_r1[2] * norm[2]) + (extr_r2[2] * tau1[2]) + (extr_r3[2] * tau2[2]);

  // Rotate extrinsic curvature tensor back to global coordinate frame.
  qglobal[17] = inv_v1[0]; qglobal[18] = inv_v1[1]; qglobal[19] = inv_v1[2];
  qglobal[20] = inv_v2[0]; qglobal[21] = inv_v2[1]; qglobal[22] = inv_v2[2];
  qglobal[23] = inv_v3[0]; qglobal[24] = inv_v3[1]; qglobal[25] = inv_v3[2];

  qglobal[26] = qlocal[26];
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(eqn, struct wv_gr_ultra_rel_euler, eqn);
  double gas_gamma = gr_ultra_rel_euler->gas_gamma;

  double sl = gkyl_gr_ultra_rel_euler_max_abs_speed(gas_gamma, ql);
  double sr = gkyl_gr_ultra_rel_euler_max_abs_speed(gas_gamma, qr);
  double amax = fmax(sl, sr);

  double fl[27], fr[27];
  gkyl_gr_ultra_rel_euler_flux(gas_gamma, ql, fl);
  gkyl_gr_ultra_rel_euler_flux(gas_gamma, qr, fr);

  bool in_excision_region_l = false;
  if (ql[26] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[26] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  double *w0 = &waves[0], *w1 = &waves[27];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 27; i++) {
      w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
      w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
    }
  }
  else {
    for (int i = 0; i < 27; i++) {
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
  const double *w0 = &waves[0], *w1 = &waves[27];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 27; i++) {
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
  const struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(eqn, struct wv_gr_ultra_rel_euler, eqn);

  double fr[27], fl[27];
  gkyl_gr_ultra_rel_euler_flux(gr_ultra_rel_euler->gas_gamma, ql, fl);
  gkyl_gr_ultra_rel_euler_flux(gr_ultra_rel_euler->gas_gamma, qr, fr);

  bool in_excision_region_l = false;
  if (ql[26] < pow(10.0, -8.0)) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[26] < pow(10.0, -8.0)) {
    in_excision_region_r = true;
  }

  if (!in_excision_region_l && !in_excision_region_r) {
    for (int m = 0; m < 27; m++) {
      flux_jump[m] = fr[m] - fl[m];
    }
  }
  else {
    for (int m = 0; m < 27; m++) {
      flux_jump[m] = 0.0;
    }
  }

  double amaxl = gkyl_gr_ultra_rel_euler_max_abs_speed(gr_ultra_rel_euler->gas_gamma, ql);
  double amaxr = gkyl_gr_ultra_rel_euler_max_abs_speed(gr_ultra_rel_euler->gas_gamma, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(eqn, struct wv_gr_ultra_rel_euler, eqn);
  double gas_gamma = gr_ultra_rel_euler->gas_gamma;

  double v[27] = { 0.0 };
  gkyl_gr_ultra_rel_euler_prim_vars(gas_gamma, q, v);

  if (v[0] < 0.0) {
    return false;
  }
  else {
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(eqn, struct wv_gr_ultra_rel_euler, eqn);
  double gas_gamma = gr_ultra_rel_euler->gas_gamma;

  return gkyl_gr_ultra_rel_euler_max_abs_speed(gas_gamma, q);
}

static inline void
gr_ultra_rel_euler_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 27; i++) {
    diag[i] = qin[i];
  }
}

static inline void
gr_ultra_rel_euler_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  // TODO: Rewrite source solver(s).
  for (int i = 0; i < 27; i++) {
    sout[i] = 0.0;
  }
}

void
gkyl_gr_ultra_rel_euler_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(base->on_dev, struct wv_gr_ultra_rel_euler, eqn);
    gkyl_cu_free(gr_ultra_rel_euler);
  }

  struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(base, struct wv_gr_ultra_rel_euler, eqn);
  gkyl_free(gr_ultra_rel_euler);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_ultra_rel_euler_new(double gas_gamma, struct gkyl_gr_spacetime* spacetime, bool use_gpu)
{
  return gkyl_wv_gr_ultra_rel_euler_inew(&(struct gkyl_wv_gr_ultra_rel_euler_inp) {
      .gas_gamma = gas_gamma,
      .spacetime = spacetime,
      .rp_type = WV_GR_ULTRA_REL_EULER_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_ultra_rel_euler_inew(const struct gkyl_wv_gr_ultra_rel_euler_inp* inp)
{
  struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = gkyl_malloc(sizeof(struct wv_gr_ultra_rel_euler));

  gr_ultra_rel_euler->eqn.type = GKYL_EQN_GR_ULTRA_REL_EULER;
  gr_ultra_rel_euler->eqn.num_equations = 27;
  gr_ultra_rel_euler->eqn.num_diag = 5;

  gr_ultra_rel_euler->gas_gamma = inp->gas_gamma;
  gr_ultra_rel_euler->spacetime = inp->spacetime;

  if (inp->rp_type == WV_GR_ULTRA_REL_EULER_RP_LAX) {
    gr_ultra_rel_euler->eqn.num_waves = 2;
    gr_ultra_rel_euler->eqn.waves_func = wave_lax_l;
    gr_ultra_rel_euler->eqn.qfluct_func = qfluct_lax_l;
  }

  gr_ultra_rel_euler->eqn.flux_jump = flux_jump;
  gr_ultra_rel_euler->eqn.check_inv_func = check_inv;
  gr_ultra_rel_euler->eqn.max_speed_func = max_speed;
  gr_ultra_rel_euler->eqn.rotate_to_local_func = rot_to_local;
  gr_ultra_rel_euler->eqn.rotate_to_global_func = rot_to_global;

  gr_ultra_rel_euler->eqn.wall_bc_func = gr_ultra_rel_euler_wall;
  gr_ultra_rel_euler->eqn.no_slip_bc_func = gr_ultra_rel_euler_no_slip;

  gr_ultra_rel_euler->eqn.cons_to_riem = cons_to_riem;
  gr_ultra_rel_euler->eqn.riem_to_cons = riem_to_cons;

  gr_ultra_rel_euler->eqn.cons_to_diag = gr_ultra_rel_euler_cons_to_diag;

  gr_ultra_rel_euler->eqn.source_func = gr_ultra_rel_euler_source;

  gr_ultra_rel_euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_ultra_rel_euler->eqn.flags);
  gr_ultra_rel_euler->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_ultra_rel_euler_free);
  gr_ultra_rel_euler->eqn.on_dev = &gr_ultra_rel_euler->eqn; // On the CPU, the equation object points to itself.

  return &gr_ultra_rel_euler->eqn;
}

double
gkyl_wv_gr_ultra_rel_euler_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(eqn, struct wv_gr_ultra_rel_euler, eqn);
  double gas_gamma = gr_ultra_rel_euler->gas_gamma;

  return gas_gamma;
}

struct gkyl_gr_spacetime*
gkyl_wv_gr_ultra_rel_euler_spacetime(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(eqn, struct wv_gr_ultra_rel_euler, eqn);
  struct gkyl_gr_spacetime *spacetime = gr_ultra_rel_euler->spacetime;

  return spacetime;
}