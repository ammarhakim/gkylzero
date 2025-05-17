#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_twofluid.h>
#include <gkyl_wv_gr_twofluid_priv.h>

void
gkyl_gr_twofluid_flux(double gas_gamma_elc, double gas_gamma_ion, double light_speed, double e_fact, double b_fact, const double q[84], double flux[84])
{
  double v[84] = { 0.0 };
  gkyl_gr_twofluid_prim_vars(gas_gamma_elc, gas_gamma_ion, q, v);
  double rho_elc = v[0];
  double vx_elc = v[1];
  double vy_elc = v[2];
  double vz_elc = v[3];
  double p_elc = v[4];

  double rho_ion = v[5];
  double vx_ion = v[6];
  double vy_ion = v[7];
  double vz_ion = v[8];
  double p_ion = v[9];

  double Dx = v[10], Dy = v[11], Dz = v[12];
  double Bx = v[13], By = v[14], Bz = v[15];

  double phi = v[16];
  double psi = v[17];

  double lapse = v[18];
  double shift_x = v[19];
  double shift_y = v[20];
  double shift_z = v[21];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[22]; spatial_metric[0][1] = v[23]; spatial_metric[0][2] = v[24];
  spatial_metric[1][0] = v[25]; spatial_metric[1][1] = v[26]; spatial_metric[1][2] = v[27];
  spatial_metric[2][0] = v[28]; spatial_metric[2][1] = v[29]; spatial_metric[2][2] = v[30];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  bool in_excision_region = false;
  if (v[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel_elc[3];
    double v_sq_elc = 0.0;
    vel_elc[0] = vx_elc; vel_elc[1] = vy_elc; vel_elc[2] = vz_elc;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_elc += spatial_metric[i][j] * vel_elc[i] * vel_elc[j];
      }
    }

    double W_elc = 1.0 / (sqrt(1.0 - v_sq_elc));
    if (v_sq_elc > 1.0 - pow(10.0, -8.0)) {
      W_elc = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double h_elc = 1.0 + ((p_elc / rho_elc) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));

    flux[0] = (lapse * sqrt(spatial_det)) * (rho_elc * W_elc * (vx_elc - (shift_x / lapse)));
    flux[1] = (lapse * sqrt(spatial_det)) * (rho_elc * h_elc * (W_elc * W_elc) * (vx_elc * (vx_elc - (shift_x / lapse))) + p_elc);
    flux[2] = (lapse * sqrt(spatial_det)) * (rho_elc * h_elc * (W_elc * W_elc) * (vy_elc * (vx_elc - (shift_x / lapse))));
    flux[3] = (lapse * sqrt(spatial_det)) * (rho_elc * h_elc * (W_elc * W_elc) * (vz_elc * (vx_elc - (shift_x / lapse))));
    flux[4] = (lapse * sqrt(spatial_det)) * (((rho_elc * h_elc * (W_elc * W_elc)) - p_elc - (rho_elc * W_elc)) * (vx_elc - (shift_x / lapse)) + (p_elc * vx_elc));

    double vel_ion[3];
    double v_sq_ion = 0.0;
    vel_ion[0] = vx_ion; vel_ion[1] = vy_ion; vel_ion[2] = vz_ion;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_ion += spatial_metric[i][j] * vel_ion[i] * vel_ion[j];
      }
    }

    double W_ion = 1.0 / (sqrt(1.0 - v_sq_ion));
    if (v_sq_ion > 1.0 - pow(10.0, -8.0)) {
      W_ion = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double h_ion = 1.0 + ((p_ion / rho_ion) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));

    flux[5] = (lapse * sqrt(spatial_det)) * (rho_ion * W_ion * (vx_ion - (shift_x / lapse)));
    flux[6] = (lapse * sqrt(spatial_det)) * (rho_ion * h_ion * (W_ion * W_ion) * (vx_ion * (vx_ion - (shift_x / lapse))) + p_ion);
    flux[7] = (lapse * sqrt(spatial_det)) * (rho_ion * h_ion * (W_ion * W_ion) * (vy_ion * (vx_ion - (shift_x / lapse))));
    flux[8] = (lapse * sqrt(spatial_det)) * (rho_ion * h_ion * (W_ion * W_ion) * (vz_ion * (vx_ion - (shift_x / lapse))));
    flux[9] = (lapse * sqrt(spatial_det)) * (((rho_ion * h_ion * (W_ion * W_ion)) - p_ion - (rho_ion * W_ion)) * (vx_ion - (shift_x / lapse)) + (p_ion * vx_ion));

    double Ex = (lapse * Dx) + ((shift_y * Bz) - (shift_z * By));
    double Ey = (lapse * Dy) - ((shift_x * Bz) - (shift_z * Bx));
    double Ez = (lapse * Dz) + ((shift_x * By) - (shift_y * Bx));

    double Hx = (lapse * Bx) - ((shift_y * Dz) - (shift_z * Dy));
    double Hy = (lapse * By) + ((shift_x * Dz) - (shift_z * Dx));
    double Hz = (lapse * Bz) - ((shift_x * Dy) - (shift_y * Dx));

    flux[10] = e_fact * (light_speed * light_speed) * phi;
    flux[11] = (light_speed * light_speed) * Hz;
    flux[12] = -(light_speed * light_speed) * Hy;
    flux[13] = b_fact * psi;
    flux[14] = -Ez;
    flux[15] = Ey;
    flux[16] = e_fact * Dx;
    flux[17] = b_fact * (light_speed * light_speed) * Bx;

    for (int i = 18; i < 84; i++) {
      flux[i] = 0.0;
    }
  }
  else {
    for (int i = 0; i < 84; i++) {
      flux[i] = 0.0;
    }
  }
}

void
gkyl_gr_twofluid_prim_vars(double gas_gamma_elc, double gas_gamma_ion, const double q[84], double v[84])
{
  double Dx = q[10], Dy = q[11], Dz = q[12];
  double Bx = q[13], By = q[14], Bz = q[15];

  double phi = q[16];
  double psi = q[17];

  double lapse = q[18];
  double shift_x = q[19];
  double shift_y = q[20];
  double shift_z = q[21];

  double spatial_metric[3][3];
  spatial_metric[0][0] = q[22]; spatial_metric[0][1] = q[23]; spatial_metric[0][2] = q[24];
  spatial_metric[1][0] = q[25]; spatial_metric[1][1] = q[26]; spatial_metric[1][2] = q[27];
  spatial_metric[2][0] = q[28]; spatial_metric[2][1] = q[29]; spatial_metric[2][2] = q[30];
  
  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = q[31]; extrinsic_curvature[0][1] = q[32]; extrinsic_curvature[0][2] = q[33];
  extrinsic_curvature[1][0] = q[34]; extrinsic_curvature[1][1] = q[35]; extrinsic_curvature[1][2] = q[36];
  extrinsic_curvature[2][0] = q[37]; extrinsic_curvature[2][1] = q[38]; extrinsic_curvature[2][2] = q[39];

  double lapse_der[3];
  lapse_der[0] = q[41];
  lapse_der[1] = q[42];
  lapse_der[2] = q[43];

  double shift_der[3][3];
  shift_der[0][0] = q[44]; shift_der[0][1] = q[45]; shift_der[0][2] = q[46];
  shift_der[1][0] = q[47]; shift_der[1][1] = q[48]; shift_der[1][2] = q[49];
  shift_der[2][0] = q[50]; shift_der[2][1] = q[51]; shift_der[2][2] = q[52];

  double spatial_metric_der[3][3][3];
  spatial_metric_der[0][0][0] = q[53]; spatial_metric_der[0][0][1] = q[54]; spatial_metric_der[0][0][2] = q[55];
  spatial_metric_der[0][1][0] = q[56]; spatial_metric_der[0][1][1] = q[57]; spatial_metric_der[0][1][2] = q[58];
  spatial_metric_der[0][2][0] = q[59]; spatial_metric_der[0][2][1] = q[60]; spatial_metric_der[0][2][2] = q[61];

  spatial_metric_der[1][0][0] = q[62]; spatial_metric_der[1][0][1] = q[63]; spatial_metric_der[1][0][2] = q[64];
  spatial_metric_der[1][1][0] = q[65]; spatial_metric_der[1][1][1] = q[66]; spatial_metric_der[1][1][2] = q[67];
  spatial_metric_der[1][2][0] = q[68]; spatial_metric_der[1][2][1] = q[69]; spatial_metric_der[1][2][2] = q[70];

  spatial_metric_der[0][0][0] = q[71]; spatial_metric_der[0][0][1] = q[72]; spatial_metric_der[0][0][2] = q[73];
  spatial_metric_der[0][1][0] = q[74]; spatial_metric_der[0][1][1] = q[75]; spatial_metric_der[0][1][2] = q[76];
  spatial_metric_der[0][2][0] = q[77]; spatial_metric_der[0][2][1] = q[78]; spatial_metric_der[0][2][2] = q[79];

  double evol_param = q[80];
  double x = q[81];
  double y = q[82];
  double z = q[83];

  bool in_excision_region = false;
  if (q[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

    double D_elc = q[0] / sqrt(spatial_det);
    double momx_elc = q[1] / sqrt(spatial_det);
    double momy_elc = q[2] / sqrt(spatial_det);
    double momz_elc = q[3] / sqrt(spatial_det);
    double Etot_elc = q[4] / sqrt(spatial_det);

    double C_elc = D_elc / sqrt(((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) + (momz_elc * momz_elc)));
    double C0_elc = (D_elc + Etot_elc) / sqrt(((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) + (momz_elc * momz_elc)));
    if (((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) + (momz_elc * momz_elc)) < pow(10.0, -8.0)) {
      C_elc = D_elc / sqrt(pow(10.0, -8.0));
      C0_elc = (D_elc + Etot_elc) / sqrt(pow(10.0, -8.0));
    }

    double alpha0_elc = -1.0 / (gas_gamma_elc * gas_gamma_elc);
    double alpha1_elc = -2.0 * C_elc * ((gas_gamma_elc - 1.0) / (gas_gamma_elc * gas_gamma_elc));
    double alpha2_elc = ((gas_gamma_elc - 2.0) / gas_gamma_elc) * ((C0_elc * C0_elc) - 1.0) + 1.0 - (C_elc * C_elc) * ((gas_gamma_elc - 1.0) / gas_gamma_elc) *
      ((gas_gamma_elc - 1.0) / gas_gamma_elc);
    double alpha4_elc = (C0_elc * C0_elc) - 1.0;
    double eta_elc = 2.0 * C_elc *((gas_gamma_elc - 1.0) / gas_gamma_elc);

    double guess_elc = 1.0;
    int iter_elc = 0;

    while (iter_elc < 100) {
      double poly_elc = (alpha4_elc * (guess_elc * guess_elc * guess_elc) * (guess_elc - eta_elc)) + (alpha2_elc * (guess_elc * guess_elc)) +
        (alpha1_elc * guess_elc) + alpha0_elc;
      double poly_der_elc = alpha1_elc + (2.0 * alpha2_elc * guess_elc) + (4.0 * alpha4_elc * (guess_elc * guess_elc * guess_elc)) -
        (3.0 * eta_elc * alpha4_elc * (guess_elc * guess_elc));

      double guess_new_elc = guess_elc - (poly_elc / poly_der_elc);

      if (fabs(guess_elc - guess_new_elc) < pow(10.0, -8.0)) {
        iter_elc = 100;
      }
      else {
        iter_elc += 1;
        guess_elc = guess_new_elc;
      }
    }

    double W_elc = 0.5 * C0_elc * guess_elc * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma_elc - 1.0) / gas_gamma_elc) * ((1.0 - (C_elc * guess_elc)) /
      ((C0_elc * C0_elc) * (guess_elc * guess_elc))))));
    double h_elc = 1.0 / (C_elc * guess_elc);

    v[0] = D_elc / W_elc; 
    v[1] = momx_elc / (v[0] * h_elc * (W_elc * W_elc));
    v[2] = momy_elc / (v[0] * h_elc * (W_elc * W_elc));
    v[3] = momz_elc / (v[0] * h_elc * (W_elc * W_elc));
    v[4] = (v[0] * h_elc * (W_elc * W_elc)) - D_elc - Etot_elc;

    if (v[0] < pow(10.0, -8.0)) {
      v[0] = pow(10.0, -8.0);
    }
    if (v[4] < pow(10.0, -8.0)) {
      v[4] = pow(10.0, -8.0);
    }

    double D_ion = q[5] / sqrt(spatial_det);
    double momx_ion = q[6] / sqrt(spatial_det);
    double momy_ion = q[7] / sqrt(spatial_det);
    double momz_ion = q[8] / sqrt(spatial_det);
    double Etot_ion = q[9] / sqrt(spatial_det);

    double C_ion = D_ion / sqrt(((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) + (momz_ion * momz_ion)));
    double C0_ion = (D_ion + Etot_ion) / sqrt(((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) + (momz_ion * momz_ion)));
    if (((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) + (momz_ion * momz_ion)) < pow(10.0, -8.0)) {
      C_ion = D_ion / sqrt(pow(10.0, -8.0));
      C0_ion = (D_ion + Etot_ion) / sqrt(pow(10.0, -8.0));
    }

    double alpha0_ion = -1.0 / (gas_gamma_ion * gas_gamma_ion);
    double alpha1_ion = -2.0 * C_ion * ((gas_gamma_ion - 1.0) / (gas_gamma_ion * gas_gamma_ion));
    double alpha2_ion = ((gas_gamma_ion - 2.0) / gas_gamma_ion) * ((C0_ion * C0_ion) - 1.0) + 1.0 - (C_ion * C_ion) * ((gas_gamma_ion - 1.0) / gas_gamma_ion) *
      ((gas_gamma_ion - 1.0) / gas_gamma_ion);
    double alpha4_ion = (C0_ion * C0_ion) - 1.0;
    double eta_ion = 2.0 * C_ion *((gas_gamma_ion - 1.0) / gas_gamma_ion);

    double guess_ion = 1.0;
    int iter_ion = 0;

    while (iter_ion < 100) {
      double poly_ion = (alpha4_ion * (guess_ion * guess_ion * guess_ion) * (guess_ion - eta_ion)) + (alpha2_ion * (guess_ion * guess_ion)) +
        (alpha1_ion * guess_ion) + alpha0_ion;
      double poly_der_ion = alpha1_ion + (2.0 * alpha2_ion * guess_ion) + (4.0 * alpha4_ion * (guess_ion * guess_ion * guess_ion)) -
        (3.0 * eta_ion * alpha4_ion * (guess_ion * guess_ion));

      double guess_new_ion = guess_ion - (poly_ion / poly_der_ion);

      if (fabs(guess_ion - guess_new_ion) < pow(10.0, -8.0)) {
        iter_ion = 100;
      }
      else {
        iter_ion += 1;
        guess_ion = guess_new_ion;
      }
    }

    double W_ion = 0.5 * C0_ion * guess_ion * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma_ion - 1.0) / gas_gamma_ion) * ((1.0 - (C_ion * guess_ion)) /
      ((C0_ion * C0_ion) * (guess_ion * guess_ion))))));
    double h_ion = 1.0 / (C_ion * guess_ion);

    v[5] = D_ion / W_ion; 
    v[6] = momx_ion / (v[5] * h_ion * (W_ion * W_ion));
    v[7] = momy_ion / (v[5] * h_ion * (W_ion * W_ion));
    v[8] = momz_ion / (v[5] * h_ion * (W_ion * W_ion));
    v[9] = (v[5] * h_ion * (W_ion * W_ion)) - D_ion - Etot_ion;

    if (v[5] < pow(10.0, -8.0)) {
      v[5] = pow(10.0, -8.0);
    }
    if (v[9] < pow(10.0, -8.0)) {
      v[9] = pow(10.0, -8.0);
    }

    v[10] = Dx; v[11] = Dy; v[12] = Dz;
    v[13] = Bx; v[14] = By; v[15] = Bz;

    v[16] = phi;
    v[17] = psi;

    v[18] = lapse;
    v[19] = shift_x;
    v[20] = shift_y;
    v[21] = shift_z;

    v[22] = spatial_metric[0][0]; v[23] = spatial_metric[0][1]; v[24] = spatial_metric[0][2];
    v[25] = spatial_metric[1][0]; v[26] = spatial_metric[1][1]; v[27] = spatial_metric[1][2];
    v[28] = spatial_metric[2][0]; v[29] = spatial_metric[2][1]; v[30] = spatial_metric[2][2];

    v[31] = extrinsic_curvature[0][0]; v[32] = extrinsic_curvature[0][1]; v[33] = extrinsic_curvature[0][2];
    v[34] = extrinsic_curvature[1][0]; v[35] = extrinsic_curvature[1][1]; v[36] = extrinsic_curvature[1][2];
    v[37] = extrinsic_curvature[2][0]; v[38] = extrinsic_curvature[2][1]; v[39] = extrinsic_curvature[2][2];

    v[40] = 1.0;

    v[41] = lapse_der[0];
    v[42] = lapse_der[1];
    v[43] = lapse_der[2];

    v[44] = shift_der[0][0]; v[45] = shift_der[0][1]; v[46] = shift_der[0][2];
    v[47] = shift_der[1][0]; v[48] = shift_der[1][1]; v[49] = shift_der[1][2];
    v[50] = shift_der[2][0]; v[51] = shift_der[2][1]; v[52] = shift_der[2][2];

    v[53] = spatial_metric_der[0][0][0]; v[54] = spatial_metric_der[0][0][1]; v[55] = spatial_metric_der[0][0][2];
    v[56] = spatial_metric_der[0][1][0]; v[57] = spatial_metric_der[0][1][1]; v[58] = spatial_metric_der[0][1][2];
    v[59] = spatial_metric_der[0][2][0]; v[60] = spatial_metric_der[0][2][1]; v[61] = spatial_metric_der[0][2][2];

    v[62] = spatial_metric_der[1][0][0]; v[63] = spatial_metric_der[1][0][1]; v[64] = spatial_metric_der[1][0][2];
    v[65] = spatial_metric_der[1][1][0]; v[66] = spatial_metric_der[1][1][1]; v[67] = spatial_metric_der[1][1][2];
    v[68] = spatial_metric_der[1][2][0]; v[69] = spatial_metric_der[1][2][1]; v[70] = spatial_metric_der[1][2][2];

    v[71] = spatial_metric_der[2][0][0]; v[72] = spatial_metric_der[2][0][1]; v[73] = spatial_metric_der[2][0][2];
    v[74] = spatial_metric_der[2][1][0]; v[75] = spatial_metric_der[2][1][1]; v[76] = spatial_metric_der[2][1][2];
    v[77] = spatial_metric_der[2][2][0]; v[78] = spatial_metric_der[2][2][1]; v[79] = spatial_metric_der[2][2][2];

    v[80] = evol_param;
    v[81] = x;
    v[82] = y;
    v[83] = z;
  }
  else {
    for (int i = 0; i < 84; i++) {
      v[i] = 0.0;
    }
    
    v[40] = -1.0;
  }
}

void 
gkyl_gr_twofluid_inv_spatial_metric(const double q[84], double ***inv_spatial_metric)
{
  double spatial_metric[3][3];
  spatial_metric[0][0] = q[22]; spatial_metric[0][1] = q[23]; spatial_metric[0][2] = q[24];
  spatial_metric[1][0] = q[25]; spatial_metric[1][1] = q[26]; spatial_metric[1][2] = q[27];
  spatial_metric[2][0] = q[28]; spatial_metric[2][1] = q[29]; spatial_metric[2][2] = q[30];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += spatial_metric[i][i];
  }

  double spatial_metric_sq[3][3];
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

  double euclidean_metric[3][3];
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
}

void
gkyl_gr_twofluid_stress_energy_tensor_elc(double gas_gamma_elc, double gas_gamma_ion, const double q[84], double ***stress_energy_elc)
{
  double v[84] = { 0.0 };
  gkyl_gr_twofluid_prim_vars(gas_gamma_elc, gas_gamma_ion, q, v);
  double rho_elc = v[0];
  double vx_elc = v[1];
  double vy_elc = v[2];
  double vz_elc = v[3];
  double p_elc = v[4];

  double lapse = v[18];
  double shift_x = v[19];
  double shift_y = v[20];
  double shift_z = v[21];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[22]; spatial_metric[0][1] = v[23]; spatial_metric[0][2] = v[24];
  spatial_metric[1][0] = v[25]; spatial_metric[1][1] = v[26]; spatial_metric[1][2] = v[27];
  spatial_metric[2][0] = v[28]; spatial_metric[2][1] = v[29]; spatial_metric[2][2] = v[30];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_twofluid_inv_spatial_metric(q, &inv_spatial_metric);

  bool in_excision_region = false;
  if (v[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel_elc[3];
    double v_sq_elc = 0.0;
    vel_elc[0] = vx_elc; vel_elc[1] = vy_elc; vel_elc[2] = vz_elc;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_elc += spatial_metric[i][j] * vel_elc[i] * vel_elc[j];
      }
    }

    double W_elc = 1.0 / sqrt(1.0 - v_sq_elc);
    if (v_sq_elc > 1.0 - pow(10.0, -8.0)) {
      W_elc = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double h_elc = 1.0 + ((p_elc / rho_elc) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));

    double spacetime_vel_elc[4];
    spacetime_vel_elc[0] = W_elc / lapse;
    spacetime_vel_elc[1] = (W_elc * vx_elc) - (shift_x * (W_elc / lapse));
    spacetime_vel_elc[2] = (W_elc * vy_elc) - (shift_y * (W_elc / lapse));
    spacetime_vel_elc[3] = (W_elc * vz_elc) - (shift_z * (W_elc / lapse));

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double inv_spacetime_metric[4][4];
    inv_spacetime_metric[0][0] = - (1.0 / (lapse * lapse));
    for (int i = 0; i < 3; i++) {
      inv_spacetime_metric[0][i] = (1.0 / (lapse * lapse)) * shift[i];
      inv_spacetime_metric[i][0] = (1.0 / (lapse * lapse)) * shift[i];
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spacetime_metric[i][j] = inv_spatial_metric[i][j] - ((1.0 / (lapse * lapse)) * shift[i] * shift[j]);
      }
    }

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        (*stress_energy_elc)[i][j] = (rho_elc * h_elc * spacetime_vel_elc[i] * spacetime_vel_elc[j]) + (p_elc * inv_spacetime_metric[i][j]);
      }
    }
  }
  else {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        (*stress_energy_elc)[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(inv_spatial_metric);
}

void
gkyl_gr_twofluid_stress_energy_tensor_ion(double gas_gamma_elc, double gas_gamma_ion, const double q[84], double ***stress_energy_ion)
{
  double v[84] = { 0.0 };
  gkyl_gr_twofluid_prim_vars(gas_gamma_elc, gas_gamma_ion, q, v);
  double rho_ion = v[5];
  double vx_ion = v[6];
  double vy_ion = v[7];
  double vz_ion = v[8];
  double p_ion = v[9];

  double lapse = v[18];
  double shift_x = v[19];
  double shift_y = v[20];
  double shift_z = v[21];

  double spatial_metric[3][3];
  spatial_metric[0][0] = v[22]; spatial_metric[0][1] = v[23]; spatial_metric[0][2] = v[24];
  spatial_metric[1][0] = v[25]; spatial_metric[1][1] = v[26]; spatial_metric[1][2] = v[27];
  spatial_metric[2][0] = v[28]; spatial_metric[2][1] = v[29]; spatial_metric[2][2] = v[30];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_gr_twofluid_inv_spatial_metric(q, &inv_spatial_metric);

  bool in_excision_region = false;
  if (v[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double vel_ion[3];
    double v_sq_ion = 0.0;
    vel_ion[0] = vx_ion; vel_ion[1] = vy_ion; vel_ion[2] = vz_ion;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq_ion += spatial_metric[i][j] * vel_ion[i] * vel_ion[j];
      }
    }

    double W_ion = 1.0 / sqrt(1.0 - v_sq_ion);
    if (v_sq_ion > 1.0 - pow(10.0, -8.0)) {
      W_ion = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double h_ion = 1.0 + ((p_ion / rho_ion) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));

    double spacetime_vel_ion[4];
    spacetime_vel_ion[0] = W_ion / lapse;
    spacetime_vel_ion[1] = (W_ion * vx_ion) - (shift_x * (W_ion / lapse));
    spacetime_vel_ion[2] = (W_ion * vy_ion) - (shift_y * (W_ion / lapse));
    spacetime_vel_ion[3] = (W_ion * vz_ion) - (shift_z * (W_ion / lapse));

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double inv_spacetime_metric[4][4];
    inv_spacetime_metric[0][0] = - (1.0 / (lapse * lapse));
    for (int i = 0; i < 3; i++) {
      inv_spacetime_metric[0][i] = (1.0 / (lapse * lapse)) * shift[i];
      inv_spacetime_metric[i][0] = (1.0 / (lapse * lapse)) * shift[i];
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spacetime_metric[i][j] = inv_spatial_metric[i][j] - ((1.0 / (lapse * lapse)) * shift[i] * shift[j]);
      }
    }

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        (*stress_energy_ion)[i][j] = (rho_ion * h_ion * spacetime_vel_ion[i] * spacetime_vel_ion[j]) + (p_ion * inv_spacetime_metric[i][j]);
      }
    }
  }
  else {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        (*stress_energy_ion)[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(inv_spatial_metric);
}