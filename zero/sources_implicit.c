#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_fv_proj.h>
#include <gkyl_sources_implicit.h>
#include <gkyl_sources_explicit.h>
#include <gkyl_source_utils.h>
#include <gkyl_mat.h>

void
implicit_em_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double *app_accel_s[GKYL_MAX_SPECIES], double* em, const double* app_current, const double* ext_em)
{
  int nfluids = mom_em->nfluids;
  double epsilon0 = mom_em->epsilon0;

  double Bx = em[3] + ext_em[3];
  double By = em[4] + ext_em[4];
  double Bz = em[5] + ext_em[5];
  double B_mag = sqrt((Bx * Bx) + (By * By) + (Bz * Bz));
  
  double bx = 0.0, by = 0.0, bz = 0.0;
  if (B_mag > 0.0) {
    bx = Bx / B_mag;
    by = By / B_mag;
    bz = Bz / B_mag;
  }

  double q_over_m[GKYL_MAX_SPECIES];
  double wc_dt[GKYL_MAX_SPECIES];
  double wp_dt_sq[GKYL_MAX_SPECIES];
  double J_old[GKYL_MAX_SPECIES][3];
  double J[GKYL_MAX_SPECIES][3];

  double scale_fact_E = 1.0;
  double scale_fact_curr = 1.0;
  double t_ramp_E = mom_em->t_ramp_E;
  double t_ramp_curr = mom_em->t_ramp_curr;

  double w0_sq = 0.0;
  double gam_sq = 0.0;
  double delta = 0.0;
  double Kx = 0.0, Ky = 0.0, Kz = 0.0;

  for (int i = 0; i < nfluids; i++) {
    double q = mom_em->param[i].charge;
    double m = mom_em->param[i].mass;
    q_over_m[i] = q / m;

    const double *f = fluid_s[i];
    const double *app_accel = app_accel_s[i];

    double rho = f[0];
    double mom_x = f[1], mom_y = f[2], mom_z = f[3];

    J_old[i][0] = mom_x * q_over_m[i];
    J_old[i][1] = mom_y * q_over_m[i];
    J_old[i][2] = mom_z * q_over_m[i];

    if (mom_em->ramp_app_E) {
      scale_fact_E = fmin(1.0, t_curr / t_ramp_E);
    }

    J[i][0] = J_old[i][0] + (0.5 * dt * q_over_m[i] * rho * ((q_over_m[i] * ext_em[0] * scale_fact_E) + app_accel[0]));
    J[i][1] = J_old[i][1] + (0.5 * dt * q_over_m[i] * rho * ((q_over_m[i] * ext_em[1] * scale_fact_E) + app_accel[1]));
    J[i][2] = J_old[i][2] + (0.5 * dt * q_over_m[i] * rho * ((q_over_m[i] * ext_em[2] * scale_fact_E) + app_accel[2]));

    wc_dt[i] = q_over_m[i] * B_mag * dt;
    wp_dt_sq[i] = (rho * (q_over_m[i] * q_over_m[i]) * (dt * dt)) / epsilon0;

    double denom = 1.0 + ((wc_dt[i] * wc_dt[i]) / 4.0);
    w0_sq += wp_dt_sq[i] / denom;
    gam_sq += (wp_dt_sq[i] * (wc_dt[i] * wc_dt[i])) / denom;
    delta += (wp_dt_sq[i] * wc_dt[i]) / denom;

    Kx -= (dt / denom) * (J[i][0] + (((wc_dt[i] * wc_dt[i]) / 4.0) * bx * ((bx * J[i][0]) + (by * J[i][1]) + (bz * J[i][2]))) -
      ((wc_dt[i] / 2.0) * ((by * J[i][2]) - (bz * J[i][1]))));
    Ky -= (dt / denom) * (J[i][1] + (((wc_dt[i] * wc_dt[i]) / 4.0) * by * ((bx * J[i][0]) + (by * J[i][1]) + (bz * J[i][2]))) -
      ((wc_dt[i] / 2.0) * ((bz * J[i][0]) - (bx * J[i][2]))));
    Kz -= (dt / denom) * (J[i][2] + (((wc_dt[i] * wc_dt[i]) / 4.0) * bz * ((bx * J[i][0]) + (by * J[i][1]) + (bz * J[i][2]))) -
      ((wc_dt[i] / 2.0) * ((bx * J[i][1]) - (by * J[i][0]))));
  }

  double Delta_sq = (delta * delta) / (1.0 + (w0_sq / 4.0));

  double Fx_old = em[0] * epsilon0;
  double Fy_old = em[1] * epsilon0;
  double Fz_old = em[2] * epsilon0;

  if (mom_em->ramp_app_curr) {
    scale_fact_curr = fmin(1.0, t_curr / t_ramp_curr);
  }

  double Fx = Fx_old - (0.5 * dt * app_current[0] * scale_fact_curr);
  double Fy = Fy_old - (0.5 * dt * app_current[1] * scale_fact_curr);
  double Fz = Fz_old - (0.5 * dt * app_current[2] * scale_fact_curr);

  double Fx_K = Fx + (0.5 * Kx);
  double Fy_K = Fy + (0.5 * Ky);
  double Fz_K = Fz + (0.5 * Kz);

  double Fx_bar = (1.0 / (1.0 + (w0_sq / 4.0) + (Delta_sq / 64.0))) * (Fx_K + ((((Delta_sq / 64.0) - (gam_sq / 16.0)) /
    (1.0 + (w0_sq / 4.0) + (gam_sq / 16.0))) * bx * ((bx * Fx_K) + (by * Fy_K) + (bz * Fz_K))) + (((delta / 8.0) /
    (1.0 + (w0_sq / 4.0))) * ((by * Fz_K) - (bz * Fy_K))));
  double Fy_bar = (1.0 / (1.0 + (w0_sq / 4.0) + (Delta_sq / 64.0))) * (Fy_K + ((((Delta_sq / 64.0) - (gam_sq / 16.0)) /
    (1.0 + (w0_sq / 4.0) + (gam_sq / 16.0))) * by * ((bx * Fx_K) + (by * Fy_K) + (bz * Fz_K))) + (((delta / 8.0) /
    (1.0 + (w0_sq / 4.0))) * ((bz * Fx_K) - (bx * Fz_K))));
  double Fz_bar = (1.0 / (1.0 + (w0_sq / 4.0) + (Delta_sq / 64.0))) * (Fz_K + ((((Delta_sq / 64.0) - (gam_sq / 16.0)) /
    (1.0 + (w0_sq / 4.0) + (gam_sq / 16.0))) * bz * ((bx * Fx_K) + (by * Fy_K) + (bz * Fz_K))) + (((delta / 8.0) /
    (1.0 + (w0_sq / 4.0))) * ((bx * Fy_K) - (by * Fx_K))));
  
  em[0] = ((2.0 * Fx_bar) - Fx_old) / epsilon0;
  em[1] = ((2.0 * Fy_bar) - Fy_old) / epsilon0;
  em[2] = ((2.0 * Fz_bar) - Fz_old) / epsilon0;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double Jx_star = J[i][0] + (Fx_bar * ((wp_dt_sq[i] / dt) / 2.0));
    double Jy_star = J[i][1] + (Fy_bar * ((wp_dt_sq[i] / dt) / 2.0));
    double Jz_star = J[i][2] + (Fz_bar * ((wp_dt_sq[i] / dt) / 2.0));

    double Jx_new = ((2.0 * (Jx_star + (((wc_dt[i] * wc_dt[i]) / 4.0) * bx * ((bx * Jx_star) + (by * Jy_star) + (bz * Jz_star))) -
      ((wc_dt[i] / 2.0) * ((by * Jz_star) - (bz * Jy_star))))) / (1.0 + ((wc_dt[i] * wc_dt[i]) / 4.0))) - J_old[i][0];
    double Jy_new = ((2.0 * (Jy_star + (((wc_dt[i] * wc_dt[i]) / 4.0) * by * ((bx * Jx_star) + (by * Jy_star) + (bz * Jz_star))) -
      ((wc_dt[i] / 2.0) * ((bz * Jx_star) - (bx * Jz_star))))) / (1.0 + ((wc_dt[i] * wc_dt[i]) / 4.0))) - J_old[i][1];
    double Jz_new = ((2.0 * (Jz_star + (((wc_dt[i] * wc_dt[i]) / 4.0) * bz * ((bx * Jx_star) + (by * Jy_star) + (bz * Jz_star))) -
      ((wc_dt[i] / 2.0) * ((bx * Jy_star) - (by * Jx_star))))) / (1.0 + ((wc_dt[i] * wc_dt[i]) / 4.0))) - J_old[i][2];
    
    f[1] = Jx_new / q_over_m[i];
    f[2] = Jy_new / q_over_m[i];
    f[3] = Jz_new / q_over_m[i];
  }
}

void
implicit_neut_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];
    const double *app_accel = app_accel_s[i];

    double rho = f[0];

    f[1] += dt * rho * app_accel[0];
    f[2] += dt * rho * app_accel[1];
    f[3] += dt * rho * app_accel[2];
  }
}

void
implicit_collision_source_update(const gkyl_moment_em_coupling* mom_em, double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;
  double nu_base[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];
  for (int i  = 0; i < nfluids; i++) {
    for (int j = 0; j < nfluids; j++) {
      nu_base[i][j] = (mom_em->nu_base)[i][j];
    }
  }

  double nu[nfluids * nfluids];
  for (int i = 0; i < nfluids; i++) {
    double *nu_i = nu + (nfluids * i);

    for (int j = 0; j < nfluids; j++) {
      double rho = fluid_s[j][0];

      nu_i[j] = nu_base[i][j] * rho;
    }
  }

  double lhs[nfluids][nfluids];
  double rhs[nfluids][3];
  for (int i = 0; i < nfluids; i++) {
    for (int j =0 ; j < nfluids; j++) {
      lhs[i][j] = 0.0;
    }

    rhs[i][0] = 0.0; rhs[i][1] = 0.0; rhs[i][2] = 0.0;
  }

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double rho = f[0];
    double mom_x = f[1], mom_y = f[2], mom_z = f[3];

    rhs[i][0] = mom_x / rho;
    rhs[i][1] = mom_y / rho;
    rhs[i][2] = mom_z / rho;

    lhs[i][i] = 1.0;
    double* nu_i = nu + (nfluids * i);

    for (int j = 0; j < nfluids; j++) {
      if (i == j) {
        lhs[i][i] += 0.5 * dt * nu_i[i];
      }
      else {
        double dt_nu_ij = 0.5 * dt * nu_i[j];
        lhs[i][i] += dt_nu_ij;
        lhs[i][j] -= dt_nu_ij;
      }
    }
  }

  struct gkyl_mat *lhs_mat = gkyl_mat_new(nfluids, nfluids, 0.0);
  struct gkyl_mat *rhs_mat = gkyl_mat_new(nfluids, 3, 0.0);

  for (int i = 0; i < nfluids; i++) {
    for (int j = 0; j < nfluids; j++) {
      gkyl_mat_set(lhs_mat, i, j, lhs[i][j]);
    }

    gkyl_mat_set(rhs_mat, i, 0, rhs[i][0]);
    gkyl_mat_set(rhs_mat, i, 1, rhs[i][1]);
    gkyl_mat_set(rhs_mat, i, 2, rhs[i][2]);
  }

  gkyl_mem_buff sol_buff = gkyl_mem_buff_new(sizeof(long[nfluids]));
  bool status = gkyl_mat_linsolve_lu(lhs_mat, rhs_mat, gkyl_mem_buff_data(sol_buff));

  for (int i = 0; i < nfluids; i++) {
    for (int j = 0; j < 3; j++) {
      rhs[i][j] = gkyl_mat_get(rhs_mat, i, j);
    }
  }

  double rhs_T[nfluids][1];
  for (int i = 0; i < nfluids; i++) {
    rhs_T[i][0] = 0.0;

    for (int j = 0; j < nfluids; j++) {
      lhs[i][j] = 0.0;
    }
  }

  double T[nfluids];
  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];
    double m = mom_em->param[i].mass;

    double rho = f[0];
    double mom_x = f[1], mom_y = f[2], mom_z = f[3];
    double E = f[4];
    double internal_energy = E - (0.5 * ((mom_x * mom_x) + (mom_y * mom_y) + (mom_z * mom_z)) / rho);

    T[i] = (internal_energy / rho) * m;
    rhs_T[i][0] = T[i];
    lhs[i][i] = 1.0;

    double *nu_i = nu + (nfluids * i);

    for (int j = 0; j < nfluids; j++) {
      if (i != j) {
        double m_j = mom_em->param[j].mass;
        
        double du_sq = ((rhs[i][0] - rhs[j][0]) * (rhs[i][0] - rhs[j][0])) + ((rhs[i][1] - rhs[j][1]) * (rhs[i][1] - rhs[j][1])) +
          ((rhs[i][2] - rhs[j][2]) * (rhs[i][2] - rhs[j][2]));
        double coeff_ij = (dt * nu_i[j] * m) / (m + m_j);

        rhs_T[i][0] += 0.5 * coeff_ij * m_j * du_sq;
        lhs[i][i] += coeff_ij;
        lhs[i][j] -= coeff_ij;
      }
    }
  }

  struct gkyl_mat *lhs_T_mat = gkyl_mat_new(nfluids, nfluids, 0.0);
  struct gkyl_mat *rhs_T_mat = gkyl_mat_new(nfluids, 1, 0.0);

  for (int i = 0; i < nfluids; i++) {
    for (int j = 0; j < nfluids; j++) {
      gkyl_mat_set(lhs_T_mat, i, j, lhs[i][j]);
    }

    gkyl_mat_set(rhs_T_mat, i, 0, rhs_T[i][0]);
  }

  status = gkyl_mat_linsolve_lu(lhs_T_mat, rhs_T_mat, gkyl_mem_buff_data(sol_buff));

  for (int i = 0; i < nfluids; i++) {
    rhs_T[i][0] = gkyl_mat_get(rhs_T_mat, i, 0);
  }

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];
    double m = mom_em->param[i].mass;

    double rho = f[0];

    f[4] = ((2.0 * rhs_T[i][0] - T[i]) * rho) / m;
  }

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double rho = f[0];
    double mom_x = f[1], mom_y = f[2], mom_z = f[3];
    double E = f[4];

    f[1] = 2.0 * rho * rhs[i][0] - mom_x;
    f[2] = 2.0 * rho * rhs[i][1] - mom_y;
    f[3] = 2.0 * rho * rhs[i][2] - mom_z;

    mom_x = f[1], mom_y = f[2], mom_z = f[3];
    f[4] = E + (0.5 * ((mom_x * mom_x) + (mom_y * mom_y) + (mom_z * mom_z)) / rho);
  }

  gkyl_mem_buff_release(sol_buff);
  gkyl_mat_release(lhs_mat);
  gkyl_mat_release(rhs_mat);
  gkyl_mat_release(rhs_T_mat);
}

void
implicit_source_coupling_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], const double* p_rhs_s[GKYL_MAX_SPECIES], double* em, const double* app_current,
  const double* ext_em, const double* nT_sources_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;
  double ke_old[GKYL_MAX_SPECIES];
  double p_tensor_old[6], p_tensor_rhs[6], p_tensor_new[GKYL_MAX_SPECIES][6];

  for (int i = 0; i < nfluids; i++) {
    const double *f = fluid_s[i];
    const double *p_rhs = p_rhs_s[i];

    double q = mom_em->param[i].charge;
    double m = mom_em->param[i].mass;
    double k0 = mom_em->param[i].k0;

    if (mom_em->param[i].type == GKYL_EQN_EULER) {
      double rho = f[0];
      double mom_x = f[1], mom_y = f[2], mom_z = f[3];

      ke_old[i] = 0.5 * (((mom_x * mom_x) + (mom_y * mom_y) + (mom_z * mom_z)) / rho);
    }
    else if (mom_em->param[i].type == GKYL_EQN_TEN_MOMENT) {
      double q_over_m = q / m;

      double rho = f[0];
      double mom_x = f[1], mom_y = f[2], mom_z = f[3];
      double p11 = f[4], p12 = f[5], p13 = f[6];
      double p22 = f[7], p23 = f[8], p33 = f[9];

      p_tensor_old[0] = p11 - ((mom_x * mom_x) / rho);
      p_tensor_old[1] = p12 - ((mom_x * mom_y) / rho);
      p_tensor_old[2] = p13 - ((mom_x * mom_z) / rho);
      p_tensor_old[3] = p22 - ((mom_y * mom_y) / rho);
      p_tensor_old[4] = p23 - ((mom_y * mom_z) / rho);
      p_tensor_old[5] = p33 - ((mom_z * mom_z) / rho);

      double p11_rhs = p_rhs[4], p12_rhs = p_rhs[5], p13_rhs = p_rhs[6];
      double p22_rhs = p_rhs[7], p23_rhs = p_rhs[8], p33_rhs = p_rhs[9];

      p_tensor_rhs[0] = p_tensor_old[0] + (0.5 * dt * p11_rhs);
      p_tensor_rhs[1] = p_tensor_old[1] + (0.5 * dt * p12_rhs);
      p_tensor_rhs[2] = p_tensor_old[2] + (0.5 * dt * p13_rhs);
      p_tensor_rhs[3] = p_tensor_old[3] + (0.5 * dt * p22_rhs);
      p_tensor_rhs[4] = p_tensor_old[4] + (0.5 * dt * p23_rhs);
      p_tensor_rhs[5] = p_tensor_old[5] + (0.5 * dt * p33_rhs);

      double p = (1.0 / 3.0) * (p_tensor_old[0] + p_tensor_old[3] + p_tensor_old[5]);
      double v_th = sqrt(p / rho);
      double nu = v_th * k0;
      double exp_nu = exp(nu * dt);

      if (mom_em->is_charged_species) {
        pressure_tensor_rotate(q_over_m, dt, em, ext_em, p_tensor_old, p_tensor_rhs, p_tensor_new[i]);
      }

      p_tensor_new[i][0] = ((p_tensor_new[i][0] - p) / exp_nu) + p;
      p_tensor_new[i][1] = p_tensor_new[i][1] / exp_nu;
      p_tensor_new[i][2] = p_tensor_new[i][2] / exp_nu;
      p_tensor_new[i][3] = ((p_tensor_new[i][3] - p) / exp_nu) + p;
      p_tensor_new[i][4] = p_tensor_new[i][4] / exp_nu;
      p_tensor_new[i][5] = ((p_tensor_new[i][5] - p) / exp_nu) + p;
    }
  }

  if (mom_em->is_charged_species) {
    implicit_em_source_update(mom_em, t_curr, dt, fluid_s, app_accel_s, em, app_current, ext_em);
  }
  else {
    implicit_neut_source_update(mom_em, t_curr, dt, fluid_s, app_accel_s);
  }

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    if (mom_em->param[i].type == GKYL_EQN_EULER) {
      double rho = f[0];
      double mom_x = f[1], mom_y = f[2], mom_z = f[3];
      
      f[4] += (0.5 * ((mom_x * mom_x) + (mom_y * mom_y) + (mom_z * mom_z)) / rho) - ke_old[i];
    }
    else {
      double rho = f[0];
      double mom_x = f[1], mom_y = f[2], mom_z = f[3];

      f[4] = ((mom_x * mom_x) / rho) + p_tensor_new[i][0];
      f[5] = ((mom_x * mom_y) / rho) + p_tensor_new[i][1];
      f[6] = ((mom_x * mom_z) / rho) + p_tensor_new[i][2];
      f[7] = ((mom_y * mom_y) / rho) + p_tensor_new[i][3];
      f[8] = ((mom_y * mom_z) / rho) + p_tensor_new[i][4];
      f[9] = ((mom_z * mom_z) / rho) + p_tensor_new[i][5];
    }
  }

  if (mom_em->has_collision) {
    implicit_collision_source_update(mom_em, dt, fluid_s);
  }

  // These terms are handled by their own specialized explicit forcing solver(s). To be revisited...
  if (mom_em->has_nT_sources) {
    explicit_nT_source_update(mom_em, dt, fluid_s, nT_sources_s);
  }
  if (mom_em->has_frictional_sources) {
    int subcycling = 1000;

    for (int i = 0; i < subcycling; i++) {
      explicit_frictional_source_update(mom_em, t_curr, (1.0 / subcycling) * dt, fluid_s);
    }
  }
}