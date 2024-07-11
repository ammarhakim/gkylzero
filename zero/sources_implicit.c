#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_fv_proj.h>
#include <gkyl_sources_implicit.h>
#include <gkyl_mat.h>

void
implicit_source_coupling_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], const double* p_rhs_s[GKYL_MAX_SPECIES], double* em, const double* app_current,
  const double* ext_em, const double* nT_source_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;
  double ke_old[GKYL_MAX_SPECIES];
  double p_inp[6], p_source[6], p_tensor[GKYL_MAX_SPECIES][6];

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

      p_inp[0] = p11 - ((mom_x * mom_x) / rho);
      p_inp[1] = p12 - ((mom_x * mom_y) / rho);
      p_inp[2] = p13 - ((mom_x * mom_z) / rho);
      p_inp[3] = p22 - ((mom_y * mom_y) / rho);
      p_inp[4] = p23 - ((mom_y * mom_z) / rho);
      p_inp[5] = p33 - ((mom_z * mom_z) / rho);

      double p11_rhs = p_rhs[4], p12_rhs = p_rhs[5], p13_rhs = p_rhs[6];
      double p22_rhs = p_rhs[7], p23_rhs = p_rhs[8], p33_rhs = p_rhs[9];

      p_source[0] = p_inp[0] + (0.5 * dt * p11_rhs);
      p_source[1] = p_inp[1] + (0.5 * dt * p12_rhs);
      p_source[2] = p_inp[2] + (0.5 * dt * p13_rhs);
      p_source[3] = p_inp[3] + (0.5 * dt * p22_rhs);
      p_source[4] = p_inp[4] + (0.5 * dt * p23_rhs);
      p_source[5] = p_inp[5] + (0.5 * dt * p33_rhs);

      double p = (1.0 / 3.0) * (p_inp[0] + p_inp[3] + p_inp[5]);
      double v_th = sqrt(p / rho);
      double nu = v_th * k0;
      double exp_nu = exp(nu * dt);

      if (mom_em->is_charged_species) {
        // TODO: Add pressure tensor rotation here.
      }

      p_tensor[i][0] = ((p_tensor[i][0] - p) / exp_nu) + p;
      p_tensor[i][1] = p_tensor[i][1] / exp_nu;
      p_tensor[i][2] = p_tensor[i][2] / exp_nu;
      p_tensor[i][3] = ((p_tensor[i][3] - p) / exp_nu) + p;
      p_tensor[i][4] = p_tensor[i][4] / exp_nu;
      p_tensor[i][5] = ((p_tensor[i][5] - p) / exp_nu) + p;
    }
  }

  if (mom_em->is_charged_species) {
    // TODO: Add implicit source updater for electromagnetic sources here.
  }
  else {
    // TODO: Add implicit source updater for neutrals here.
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

      f[4] = ((mom_x * mom_x) / rho) + p_tensor[i][0];
      f[5] = ((mom_x * mom_y) / rho) + p_tensor[i][1];
      f[6] = ((mom_x * mom_z) / rho) + p_tensor[i][2];
      f[7] = ((mom_y * mom_y) / rho) + p_tensor[i][3];
      f[8] = ((mom_y * mom_z) / rho) + p_tensor[i][4];
      f[9] = ((mom_z * mom_z) / rho) + p_tensor[i][5];
    }
  }

  if (mom_em->has_collision) {
    // TODO: Add implicit source updater for collisions.
  }

  if (mom_em->has_nT_sources) {
    // TODO: Add implicit source updater for number density and temperature.
  }
}