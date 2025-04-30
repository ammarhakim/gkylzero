#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_fv_proj.h>
#include <gkyl_sources_explicit.h>
#include <gkyl_moment_em_coupling_priv.h>
#include <gkyl_mat.h>

void
explicit_nT_source_update_euler(const double mass, const double dt, double* fluid_old, double* fluid_new, const double* nT_sources)
{
  double rho_old = fluid_old[0];
  double n_old = rho_old / mass;

  double u = fluid_old[1] / rho_old;
  double v = fluid_old[2] / rho_old;
  double w = fluid_old[3] / rho_old;
  double v_sq = (u * u) + (v * v) + (w * w);
  double TT_old = (fluid_old[4] - (0.5 * rho_old * v_sq)) / n_old;

  double n_new = n_old + (dt * nT_sources[0]);
  double TT_new = TT_old + (dt * nT_sources[1]);

  double rho_new = n_new * mass;
  fluid_new[0] = rho_new;
  fluid_new[1] = rho_new * u;
  fluid_new[2] = rho_new * v;
  fluid_new[3] = rho_new * w;
  fluid_new[4] = (n_new * TT_new) + (0.5 * rho_new * v_sq);
}

void
explicit_nT_source_update(const gkyl_moment_em_coupling* mom_em, const double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* nT_sources_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];
    const double *nT_sources = nT_sources_s[i];

    double mass = mom_em->param[i].mass;

    explicit_nT_source_update_euler(mass, dt, f, f, nT_sources);
  }
}

void
explicit_frictional_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double Z, const double T_elc, const double Lambda_ee,
  double t_curr, const double dt, double* f_elc_old, double* f_ion_old, double* f_elc_new, double* f_ion_new)
{
  int nfluids = mom_em->nfluids;
  double pi = M_PI;

  if (nfluids == 2) {
    double mass_elc = mom_em->param[0].mass;
    double epsilon0 = mom_em->epsilon0;
    
    double rho_elc = f_elc_old[0];
    double rho_ion = f_ion_old[0];

    double u_elc = f_elc_old[1], v_elc = f_elc_old[2], w_elc = f_elc_old[3];
    double u_ion = f_ion_old[1], v_ion = f_ion_old[2], w_ion = f_ion_old[3];

    double n_elc = rho_elc / mass_elc;

    double tau_ei = (1.0 / Z) * ((3.0 * sqrt(mass_elc) * ((4.0 * pi * epsilon0) * (4.0 * pi * epsilon0)) * pow(T_elc, 3.0 / 2.0)) /
      (4.0 * sqrt(2.0 * pi) * n_elc * exp(4.0) * log(Lambda_ee)));
    double alpha_par = 1.0 - (pow(Z, 2.0 / 3.0) / ((1.46 * pow(Z, 2.0 / 3.0)) - (0.33 * pow (Z, 1.0 / 3.0)) + 0.888));

    double mom_src_x = -(rho_elc / tau_ei) * (alpha_par * (u_elc - u_ion));
    double mom_src_y = -(rho_elc / tau_ei) * (alpha_par * (v_elc - v_ion));
    double mom_src_z = -(rho_elc / tau_ei) * (alpha_par * (w_elc - w_ion));

    double E_src = (mom_src_x * u_ion) + (mom_src_y * v_ion) + (mom_src_z * w_ion);

    f_elc_new[1] = f_elc_old[1] + (dt * mom_src_x);
    f_elc_new[2] = f_elc_old[2] + (dt * mom_src_y);
    f_elc_new[3] = f_elc_old[3] + (dt * mom_src_z);

    f_ion_new[1] = f_ion_old[1] - (dt * mom_src_x);
    f_ion_new[2] = f_ion_old[2] - (dt * mom_src_y);
    f_ion_new[3] = f_ion_old[3] - (dt * mom_src_z);

    if (mom_em->param[0].type == GKYL_EQN_EULER) {
      f_elc_new[4] = f_elc_old[4] + (dt * E_src);
    }
    if (mom_em->param[1].type == GKYL_EQN_EULER) {
      f_ion_new[4] = f_ion_old[4] - (dt * E_src);
    }

    f_elc_new[0] = f_elc_old[0];
    f_ion_new[0] = f_ion_old[0];
  }
}

void
explicit_frictional_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  if (nfluids == 2) {
    double *f_elc = fluid_s[0];
    double *f_ion = fluid_s[1];

    double Z = mom_em->friction_Z;
    double T_elc = mom_em->friction_T_elc;
    double Lambda_ee = mom_em->friction_Lambda_ee;

    int elc_num_equations = 0;
    int ion_num_equations = 0;

    if (mom_em->param[0].type == GKYL_EQN_EULER) {
      elc_num_equations = 5;
    }
    else if (mom_em->param[0].type == GKYL_EQN_ISO_EULER) {
      elc_num_equations = 4;
    }

    // Right now, we allocate these arrays on the stack with size 5 (current maximum supported number of equations).
    double f_elc_new[5];
    double f_elc_stage1[5];
    double f_elc_stage2[5];
    double f_elc_old[5];

    if (mom_em->param[1].type == GKYL_EQN_EULER) {
      ion_num_equations = 5;
    }
    else if (mom_em->param[1].type == GKYL_EQN_ISO_EULER) {
      ion_num_equations = 4;
    }

    // Right now, we allocate these arrays on the stack with size 5 (current maximum supported number of equations).
    double f_ion_new[5];
    double f_ion_stage1[5];
    double f_ion_stage2[5];
    double f_ion_old[5];

    for (int i = 0; i < elc_num_equations; i++) {
      f_elc_old[i] = f_elc[i];
    }
    for (int i = 0; i < ion_num_equations; i++) {
      f_ion_old[i] = f_ion[i];
    }

    explicit_frictional_source_update_euler(mom_em, Z, T_elc, Lambda_ee, t_curr, dt, f_elc_old, f_ion_old, f_elc_new, f_ion_new);
    for (int i = 0; i < elc_num_equations; i++) {
      f_elc_stage1[i] = f_elc_new[i];
    }
    for (int i = 0; i < ion_num_equations; i++) {
      f_ion_stage1[i] = f_ion_new[i];
    }

    explicit_frictional_source_update_euler(mom_em, Z, T_elc, Lambda_ee, t_curr + dt, dt, f_elc_stage1, f_ion_stage1, f_elc_new, f_ion_new);
    for (int i = 0; i < elc_num_equations; i++) {
      f_elc_stage2[i] = (0.75 * f_elc_old[i]) + (0.25 * f_elc_new[i]);
    }
    for (int i = 0; i < ion_num_equations; i++) {
      f_ion_stage2[i] = (0.75 * f_ion_old[i]) + (0.25 * f_ion_new[i]);
    }

    explicit_frictional_source_update_euler(mom_em, Z, T_elc, Lambda_ee, t_curr + (0.5 * dt), dt, f_elc_stage2, f_ion_stage2, f_elc_new, f_ion_new);
    for (int i = 0; i < elc_num_equations; i++) {
      f_elc[i] = ((1.0 / 3.0) * f_elc_old[i]) + ((2.0 / 3.0) * f_elc_new[i]);
    }
    for (int i = 0; i < ion_num_equations; i++) {
      f_ion[i] = ((1.0 / 3.0) * f_ion_old[i]) + ((2.0 / 3.0) * f_ion_new[i]);
    }
  }
}

void
explicit_volume_source_5m_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, const double U0, const double R0,
  double t_curr, const double dt, double* fluid_old, double* fluid_new)
{
  double rho = fluid_old[0];
  double vx = fluid_old[1] / rho;
  double vy = fluid_old[2] / rho;
  double vz = fluid_old[3] / rho;
  double p = (gas_gamma - 1.0) * (fluid_old[4] - (0.5 * rho * ((vx * vx) + (vy * vy) + (vz * vz))));

  double a = 1.0 + ((U0 * t_curr) / R0);

  for (int i = 0; i < 5; i++) {
    fluid_new[i] = fluid_old[i];
  }

  fluid_new[0] -= dt * ((2.0 * U0) / (a * R0)) * rho;

  fluid_new[2] -= dt * (U0 / (a * R0)) * rho * vy;
  fluid_new[3] -= dt * (U0 / (a * R0)) * rho * vz;

  double p_new = p;
  p_new -= dt * (gas_gamma * ((2.0 * U0) / (a * R0)) * p);

  fluid_new[4] = (p_new / (gas_gamma - 1.0)) + (0.5 * rho * (vx * vx) + (vy * vy) + (vz * vz));
}

void
explicit_volume_source_10m_update_euler(const gkyl_moment_em_coupling* mom_em, const double U0, const double R0, double t_curr, const double dt,
  double* fluid_old, double* fluid_new)
{
  double rho = fluid_old[0];
  double vx = fluid_old[1] / rho;
  double vy = fluid_old[2] / rho;
  double vz = fluid_old[3] / rho;
  
  double Pxx = fluid_old[4], Pxy = fluid_old[5], Pxz = fluid_old[6];
  double Pyy = fluid_old[7], Pyz = fluid_old[8], Pzz = fluid_old[9];

  double a = 1.0 + ((U0 * t_curr) / R0);

  for (int i = 0; i < 10; i++) {
    fluid_new[i] = fluid_old[i];
  }

  fluid_new[0] -= dt * ((2.0 * U0) / (a * R0)) * rho;

  fluid_new[2] -= dt * (U0 / (a * R0)) * rho * vy;
  fluid_new[3] -= dt * (U0 / (a * R0)) * rho * vz;

  fluid_new[4] -= dt * ((2.0 * U0) / (a * R0)) * Pxx;
  fluid_new[5] -= dt * ((((2.0 * U0) / (a * R0)) * Pxy) + ((U0 / (a * R0)) * Pxy));
  fluid_new[6] -= dt * ((((2.0 * U0) / (a * R0)) * Pxz) + ((U0 / (a * R0)) * Pxz));
  fluid_new[7] -= dt * ((((2.0 * U0) / (a * R0)) * Pyy) + (((2.0 * U0) / (a * R0)) * Pyy));
  fluid_new[8] -= dt * ((((2.0 * U0) / (a * R0)) * Pyz) + (((2.0 * U0) / (a * R0)) * Pyz));
  fluid_new[9] -= dt * ((((2.0 * U0) / (a * R0)) * Pzz) + (((2.0 * U0) / (a * R0)) * Pzz));
}

void
explicit_volume_source_maxwell_update_euler(const gkyl_moment_em_coupling* mom_em, const double U0, const double R0, double t_curr,
  const double dt, double* em_old, double* em_new, const double* ext_em)
{
  double a = 1.0 + ((U0 * t_curr) / R0);

  double Ex = em_old[0] + ext_em[0];
  double Ey = em_old[1] + ext_em[1];
  double Ez = em_old[2] + ext_em[2];

  double Bx = em_old[3] + ext_em[3];
  double By = em_old[4] + ext_em[4];
  double Bz = em_old[5] + ext_em[5];

  for (int i = 0; i < 8; i++) {
    em_new[i] = em_old[i];
  }

  em_new[0] -= dt * (U0 / (a * R0)) * 2.0 * Ex;
  em_new[1] -= dt * (U0 / (a * R0)) * Ey;
  em_new[2] -= dt * (U0 / (a * R0)) * Ez;

  em_new[3] -= dt * (U0 / (a * R0)) * 2.0 * Bx;
  em_new[4] -= dt * (U0 / (a * R0)) * By;
  em_new[5] -= dt * (U0 / (a * R0)) * Bz;
}

void
explicit_volume_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES],
  double* em, const double* ext_em)
{
  int nfluids = mom_em->nfluids;

  double gas_gamma = mom_em->volume_gas_gamma;
  double U0 = mom_em->volume_U0;
  double R0 = mom_em->volume_R0;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    if (mom_em->param[i].type == GKYL_EQN_EULER) {
      double f_new[5], f_stage1[5], f_stage2[5], f_old[5];

      for (int j = 0; j < 5; j++) {
        f_old[j] = f[j];
      }

      explicit_volume_source_5m_update_euler(mom_em, gas_gamma, U0, R0, t_curr, dt, f_old, f_new);
      for (int j = 0; j < 5; j++) {
        f_stage1[j] = f_new[j];
      }

      explicit_volume_source_5m_update_euler(mom_em, gas_gamma, U0, R0, t_curr + dt, dt, f_stage1, f_new);
      for (int j = 0; j < 5; j++) {
        f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
      }

      explicit_volume_source_5m_update_euler(mom_em, gas_gamma, U0, R0, t_curr + (0.5 * dt), dt, f_stage2, f_new);
      for (int j = 0; j < 5; j++) {
        f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
      }
    }
    else if (mom_em->param[i].type == GKYL_EQN_TEN_MOMENT) {
      double f_new[10], f_stage1[10], f_stage2[10], f_old[10];

      for (int j = 0; j < 10; j++) {
        f_old[j] = f[j];
      }

      explicit_volume_source_10m_update_euler(mom_em, U0, R0, t_curr, dt, f_old, f_new);
      for (int j = 0; j < 10; j++) {
        f_stage1[j] = f_new[j];
      }

      explicit_volume_source_10m_update_euler(mom_em, U0, R0, t_curr + dt, dt, f_stage1, f_new);
      for (int j = 0; j < 10; j++) {
        f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
      }

      explicit_volume_source_10m_update_euler(mom_em, U0, R0, t_curr + (0.5 * dt), dt, f_stage2, f_new);
      for (int j = 0; j < 10; j++) {
        f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
      }
    }
  }

  double em_new[8], em_stage1[8], em_stage2[8], em_old[8];

  for (int i = 0; i < 8; i++) {
    em_old[i] = em[i];
  }

  explicit_volume_source_maxwell_update_euler(mom_em, U0, R0, t_curr, dt, em_old, em_new, ext_em);
  for (int i = 0; i < 8; i++) {
    em_stage1[i] = em_new[i];
  }

  explicit_volume_source_maxwell_update_euler(mom_em, U0, R0, t_curr + dt, dt, em_stage1, em_new, ext_em);
  for (int i = 0; i < 8; i++) {
    em_stage2[i] = (0.75 * em_old[i]) + (0.25 * em_new[i]);
  }

  explicit_volume_source_maxwell_update_euler(mom_em, U0, R0, t_curr + (0.5 * dt), dt, em_stage2, em_new, ext_em);
  for (int i = 0; i < 8; i++) {
    em[i] = ((1.0 / 3.0) * em_old[i]) + ((2.0 / 3.0) * em_new[i]);
  }
}

void
explicit_reactive_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, const double specific_heat_capacity,
  const double energy_of_formation, const double ignition_temperature, const double reaction_rate, double t_curr, const double dt,
  double* fluid_old, double* fluid_new)
{
  double rho = fluid_old[0];
  double vx = fluid_old[1] / rho;
  double vy = fluid_old[2] / rho;
  double vz = fluid_old[3] / rho;
  double reaction_progress = fluid_old[5] / rho;

  double specific_internal_energy = (fluid_old[4] / rho) - (0.5 * ((vx * vx) + (vy * vy) + (vz * vz))) -
    (energy_of_formation * (reaction_progress - 1.0));
  double temperature = specific_internal_energy / specific_heat_capacity;

  for (int i = 0; i < 6; i++) {
    fluid_new[i] = fluid_old[i];
  }

  if (temperature > ignition_temperature) {
    fluid_new[5] -= dt * (rho * reaction_progress * reaction_rate);
  }
}

void
explicit_reactive_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  double gas_gamma = mom_em->reactivity_gas_gamma;
  double specific_heat_capacity = mom_em->reactivity_specific_heat_capacity;
  double energy_of_formation = mom_em->reactivity_energy_of_formation;
  double ignition_temperature = mom_em->reactivity_ignition_temperature;
  double reaction_rate = mom_em->reactivity_reaction_rate;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double f_new[6], f_stage1[6], f_stage2[6], f_old[6];

    for (int j = 0; j < 6; j++) {
      f_old[j] = f[j];
    }

    explicit_reactive_source_update_euler(mom_em, gas_gamma, specific_heat_capacity, energy_of_formation, ignition_temperature, reaction_rate,
      t_curr, dt, f_old, f_new);
    for (int j = 0; j < 6; j++) {
      f_stage1[j] = f_new[j];
    }

    explicit_reactive_source_update_euler(mom_em, gas_gamma, specific_heat_capacity, energy_of_formation, ignition_temperature, reaction_rate,
      t_curr + dt, dt, f_stage1, f_new);
    for (int j = 0; j < 6; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
    }

    explicit_reactive_source_update_euler(mom_em, gas_gamma, specific_heat_capacity, energy_of_formation, ignition_temperature, reaction_rate,
      t_curr + (0.5 * dt), dt, f_stage2, f_new);
    for (int j = 0; j < 6; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
    }
  }
}

void
explicit_medium_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, const double kappa, double t_curr,
  const double dt, double* fluid_old, double* fluid_new)
{
  double exp_2a = fluid_old[0];

  double a_dt = fluid_old[1], a_dx = fluid_old[2];
  double b_dt = fluid_old[3], b_dx = fluid_old[4];
  double c_dt = fluid_old[5], c_dx = fluid_old[6];
  
  double a_dt_dx = fluid_old[7], a_dx_dx = fluid_old[8];
  double b_dt_dx = fluid_old[9], b_dx_dx = fluid_old[10];
  double c_dt_dx = fluid_old[11], c_dx_dx = fluid_old[12];

  double Etot = fluid_old[13];
  double mom = fluid_old[14];

  double rho = (1.0 / (gas_gamma - 1.0)) * ((-0.5 * (2.0 - gas_gamma) * Etot) + sqrt((0.25 * (2.0 - gas_gamma) * (2.0 - gas_gamma) * Etot * Etot) +
    ((gas_gamma - 1.0) * ((Etot * Etot) - (mom * mom)))));

  double vel = 0.0;
  if (fabs(mom) > pow(10.0, -8.0)) {
    vel = ((gas_gamma * rho) / (2.0 * mom)) * (sqrt(1.0 + ((4 * mom * mom) / ((gas_gamma * gas_gamma) * (rho * rho)))) - 1.0);
  }

  double p = (gas_gamma - 1.0) * rho;

  double W = 1.0 / sqrt(1.0 - (vel * vel));
  if (vel * vel > 1.0 - pow(1.0, -8.0)) {
    W = 1.0 / sqrt(pow(1.0, -8.0));
  }

  for (int i = 0; i < 15; i++) {
    fluid_new[i] = fluid_old[i];
  }

  fluid_new[0] += dt * (2.0 * a_dt * exp_2a);

  fluid_new[1] += dt * (a_dx_dx + (b_dt * b_dt) - (b_dx * b_dx) - (c_dt * c_dt) + (c_dx * c_dx) - (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p))));
  fluid_new[2] += dt * (a_dt_dx);
  fluid_new[3] += dt * (b_dx_dx - (2.0 * (b_dt * b_dt)) + (2.0 * (b_dx * b_dx)) + (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p))));
  fluid_new[4] += dt * (b_dt_dx);
  fluid_new[5] += dt * (c_dx_dx - (2.0 * ((b_dt * c_dt) - (b_dx * c_dx))));
  fluid_new[6] += dt * (c_dt_dx);

  fluid_new[7] += dt * ((2.0 * (b_dt * b_dt_dx)) - (2.0 * (b_dx * b_dx_dx)) - (2.0 * (c_dt * c_dt_dx)) + (2.0 * (c_dx * c_dx_dx)));
  fluid_new[8] += 0.0;
  fluid_new[9] += dt * (-(4.0 * (b_dt * b_dt_dx)) + (4.0 * (b_dx * b_dx_dx)));
  fluid_new[10] += 0.0;
  fluid_new[11] += dt * (-2.0 * ((b_dt * c_dt_dx) - (b_dx * c_dx_dx) + (b_dt_dx * c_dt) - (b_dx_dx * c_dx)));
  fluid_new[12] += 0.0;

  fluid_new[13] += dt * ((-Etot * (a_dt + (2.0 * b_dt))) - (2.0 * mom * (a_dx + b_dx)) - (((mom * vel) + p) * a_dt) - (2.0 * p * b_dt));
  fluid_new[14] += dt * ((-Etot * a_dx) - (2.0 * mom * (a_dt + b_dt)) - (((mom * vel) + p) * (a_dx + (2.0 * b_dx))) + (2.0 * p * b_dx));
}

void
explicit_medium_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  double gas_gamma = mom_em->medium_gas_gamma;
  double kappa = mom_em->medium_kappa;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double f_new[15], f_stage1[15], f_stage2[15], f_old[15];

    for (int j = 0; j < 15; j++) {
      f_old[j] = f[j];
    }

    explicit_medium_source_update_euler(mom_em, gas_gamma, kappa, t_curr, dt, f_old, f_new);
    for (int j = 0; j < 15; j++) {
      f_stage1[j] = f_new[j];
    }

    explicit_medium_source_update_euler(mom_em, gas_gamma, kappa, t_curr + dt, dt, f_stage1, f_new);
    for (int j = 0; j < 15; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
    }

    explicit_medium_source_update_euler(mom_em, gas_gamma, kappa, t_curr + (0.5 * dt), dt, f_stage2, f_new);
    for (int j = 0; j < 15; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
    }
  }
}

void
explicit_gr_ultra_rel_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, double t_curr, const double dt,
  double* fluid_old, double* fluid_new)
{
  double lapse = fluid_old[4];
  double shift_x = fluid_old[5];
  double shift_y = fluid_old[6];
  double shift_z = fluid_old[7];

  double spatial_metric[3][3];
  spatial_metric[0][0] = fluid_old[8]; spatial_metric[0][1] = fluid_old[9]; spatial_metric[0][2] = fluid_old[10];
  spatial_metric[1][0] = fluid_old[11]; spatial_metric[1][1] = fluid_old[12]; spatial_metric[1][2] = fluid_old[13];
  spatial_metric[2][0] = fluid_old[14]; spatial_metric[2][1] = fluid_old[15]; spatial_metric[2][2] = fluid_old[16];

  double inv_spatial_metric[3][3];
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
      inv_spatial_metric[i][j] = (1.0 / spatial_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);
    }
  }

  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = fluid_old[17]; extrinsic_curvature[0][1] = fluid_old[18]; extrinsic_curvature[0][2] = fluid_old[19];
  extrinsic_curvature[1][0] = fluid_old[20]; extrinsic_curvature[1][1] = fluid_old[21]; extrinsic_curvature[1][2] = fluid_old[22];
  extrinsic_curvature[2][0] = fluid_old[23]; extrinsic_curvature[2][1] = fluid_old[24]; extrinsic_curvature[2][2] = fluid_old[25];

  bool in_excision_region = false;
  if (fluid_old[26] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

    double Etot = fluid_old[0] / sqrt(spatial_det);
    double momx = fluid_old[1] / sqrt(spatial_det);
    double momy = fluid_old[2] / sqrt(spatial_det);
    double momz = fluid_old[3] / sqrt(spatial_det);

    double mom[3];
    double mom_sq = 0.0;
    mom[0] = momx; mom[1] = momy; mom[2] = momz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        mom_sq += inv_spatial_metric[i][j] * mom[i] * mom[j];
      }
    }

    double beta = 0.25 * (2.0 - gas_gamma);
    double p = -(2.0 * beta * Etot) + sqrt((4.0 * (beta * beta) * (Etot * Etot)) + ((gas_gamma - 1.0) * ((Etot * Etot) - mom_sq)));
    if (p < pow(10.0, -8.0)) {
      p = pow(10.0, -8.0);
    }
    double rho = p / (gas_gamma - 1.0);

    double cov_vx = momx / (Etot + p);
    double cov_vy = momy / (Etot + p);
    double cov_vz = momz / (Etot + p);

    if (Etot + p < pow(10.0, -8.0)) {
      cov_vx = momx / pow(10.0, -8.0);
      cov_vy = momy / pow(10.0, -8.0);
      cov_vz = momz / pow(10.0, -8.0);
    }

    if (cov_vx > 1.0 - pow(10.0, -8.0)) {
      cov_vx = 1.0 - pow(10.0, -8.0);
    }
    if (cov_vy > 1.0 - pow(10.0, -8.0)) {
      cov_vy = 1.0 - pow(10.0, -8.0);
    }
    if (cov_vz > 1.0 - pow(10.0, -8.0)) {
      cov_vz = 1.0 - pow(10.0, -8.0);
    }

    double cov_vel[3];
    cov_vel[0] = cov_vx; cov_vel[1] = cov_vy; cov_vel[2] = cov_vz;

    double vel[3];
    for (int i = 0; i < 3; i++) {
      vel[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        vel[i] += inv_spatial_metric[i][j] * cov_vel[j];
      }
    }

    double vx = vel[0];
    double vy = vel[1];
    double vz = vel[2];

    double v_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W = 1.0 / sqrt(1.0 - v_sq);
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

    double spacetime_vel[4];
    spacetime_vel[0] = W / lapse;
    spacetime_vel[1] = (W * vx) - (shift_x * (W / lapse));
    spacetime_vel[2] = (W * vy) - (shift_y * (W / lapse));
    spacetime_vel[3] = (W * vz) - (shift_z * (W / lapse));

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

    double stress_energy[4][4];
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        stress_energy[i][j] = (rho * h * spacetime_vel[i] * spacetime_vel[j]) + (p * inv_spacetime_metric[i][j]);
      }
    }

    double lapse_der[3];
    lapse_der[0] = fluid_old[27];
    lapse_der[1] = fluid_old[28];
    lapse_der[2] = fluid_old[29];

    double shift_der[3][3];
    shift_der[0][0] = fluid_old[30]; shift_der[0][1] = fluid_old[31]; shift_der[0][2] = fluid_old[32];
    shift_der[1][0] = fluid_old[33]; shift_der[1][1] = fluid_old[34]; shift_der[1][2] = fluid_old[35];
    shift_der[2][0] = fluid_old[36]; shift_der[2][1] = fluid_old[37]; shift_der[2][2] = fluid_old[38];

    double spatial_metric_der[3][3][3];
    spatial_metric_der[0][0][0] = fluid_old[39]; spatial_metric_der[0][0][1] = fluid_old[40]; spatial_metric_der[0][0][2] = fluid_old[41];
    spatial_metric_der[0][1][0] = fluid_old[42]; spatial_metric_der[0][1][1] = fluid_old[43]; spatial_metric_der[0][1][2] = fluid_old[44];
    spatial_metric_der[0][2][0] = fluid_old[45]; spatial_metric_der[0][2][1] = fluid_old[46]; spatial_metric_der[0][2][2] = fluid_old[47];

    spatial_metric_der[1][0][0] = fluid_old[48]; spatial_metric_der[1][0][1] = fluid_old[49]; spatial_metric_der[1][0][2] = fluid_old[50];
    spatial_metric_der[1][1][0] = fluid_old[51]; spatial_metric_der[1][1][1] = fluid_old[52]; spatial_metric_der[1][1][2] = fluid_old[53];
    spatial_metric_der[1][2][0] = fluid_old[54]; spatial_metric_der[1][2][1] = fluid_old[55]; spatial_metric_der[1][2][2] = fluid_old[56];

    spatial_metric_der[0][0][0] = fluid_old[57]; spatial_metric_der[0][0][1] = fluid_old[58]; spatial_metric_der[0][0][2] = fluid_old[59];
    spatial_metric_der[0][1][0] = fluid_old[60]; spatial_metric_der[0][1][1] = fluid_old[61]; spatial_metric_der[0][1][2] = fluid_old[62];
    spatial_metric_der[0][2][0] = fluid_old[63]; spatial_metric_der[0][2][1] = fluid_old[64]; spatial_metric_der[0][2][2] = fluid_old[65];

    mom[0] = (rho + p) * (W * W) * vx;
    mom[1] = (rho + p) * (W * W) * vy;
    mom[2] = (rho + p) * (W * W) * vz;

    for (int i = 0; i < 70; i++) {
      fluid_new[i] = fluid_old[i];
    }

    // Energy density source.
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        fluid_new[0] += dt * (stress_energy[0][0] * shift[i] * shift[j] * extrinsic_curvature[i][j]);
        fluid_new[0] += dt * (2.0 * stress_energy[0][i + 1] * shift[j] * extrinsic_curvature[i][j]);
        fluid_new[0] += dt * (stress_energy[i + 1][j + 1] * extrinsic_curvature[i][j]);
      }

      fluid_new[0] -= dt * (stress_energy[0][0] * shift[i] * lapse_der[i]);
      fluid_new[0] -= dt * (stress_energy[0][i + 1] * lapse_der[i]);
    }

    // Momentum density sources.
    for (int j = 0; j < 3; j++) {
      fluid_new[1 + j] -= dt * (stress_energy[0][0] * lapse * lapse_der[j]);

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          fluid_new[1 + j] += dt * (0.5 * stress_energy[0][0] * shift[k] * shift[l] * spatial_metric_der[j][k][l]);
          fluid_new[1 + j] += dt * (0.5 * stress_energy[k + 1][l + 1] * spatial_metric_der[j][k][l]);
        }

        fluid_new[1 + j] += dt * ((mom[k] / lapse) * shift_der[j][k]);

        for (int i = 0; i < 3; i++) {
          fluid_new[1 + j] += dt * (stress_energy[0][i + 1] * shift[k] * spatial_metric_der[j][i][k]);
        }
      }
    }
  }
  else {
    for (int i = 0; i < 70; i++) {
      fluid_new[i] = fluid_old[i];
    }
  }
}

void
explicit_gr_ultra_rel_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  double gas_gamma = mom_em->gr_ultra_rel_gas_gamma;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double f_new[70], f_stage1[70], f_stage2[70], f_old[70];

    for (int j = 0; j < 70; j++) {
      f_old[j] = f[j];
    }

    explicit_gr_ultra_rel_source_update_euler(mom_em, gas_gamma, t_curr, dt, f_old, f_new);
    for (int j = 0; j < 70; j++) {
      f_stage1[j] = f_new[j];
    }

    explicit_gr_ultra_rel_source_update_euler(mom_em, gas_gamma, t_curr + dt, dt, f_stage1, f_new);
    for (int j = 0; j < 70; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
    }

    explicit_gr_ultra_rel_source_update_euler(mom_em, gas_gamma, t_curr + (0.5 * dt), dt, f_stage2, f_new);
    for (int j = 0; j < 70; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
    }
  }
}

void
explicit_e_field_source_update_euler(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double e_field_old[3], double* e_field_new,
  double* fluid_s[GKYL_MAX_SPECIES], const double* app_current)
{
  int nfluids = mom_em->nfluids;

  double epsilon0 = mom_em->epsilon0;
  double mu0 = mom_em->mu0;
  double c = 1.0 / sqrt(mu0 * epsilon0);

  e_field_new[0] = e_field_old[0] + (dt * (-(1.0 / epsilon0) * app_current[0]));
  e_field_new[1] = e_field_old[1] + (dt * (-(1.0 / epsilon0) * app_current[1]));
  e_field_new[2] = e_field_old[2] + (dt * (-(1.0 / epsilon0) * app_current[2]));

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];
    const double q = mom_em->param[i].charge;

    double rho = f[0];

    if (rho > 0.0) {
      double ux = f[1] / rho;
      double uy = f[2] / rho;
      double uz = f[3] / rho;

      double gamma = sqrt(1.0 + (((ux * ux) + (uy * uy) + (uz * uz)) / (c * c)));

      double vx = ux / gamma;
      double vy = uy / gamma;
      double vz = uz / gamma;

      e_field_new[0] += dt * (-(1.0 / epsilon0) * q * rho * vx);
      e_field_new[1] += dt * (-(1.0 / epsilon0) * q * rho * vy);
      e_field_new[2] += dt * (-(1.0 / epsilon0) * q * rho * vz);
    }
  }
}

void
explicit_e_field_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES], double* em,
  const double* app_current, const double* app_current1, const double* app_current2, const double* ext_em)
{
  double e_field_new[3], e_field_stage1[3], e_field_stage2[3];
  double e_field_old[3];
  e_field_old[0] = em[0]; e_field_old[1] = em[1]; e_field_old[2] = em[2];

  explicit_e_field_source_update_euler(mom_em, t_curr, dt, e_field_old, e_field_new, fluid_s, app_current);
  e_field_stage1[0] = e_field_new[0];
  e_field_stage1[1] = e_field_new[1];
  e_field_stage1[2] = e_field_new[2];

  explicit_e_field_source_update_euler(mom_em, t_curr + dt, dt, e_field_stage1, e_field_new, fluid_s, app_current1);
  e_field_stage2[0] = (0.75 * e_field_old[0]) + (0.25 * e_field_new[0]);
  e_field_stage2[1] = (0.75 * e_field_old[1]) + (0.25 * e_field_new[1]);
  e_field_stage2[2] = (0.75 * e_field_old[2]) + (0.25 * e_field_new[2]);

  explicit_e_field_source_update_euler(mom_em, t_curr + (0.5 * dt), dt, e_field_stage2, e_field_new, fluid_s, app_current2);
  em[0] = ((1.0 / 3.0) * e_field_old[0]) + ((2.0 / 3.0) * e_field_new[0]);
  em[1] = ((1.0 / 3.0) * e_field_old[1]) + ((2.0 / 3.0) * e_field_new[1]);
  em[2] = ((1.0 / 3.0) * e_field_old[2]) + ((2.0 / 3.0) * e_field_new[2]);
}

void
explicit_higuera_cary_push(double* vel, const double q, const double m, const double dt, const double c, const double e_field[3],
  const double b_field[3])
{
  const double q_over_m = (0.5 * dt) * (q / m);
  const double Ex = q_over_m * e_field[0];
  const double Ey = q_over_m * e_field[1];
  const double Ez = q_over_m * e_field[2];
  const double Bx = q_over_m * b_field[0];
  const double By = q_over_m * b_field[1];
  const double Bz = q_over_m * b_field[2];

  const double vel_x_minus = vel[0] + Ex;
  const double vel_y_minus = vel[1] + Ey;
  const double vel_z_minus = vel[2] + Ez;

  const double vel_star = (vel_x_minus * (Bx / c)) + (vel_y_minus * (By / c)) + (vel_z_minus * (Bz / c));
  const double gamma_minus = sqrt(1.0 + (((vel_x_minus * vel_x_minus) + (vel_y_minus * vel_y_minus) + (vel_z_minus * vel_z_minus)) / (c * c)));
  const double dot_tau_tau = (Bx * Bx) + (By * By) + (Bz * Bz);
  const double sigma = (gamma_minus * gamma_minus) - dot_tau_tau;
  const double gamma_new = sqrt(0.5 * (sigma + sqrt((sigma * sigma) + (4.0 * (dot_tau_tau + (vel_star * vel_star))))));

  const double tx = Bx / gamma_new;
  const double ty = By / gamma_new;
  const double tz = Bz / gamma_new;
  const double s =  1.0 / (1.0 + ((tx * tx) + (ty * ty) + (tz * tz)));

  const double t_vel_minus = (vel_x_minus * tx) + (vel_y_minus * ty) + (vel_z_minus * tz);
  const double vel_x_plus = s * (vel_x_minus + (t_vel_minus * tx) + ((vel_y_minus * tz) - (vel_z_minus * ty)));
  const double vel_y_plus = s * (vel_y_minus + (t_vel_minus * ty) + ((vel_z_minus * tx) - (vel_x_minus * tz)));
  const double vel_z_plus = s * (vel_z_minus + (t_vel_minus * tz) + ((vel_x_minus * ty) - (vel_y_minus * tx)));

  vel[0] = vel_x_plus + Ex + ((vel_y_plus * tz) - (vel_z_plus * ty));
  vel[1] = vel_y_plus + Ey + ((vel_z_plus * tx) - (vel_x_plus * tz));
  vel[2] = vel_z_plus + Ez + ((vel_x_plus * ty) - (vel_y_plus * tx));
}

void
explicit_higuera_cary_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], double* em, const double* ext_em)
{
  int nfluids = mom_em->nfluids;

  double epsilon0 = mom_em->epsilon0;
  double mu0 = mom_em->mu0;
  double c = 1.0 / sqrt(mu0 * epsilon0);

  double Ex = em[0] + ext_em[0];
  double Ey = em[1] + ext_em[1];
  double Ez = em[2] + ext_em[2];
  double Bx = em[3] + ext_em[3];
  double By = em[4] + ext_em[4];
  double Bz = em[5] + ext_em[5];

  double e_field[3], b_field[3];
  e_field[0] = Ex; e_field[1] = Ey; e_field[2] = Ez;
  b_field[0] = Bx; b_field[1] = By; b_field[2] = Bz;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];
    const double q = mom_em->param[i].charge;
    const double m = mom_em->param[i].mass;

    double rho = f[0];

    if (rho > 0.0) {
      double vel[3];
      vel[0] = f[1] / rho;
      vel[1] = f[2] / rho;
      vel[2] = f[3] / rho;

      explicit_higuera_cary_push(vel, q, m, dt, c, e_field, b_field);

      f[1] = rho * vel[0];
      f[2] = rho * vel[1];
      f[3] = rho * vel[2];
    }
  }
}

void
explicit_source_coupling_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], double* em, const double* app_current, const double* app_current1, const double* app_current2,
  const double* ext_em, int nstrang)
{
  if (nstrang == 0) {
    explicit_e_field_source_update(mom_em, t_curr, dt, fluid_s, em, app_current, app_current1, app_current2, ext_em);
  }
  else if (nstrang == 1) {
    explicit_higuera_cary_update(mom_em, t_curr, dt, fluid_s, app_accel_s, em, ext_em);
  }
}