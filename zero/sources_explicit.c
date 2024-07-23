#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_fv_proj.h>
#include <gkyl_sources_explicit.h>
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

    double E_elc = f_elc_old[4];
    double E_ion = f_ion_old[4];

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

    f_elc_new[4] = f_elc_old[4] + (dt * E_src);
    f_ion_new[4] = f_ion_old[4] - (dt * E_src);

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

    double f_elc_new[5], f_elc_stage1[5], f_elc_stage2[5];
    double f_ion_new[5], f_ion_stage1[5], f_ion_stage2[5];
    double f_elc_old[5], f_ion_old[5];

    for (int i = 0; i < 5; i++) {
      f_elc_old[i] = f_elc[i];
      f_ion_old[i] = f_ion[i];
    }

    explicit_frictional_source_update_euler(mom_em, Z, T_elc, Lambda_ee, t_curr, dt, f_elc_old, f_ion_old, f_elc_new, f_ion_new);
    for (int i = 0; i < 5; i++) {
      f_elc_stage1[i] = f_elc_new[i];
      f_ion_stage1[i] = f_ion_new[i];
    }

    explicit_frictional_source_update_euler(mom_em, Z, T_elc, Lambda_ee, t_curr + dt, dt, f_elc_stage1, f_ion_stage1, f_elc_new, f_ion_new);
    for (int i = 0; i < 5; i++) {
      f_elc_stage2[i] = (0.75 * f_elc_old[i]) + (0.25 * f_elc_new[i]);
      f_ion_stage2[i] = (0.75 * f_ion_old[i]) + (0.25 * f_ion_new[i]);
    }

    explicit_frictional_source_update_euler(mom_em, Z, T_elc, Lambda_ee, t_curr + (0.5 * dt), dt, f_elc_stage2, f_ion_stage2, f_elc_new, f_ion_new);
    for (int i = 0; i < 5; i++) {
      f_elc[i] = ((1.0 / 3.0) * f_elc_old[i]) + ((2.0 / 3.0) * f_elc_new[i]);
      f_ion[i] = ((1.0 / 3.0) * f_ion_old[i]) + ((2.0 / 3.0) * f_ion_new[i]);
    }
  }
}

void
explicit_volume_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, const double U0, const double R0,
  double t_curr, const double dt, double* fluid_old, double* fluid_new, double* em, const double* ext_em)
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
explicit_volume_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES],
  double* em, const double* ext_em)
{
  int nfluids = mom_em->nfluids;

  double gas_gamma = mom_em->volume_gas_gamma;
  double U0 = mom_em->volume_U0;
  double R0 = mom_em->volume_R0;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double f_new[5], f_stage1[5], f_stage2[5], f_old[5];

    for (int j = 0; j < 5; j++) {
      f_old[j] = f[j];
    }

    explicit_volume_source_update_euler(mom_em, gas_gamma, U0, R0, t_curr, dt, f_old, f_new, em, ext_em);
    for (int j = 0; j < 5; j++) {
      f_stage1[j] = f_new[j];
    }

    explicit_volume_source_update_euler(mom_em, gas_gamma, U0, R0, t_curr + dt, dt, f_stage1, f_new, em, ext_em);
    for (int j = 0; j < 5; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
    }

    explicit_volume_source_update_euler(mom_em, gas_gamma, U0, R0, t_curr + (0.5 * dt), dt, f_stage2, f_new, em, ext_em);
    for (int j = 0; j < 5; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
    }
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