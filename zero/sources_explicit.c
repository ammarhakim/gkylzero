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
explicit_source_coupling_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], double* em, const double* app_current, const double* app_current1, const double* app_current2,
  const double* ext_em, int nstrang)
{
  if (nstrang == 0) {
    explicit_e_field_source_update(mom_em, t_curr, dt, fluid_s, em, app_current, app_current1, app_current2, ext_em);
  }
  else if (nstrang == 1) {
    // TODO: Add explicit Higuera-Cary update step.
  }
}