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
explicit_source_coupling_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], double* em, const double* app_current, const double* app_current1, const double* app_current2,
  const double* ext_em, int nstrang)
{
  if (nstrang == 0) {
    // TODO: Add explicit source updater for electric fields.
  }
  else if (nstrang == 1) {
    // TODO: Add explicit Higuera-Cary update step.
  }
}