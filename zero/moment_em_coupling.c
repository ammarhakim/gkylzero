#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_fv_proj.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_moment_em_coupling_priv.h>
#include <gkyl_mat.h>

gkyl_moment_em_coupling*
gkyl_moment_em_coupling_new(struct gkyl_moment_em_coupling_inp inp)
{
  gkyl_moment_em_coupling *mom_em = gkyl_malloc(sizeof(gkyl_moment_em_coupling));

  mom_em->grid = *(inp.grid);
  mom_em->ndim = mom_em->grid.ndim;
  mom_em->nfluids = inp.nfluids;
  for (int i = 0; i < mom_em->nfluids; i++) {
    mom_em->param[i] = inp.param[i];
  }

  mom_em->epsilon0 = inp.epsilon0;
  mom_em->mu0 = inp.mu0;
  if (mom_em->epsilon0 != 0.0) {
    mom_em->is_charged_species = true;
  }
  else {
    mom_em->is_charged_species = false;
  }

  mom_em->t_ramp_E = inp.t_ramp_E;
  if (mom_em->t_ramp_E != 0.0) {
    mom_em->ramp_app_E = true;
  }
  else {
    mom_em->ramp_app_E = false;
  }
  mom_em->t_ramp_curr = inp.t_ramp_curr;
  if (mom_em->t_ramp_curr != 0.0) {
    mom_em->ramp_app_curr = true;
  }
  else {
    mom_em->ramp_app_curr = false;
  }

  mom_em->has_collision = inp.has_collision;
  if (mom_em->has_collision) {
    for (int i = 0; i < mom_em->nfluids; i++) {
      for (int j = 0; j < mom_em->nfluids; j++) {
        mom_em->nu_base[i][j] = inp.nu_base[i][j];
      }
    }
  }
  else {
    for (int i = 0; i < mom_em->nfluids; i++) {
      for (int j = 0; j < mom_em->nfluids; j++) {
        mom_em->nu_base[i][j] = 0.0;
      }
    }
  }

  mom_em->use_explicit_em_coupling = inp.use_explicit_em_coupling;

  mom_em->has_nT_sources = inp.has_nT_sources;

  return mom_em;
}

void
gkyl_moment_em_coupling_implicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], const struct gkyl_array* p_rhs[GKYL_MAX_SPECIES],
  struct gkyl_array* em, const struct gkyl_array* app_current, const struct gkyl_array* ext_em, const struct gkyl_array* nT_sources[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;
  double *fluid_s[GKYL_MAX_SPECIES];
  const double *app_accel_s[GKYL_MAX_SPECIES];
  const double *p_rhs_s[GKYL_MAX_SPECIES];
  const double *nT_sources_s[GKYL_MAX_SPECIES];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long cell_idx = gkyl_range_idx(update_range, iter.idx);

    for (int i = 0; i < nfluids; i++) {
      fluid_s[i] = gkyl_array_fetch(fluid[i], cell_idx);
      app_accel_s[i] = gkyl_array_cfetch(app_accel[i], cell_idx);
      p_rhs_s[i] = gkyl_array_cfetch(p_rhs[i], cell_idx);
      nT_sources_s[i] = gkyl_array_cfetch(nT_sources[i], cell_idx);
    }

    double *em_arr = gkyl_array_fetch(em, cell_idx);
    const double *app_current_arr = gkyl_array_cfetch(app_current, cell_idx);
    const double *ext_em_arr = gkyl_array_cfetch(ext_em, cell_idx);

    implicit_source_coupling_update(mom_em, t_curr, dt, fluid_s, app_accel_s, p_rhs_s, em_arr, app_current_arr, ext_em_arr, nT_sources_s);
  }
}

void
gkyl_moment_em_coupling_explicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], const struct gkyl_array* p_rhs[GKYL_MAX_SPECIES],
  struct gkyl_array* em, const struct gkyl_array* app_current, const struct gkyl_array* app_current1, const struct gkyl_array* app_current2,
  const struct gkyl_array* ext_em, const struct gkyl_array* nT_sources[GKYL_MAX_SPECIES], gkyl_fv_proj *proj_app_curr, int nstrang)
{
  int nfluids = mom_em->nfluids;
  double *fluid_s[GKYL_MAX_SPECIES];
  const double *app_accel_s[GKYL_MAX_SPECIES];
  const double *p_rhs_s[GKYL_MAX_SPECIES];
  const double *nT_sources_s[GKYL_MAX_SPECIES];

  double dt_local = 2.0 * dt;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long cell_idx = gkyl_range_idx(update_range, iter.idx);

    for (int i = 0; i < nfluids; i++) {
      fluid_s[i] = gkyl_array_fetch(fluid[i], cell_idx);
      app_accel_s[i] = gkyl_array_cfetch(app_accel[i], cell_idx);
    }

    double *em_arr = gkyl_array_fetch(em, cell_idx);
    const double *app_current_arr = gkyl_array_cfetch(app_current, cell_idx);
    const double *app_current1_arr = gkyl_array_cfetch(app_current1, cell_idx);
    const double *app_current2_arr = gkyl_array_cfetch(app_current2, cell_idx);
    const double *ext_em_arr = gkyl_array_cfetch(ext_em, cell_idx);

    explicit_source_coupling_update(mom_em, t_curr, dt, fluid_s, app_accel_s, em_arr, app_current_arr, app_current1_arr, app_current2_arr,
      ext_em_arr, nstrang);
  }
}

void
gkyl_moment_em_coupling_release(gkyl_moment_em_coupling* mom_em)
{
  gkyl_free(mom_em);
}