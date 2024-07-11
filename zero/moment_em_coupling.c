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
gkyl_moment_em_coupling_release(gkyl_moment_em_coupling* mom_em)
{
  gkyl_free(mom_em);
}