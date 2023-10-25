#include <assert.h>
#include <gkyl_vlasov_priv.h>

// initialize species moment object
void
vm_species_moment_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_species_moment *sm, const char *nm)
{
  assert(is_moment_name_valid(nm));

  bool is_integrated = strcmp(nm, "Integrated") == 0;

  if (s->model_id == GKYL_MODEL_SR) {
    struct gkyl_mom_vlasov_sr_auxfields sr_inp = {.p_over_gamma = s->p_over_gamma, 
      .gamma = s->gamma, .gamma_inv = s->gamma_inv, .V_drift = s->V_drift, 
      .GammaV2 = s->GammaV2, .GammaV_inv = s->GammaV_inv};
    sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &app->confBasis, 
      &app->basis, &app->local, &s->local_vel, s->model_id, &sr_inp, 
      nm, is_integrated, s->info.mass, app->use_gpu);
  }
  else {
    // No auxiliary fields for moments if not SR 
    sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &app->confBasis, 
      &app->basis, &app->local, &s->local_vel, s->model_id, 0, 
      nm, is_integrated, s->info.mass, app->use_gpu);    
  }

  int num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc);

  if (is_integrated) {
    sm->marr = mkarr(app->use_gpu, num_mom, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu)
      sm->marr_host = mkarr(false, num_mom, app->local_ext.volume);      
  }
  else {
    sm->marr = mkarr(app->use_gpu, num_mom*app->confBasis.num_basis,
      app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu)
      sm->marr_host = mkarr(false, num_mom*app->confBasis.num_basis,
        app->local_ext.volume);
  }
}

void
vm_species_moment_calc(const struct vm_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  gkyl_dg_updater_moment_advance(sm->mcalc, &phase_rng, &conf_rng, fin, sm->marr);
}

// release memory for moment data object
void
vm_species_moment_release(const struct gkyl_vlasov_app *app, const struct vm_species_moment *sm)
{
  if (app->use_gpu)
    gkyl_array_release(sm->marr_host);

  gkyl_dg_updater_moment_release(sm->mcalc);
  gkyl_array_release(sm->marr);
}
