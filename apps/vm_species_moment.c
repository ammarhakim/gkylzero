#include <assert.h>
#include <gkyl_vlasov_priv.h>

// initialize species moment object
void
vm_species_moment_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_species_moment *sm, const char *nm)
{
  assert(is_moment_name_valid(nm));

  bool is_integrated = strcmp(nm, "Integrated") == 0;
  sm->nm = nm;
  sm->vdim = app->vdim;
  int num_mom;
  sm->is_sr_five_moments = strcmp("SRFiveMoments", nm) == 0;

  if (s->model_id == GKYL_MODEL_SR) {
    if (sm->is_sr_five_moments){
      sm->n = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
      sm->vbi = mkarr(false, app->vdim*app->confBasis.num_basis, app->local_ext.volume);
      sm->T = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
      sm->mj_moms = gkyl_mj_moments_new(&s->grid, &app->confBasis,
        &app->basis, &app->local, &s->local_vel, app->local.volume, app->local_ext.volume, 
        s->p_over_gamma, s->gamma, s->gamma_inv, false);
      sm->num_basis = app->confBasis.num_basis;
      num_mom = 2 + sm->vdim;
    }
    else {
      struct gkyl_mom_vlasov_sr_auxfields sr_inp = {.p_over_gamma = s->p_over_gamma, 
        .gamma = s->gamma, .gamma_inv = s->gamma_inv, .V_drift = s->V_drift, 
        .GammaV2 = s->GammaV2, .GammaV_inv = s->GammaV_inv};
      sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &app->confBasis, 
        &app->basis, &app->local, &s->local_vel, s->model_id, &sr_inp, 
        nm, is_integrated, s->info.mass, app->use_gpu);
      num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc);
    }
  }
  else {
    // No auxiliary fields for moments if not SR 
    sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &app->confBasis, 
      &app->basis, &app->local, &s->local_vel, s->model_id, 0, 
      nm, is_integrated, s->info.mass, app->use_gpu);   
    num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc); 
  }

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
  if (sm->is_sr_five_moments) {
    //gkyl_mj_moments_advance(sm->mj_moms, &phase_rng, &conf_rng, fin, sm->marr);
    gkyl_mj_moments_advance(sm->mj_moms, fin, sm->n, sm->vbi, sm->T, &phase_rng, &conf_rng);

    // Save the outputs to n vbi T:
    gkyl_array_set_offset(sm->marr, 1.0, sm->n, 0*sm->num_basis);
    gkyl_array_set_offset(sm->marr, 1.0, sm->vbi, 1*sm->num_basis);
    gkyl_array_set_offset(sm->marr, 1.0, sm->T, (1 + sm->vdim)*sm->num_basis);
  }
  else {
    gkyl_dg_updater_moment_advance(sm->mcalc, &phase_rng, &conf_rng, fin, sm->marr);
  }
}

// release memory for moment data object
void
vm_species_moment_release(const struct gkyl_vlasov_app *app, const struct vm_species_moment *sm)
{
  if (app->use_gpu) {
    gkyl_array_release(sm->marr_host);
  }

  gkyl_array_release(sm->marr);

  if(sm->is_sr_five_moments) {
    gkyl_mj_moments_release(sm->mj_moms);
    gkyl_array_release(sm->n);
    gkyl_array_release(sm->vbi);
    gkyl_array_release(sm->T);
  }
  else {
    gkyl_dg_updater_moment_release(sm->mcalc);
  }
}
