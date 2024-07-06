#include <assert.h>
#include <gkyl_vlasov_priv.h>

// initialize species moment object
void
vm_species_moment_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_species_moment *sm, const char *nm)
{
  assert(is_moment_name_valid(nm));

  bool is_integrated = strcmp(nm, "Integrated") == 0;
  int num_mom;

  sm->is_vlasov_lte_moms = strcmp("LTEMoments", nm) == 0;
  if (sm->is_vlasov_lte_moms) {
    struct gkyl_vlasov_lte_moments_inp inp_mom = {
      .phase_grid = &s->grid,
      .vel_grid = &s->grid_vel, 
      .conf_basis = &app->confBasis,
      .vel_basis = &app->velBasis, 
      .phase_basis = &app->basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .vel_range = &s->local_vel,
      .gamma = s->gamma,
      .gamma_inv = s->gamma_inv,
      .h_ij_inv = s->h_ij_inv,
      .det_h = s->det_h,
      .model_id = s->model_id,
      .use_gpu = app->use_gpu,
    };
    // Compute (n, V_drift, T/m)
    sm->vlasov_lte_moms = gkyl_vlasov_lte_moments_inew(&inp_mom);
    // If relativistic, V_drift is the bulk four-velocity (GammaV, GammaV*V_drift)
    if (s->model_id == GKYL_MODEL_SR) {
      num_mom = app->vdim + 3; 
    }
    else {
      num_mom = app->vdim + 2;
    }
  }
  else {
    if (s->model_id == GKYL_MODEL_SR) {
      struct gkyl_mom_vlasov_sr_auxfields sr_inp = {.gamma = s->gamma};
      sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &app->confBasis, 
        &app->basis, &app->local, &s->local_vel, s->model_id, &sr_inp, 
        nm, is_integrated, app->use_gpu);
      num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc);
    }
    else {
      // No auxiliary fields for moments if not SR 
      sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &app->confBasis, 
        &app->basis, &app->local, &s->local_vel, s->model_id, 0, 
        nm, is_integrated, app->use_gpu);   
      num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc); 
    }
  }

  if (is_integrated) {
    sm->marr = mkarr(app->use_gpu, num_mom, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu) {
      sm->marr_host = mkarr(false, num_mom, app->local_ext.volume);      
    }
  }
  else {
    sm->marr = mkarr(app->use_gpu, num_mom*app->confBasis.num_basis,
      app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu) {
      sm->marr_host = mkarr(false, num_mom*app->confBasis.num_basis,
        app->local_ext.volume);
    }
  }
}

void
vm_species_moment_calc(const struct vm_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  if (sm->is_vlasov_lte_moms) {
    gkyl_vlasov_lte_moments_advance(sm->vlasov_lte_moms, 
      &phase_rng, &conf_rng, fin, sm->marr);
  }
  else {
    gkyl_dg_updater_moment_advance(sm->mcalc, 
      &phase_rng, &conf_rng, fin, sm->marr);
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

  if(sm->is_vlasov_lte_moms) {
    gkyl_vlasov_lte_moments_release(sm->vlasov_lte_moms);
  }
  else {
    gkyl_dg_updater_moment_release(sm->mcalc);
  }
}
