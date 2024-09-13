#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

// initialize neutral species moment object
void
gk_neut_species_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_species_moment *sm, const char *nm)
{
  assert(is_neut_moment_name_valid(nm));

  sm->is_integrated = strcmp(nm, "Integrated") == 0;

  sm->is_maxwellian_moms = strcmp("LTEMoments", nm) == 0;
  if (sm->is_maxwellian_moms) {
    struct gkyl_vlasov_lte_moments_inp inp_mom = {
      .phase_grid = &s->grid,
      .vel_grid = &s->grid_vel, 
      .conf_basis = &app->confBasis,
      .phase_basis = &app->neut_basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .vel_range = &s->local_vel,
      .use_gpu = app->use_gpu,
    };
    // Compute (n, ux, uy, uz, T/m)
    sm->vlasov_lte_moms = gkyl_vlasov_lte_moments_inew(&inp_mom);
    sm->num_mom = 5; 
  }
  else {
    sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &app->confBasis, 
      &app->neut_basis, &app->local, &s->local_vel, s->model_id, 0,
      nm, sm->is_integrated, app->use_gpu);    

    sm->num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc);
  }

  if (sm->is_integrated) {
    sm->marr = mkarr(app->use_gpu, sm->num_mom, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu) {
      sm->marr_host = mkarr(false, sm->num_mom, app->local_ext.volume); 
    }
  }
  else {
    sm->marr = mkarr(app->use_gpu, sm->num_mom*app->confBasis.num_basis, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu)
      sm->marr_host = mkarr(false, sm->num_mom*app->confBasis.num_basis, app->local_ext.volume);
    // Bin Op memory for rescaling moment by inverse of Jacobian
    if (app->use_gpu) {
      sm->mem_geo = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
    }
    else {
      sm->mem_geo = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);
    }
  }
}

void
gk_neut_species_moment_calc(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  if (sm->is_maxwellian_moms) {
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
gk_neut_species_moment_release(const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm)
{
  if (app->use_gpu) {
    gkyl_array_release(sm->marr_host);
  }
  gkyl_array_release(sm->marr);

  if(sm->is_maxwellian_moms) {
    gkyl_vlasov_lte_moments_release(sm->vlasov_lte_moms);
  }
  else {
    gkyl_dg_updater_moment_release(sm->mcalc);
  }

  // free the weak division memory if not computing integrated moments
  if (!sm->is_integrated) {
    gkyl_dg_bin_op_mem_release(sm->mem_geo);
  }
}
