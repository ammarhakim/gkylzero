#include <assert.h>
#include <gkyl_vlasov_poisson_priv.h>

void
vp_species_moment_init(struct gkyl_vlasov_poisson_app *app, struct vp_species *vps,
  struct vp_species_moment *sm, const char *nm)
{
  assert(is_moment_name_valid(nm));

  sm->is_integrated = strcmp(nm, "Integrated") == 0;

  sm->is_vlasov_lte_moms = strcmp("LTEMoments", nm) == 0;
  if (sm->is_vlasov_lte_moms) {
    struct gkyl_vlasov_lte_moments_inp inp_mom = {
      .phase_grid = &vps->grid,
      .conf_basis = &app->confBasis,
      .phase_basis = &app->basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .vel_range = &vps->local_vel,
      .model_id = vps->model_id,
      .mass = vps->info.mass,
      .use_gpu = app->use_gpu,
    };
    sm->vlasov_lte_moms = gkyl_vlasov_lte_moments_inew(  &inp_mom  );
    sm->num_mom = app->vdim + 2;
  }
  else {
    sm->mcalc = gkyl_dg_updater_moment_new(&vps->grid, &app->confBasis, 
      &app->basis, &app->local, &vps->local_vel, vps->model_id, 0, 
      nm, sm->is_integrated, app->use_gpu);   
    sm->num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc); 
  }

  if (sm->is_integrated) {
    sm->marr = mkarr(app->use_gpu, sm->num_mom, app->local_ext.volume);
    sm->marr_host = app->use_gpu? mkarr(false, sm->num_mom, app->local_ext.volume)
	                        : gkyl_array_acquire(sm->marr);
  }
  else {
    sm->marr = mkarr(app->use_gpu, sm->num_mom*app->confBasis.num_basis, app->local_ext.volume);
    sm->marr_host = app->use_gpu? mkarr(false, sm->num_mom*app->confBasis.num_basis, app->local_ext.volume)
	                        : gkyl_array_acquire(sm->marr);
  }
}

void
vp_species_moment_calc(const struct vp_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  if (sm->is_vlasov_lte_moms) {
    gkyl_vlasov_lte_moments_advance(sm->vlasov_lte_moms, 
      &phase_rng, &conf_rng, fin, sm->marr);
  }
  else {
    gkyl_dg_updater_moment_advance(sm->mcalc, &phase_rng, &conf_rng, fin, sm->marr);
  }
}

void
vp_species_moment_release(const struct gkyl_vlasov_poisson_app *app, const struct vp_species_moment *sm)
{
  // Release memory for moment data object.
  gkyl_array_release(sm->marr_host);
  gkyl_array_release(sm->marr);

  if(sm->is_vlasov_lte_moms) {
    gkyl_vlasov_lte_moments_release(sm->vlasov_lte_moms);
  }
  else {
    gkyl_dg_updater_moment_release(sm->mcalc);
  }
}
