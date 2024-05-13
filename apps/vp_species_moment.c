#include <assert.h>
#include <gkyl_vlasov_poisson_priv.h>

void
vp_species_moment_init(struct gkyl_vlasov_poisson_app *app, struct vp_species *vps,
  struct vp_species_moment *sm, const char *nm)
{
  // Initialize species moment object.

  assert(is_moment_name_valid(nm));

  sm->is_integrated = strcmp(nm, "Integrated") == 0;

  sm->mcalc = gkyl_dg_updater_moment_new(&vps->grid, &app->confBasis, 
    &app->basis, &app->local, &vps->local_vel, vps->model_id, 0, 
    nm, sm->is_integrated, vps->info.mass, app->use_gpu);   

  sm->num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc); 

  if (sm->is_integrated) {
    sm->marr = mkarr(app->use_gpu, sm->num_mom, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu)
      sm->marr_host = mkarr(false, sm->num_mom, app->local_ext.volume);      
  }
  else {
    sm->marr = mkarr(app->use_gpu, sm->num_mom*app->confBasis.num_basis, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu)
      sm->marr_host = mkarr(false, sm->num_mom*app->confBasis.num_basis, app->local_ext.volume);
  }
}

void
vp_species_moment_calc(const struct vp_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  gkyl_dg_updater_moment_advance(sm->mcalc, &phase_rng, &conf_rng, fin, sm->marr);
}

void
vp_species_moment_release(const struct gkyl_vlasov_poisson_app *app, const struct vp_species_moment *sm)
{
  // Release memory for moment data object.
  if (app->use_gpu)
    gkyl_array_release(sm->marr_host);

  gkyl_array_release(sm->marr);

  gkyl_dg_updater_moment_release(sm->mcalc);
}
