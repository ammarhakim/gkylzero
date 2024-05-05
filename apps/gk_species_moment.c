#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

// initialize species moment object
void
gk_species_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_species_moment *sm, const char *nm)
{
  assert(is_moment_name_valid(nm));

  sm->is_integrated = strcmp(nm, "Integrated") == 0;

  sm->mcalc = gkyl_dg_updater_moment_gyrokinetic_new(&s->grid, &app->confBasis, 
    &app->basis, &app->local, s->info.mass, s->vel_map, app->gk_geom, nm,
    sm->is_integrated, app->use_gpu);    

  sm->num_mom = gkyl_dg_updater_moment_gyrokinetic_num_mom(sm->mcalc);

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
    // Bin Op memory for rescaling moment by inverse of Jacobian
    if (app->use_gpu)
      sm->mem_geo = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
    else
      sm->mem_geo = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);
  }
}

void
gk_species_moment_calc(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  gkyl_dg_updater_moment_gyrokinetic_advance(sm->mcalc, &phase_rng, &conf_rng, fin, sm->marr);
}

// release memory for moment data object
void
gk_species_moment_release(const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm)
{
  if (app->use_gpu)
    gkyl_array_release(sm->marr_host);

  gkyl_dg_updater_moment_gyrokinetic_release(sm->mcalc);
  gkyl_array_release(sm->marr);

  // free the weak division memory if not computing integrated moments
  if (!sm->is_integrated)
    gkyl_dg_bin_op_mem_release(sm->mem_geo);
}
