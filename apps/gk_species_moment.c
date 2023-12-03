#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

// initialize species moment object
void
gk_species_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_species_moment *sm, const char *nm)
{
  assert(is_moment_name_valid(nm));

  bool is_integrated = strcmp(nm, "Integrated") == 0;

  // Acquire a pointer to the geometry for future rescaling of moments by Jacobian
  sm->gk_geom = gkyl_gk_geometry_acquire(app->gk_geom);
  sm->mem_geo = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);

  sm->mcalc = gkyl_dg_updater_moment_gyrokinetic_new(&s->grid, &app->confBasis, 
    &app->basis, &app->local, &s->local_vel, s->info.mass, sm->gk_geom,
    nm, is_integrated, app->use_gpu);    

  int num_mom = gkyl_dg_updater_moment_gyrokinetic_num_mom(sm->mcalc);

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
gk_species_moment_calc(struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  gkyl_dg_updater_moment_gyrokinetic_advance(sm->mcalc, &phase_rng, &conf_rng, fin, sm->marr);
  // Rescale moment by inverse of Jacobian
  gkyl_dg_div_op_range(sm->mem_geo, app->confBasis, 0, sm->marr, 0, sm->marr, 0, sm->gk_geom->jacobgeo, &app->local);
}

// release memory for moment data object
void
gk_species_moment_release(const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm)
{
  gkyl_gk_geometry_release(sm->gk_geom);
  gkyl_dg_bin_op_mem_release(sm->mem_geo);

  if (app->use_gpu)
    gkyl_array_release(sm->marr_host);

  gkyl_dg_updater_moment_gyrokinetic_release(sm->mcalc);
  gkyl_array_release(sm->marr);
}
