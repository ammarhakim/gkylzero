#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_radiation_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_rad_drag *rad)
{
  // allocate drag coefficients in vparallel and mu
  // vnu = 2/pi*|v|*nu(v)
  // vsqnu = 1/2*(m/B)^(3/2)*sqrt(mu)*|v|^2*nu(v)
  // where |v| = sqrt(v_par^2 + 2 mu B/m)
  // Note that through the spatial variation of B, both these drag coefficients depend on the full phase space
  rad->vnu = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  rad->vsqnu = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);

  rad->vnu_host = rad->vnu;
  rad->vsqnu_host = rad->vsqnu;
  if (app->use_gpu) {
    rad->vnu_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);
    rad->vsqnu_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);
  }

  // initialize drag coefficients
  rad->calc_gk_rad_vars = gkyl_dg_calc_gk_rad_vars_new(&s->grid, &app->confBasis, &app->basis, 
    s->info.charge, s->info.mass, app->gk_geom);
  gkyl_dg_calc_gk_rad_vars_advance(rad->calc_gk_rad_vars, &app->local, &s->local, rad->vnu_host, rad->vsqnu_host);
  if (app->use_gpu) {
    gkyl_array_copy(rad->vnu, rad->vnu_host);
    gkyl_array_copy(rad->vsqnu_host, rad->vsqnu);
  }

  // allocate moments needed for radiation update
  gk_species_moment_init(app, s, &rad->moms, "M0");

  // Radiation updater
  struct gkyl_dg_rad_gyrokinetic_drag_auxfields drag_inp = { .nI = rad->moms.marr, 
    .vnu = rad->vnu, .vsqnu = rad->vsqnu };
  rad->drag_slvr = gkyl_dg_updater_rad_gyrokinetic_new(&s->grid, 
    &app->confBasis, &app->basis, &app->local, &s->local, &drag_inp, app->use_gpu);
}

// updates the collision terms in the rhs
void
gk_species_radiation_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_rad_drag *rad, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
    
  // accumulate update due to collisions onto rhs
  gkyl_dg_updater_rad_gyrokinetic_advance(rad->drag_slvr, &species->local,
    fin, species->cflrate, rhs);
  
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void 
gk_species_radiation_release(const struct gkyl_gyrokinetic_app *app, const struct gk_rad_drag *rad)
{
  gkyl_array_release(rad->vnu);
  gkyl_array_release(rad->vsqnu);
  gkyl_dg_calc_gk_rad_vars_release(rad->calc_gk_rad_vars);

  gk_species_moment_release(app, &rad->moms);

  gkyl_dg_updater_rad_gyrokinetic_release(rad->drag_slvr);
 }
