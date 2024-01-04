#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_radiation_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_rad_drag *rad)
{
  // Fitting parameters
  double a, alpha, beta, gamma, v0;
  rad->num_cross_collisions = s->info.radiation.num_cross_collisions;
  // initialize drag coefficients
  for (int i=0; i<rad->num_cross_collisions; ++i) {
    // allocate drag coefficients in vparallel and mu for each collision
    // vnu = 2/pi*|v|*nu(v)
    // vsqnu = 1/2*(m/B)^(3/2)*sqrt(mu)*|v|^2*nu(v)
    // where |v| = sqrt(v_par^2 + 2 mu B/m)
    // Note that through the spatial variation of B, both these drag coefficients depend on the full phase space
    rad->vnu[i] = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
    rad->vsqnu[i] = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);

    rad->vnu_host[i] = rad->vnu[i];
    rad->vsqnu_host[i] = rad->vsqnu[i];
    if (app->use_gpu) {
      rad->vnu_host[i] = mkarr(false, app->basis.num_basis, s->local_ext.volume);
      rad->vsqnu_host[i] = mkarr(false, app->basis.num_basis, s->local_ext.volume);
    }

    // Fetch the species we are colliding with and the fitting parameters for that species
    rad->collide_with[i] = gk_find_species(app, s->info.radiation.collide_with[i]);
    rad->collide_with_idx[i] = gk_find_species_idx(app, s->info.radiation.collide_with[i]);
    a = s->info.radiation.a[i];
    alpha = s->info.radiation.alpha[i];
    beta = s->info.radiation.beta[i];
    gamma = s->info.radiation.gamma[i];
    v0 = s->info.radiation.v0[i];

    rad->calc_gk_rad_vars[i] = gkyl_dg_calc_gk_rad_vars_new(&s->grid, &app->confBasis, &app->basis, 
      s->info.charge, s->info.mass, app->gk_geom, 
      a, alpha, beta, gamma, v0);

    gkyl_dg_calc_gk_rad_vars_advance(rad->calc_gk_rad_vars[i], &app->local, &s->local, rad->vnu_host[i], rad->vsqnu_host[i]);
    if (app->use_gpu) {
      gkyl_array_copy(rad->vnu[i], rad->vnu_host[i]);
      gkyl_array_copy(rad->vsqnu[i], rad->vsqnu_host[i]);
    }
    // allocate density calculation needed for radiation update
    gk_species_moment_init(app, rad->collide_with[i], &rad->moms[i], "M0");
  }

  // Radiation updater
  struct gkyl_dg_rad_gyrokinetic_drag_auxfields drag_inp = { .nI = rad->moms[0].marr, 
    .vnu = rad->vnu[0], .vsqnu = rad->vsqnu[0] };
  rad->drag_slvr = gkyl_dg_updater_rad_gyrokinetic_new(&s->grid, 
    &app->confBasis, &app->basis, &app->local, &s->local, &drag_inp, app->use_gpu);
}

// // computes density for computation of total radiation drag
// void
// gk_species_radiation_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
//   struct gk_rad_drag *rad, const struct gkyl_array *fin)
// {
//   gkyl_array_clear(rad->nvnu_sum, 0.0);
//   gkyl_array_clear(rad->nvsqnu_sum, 0.0);
//   for (int i=0; i<rad->num_cross_collisions; ++i) {
//     // compute needed moments
//     gk_species_moment_calc(&rad->moms[i], species->local, app->local, fin[collide_with_idx[i]]);

//     gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, rad->nvnu[i], rad->moms.marr, rad->vnu[i],
//       &app->local, &species->local);
//     gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, rad->nvsqnu[i], rad->moms.marr, rad->vsqnu[i],
//       &app->local, &species->local);

//     gkyl_array_accumulate(rad->nvnu_sum, 1.0, rad->nvnu[i]);
//     gkyl_array_accumulate(rad->nvsqnu_sum, 1.0, rad->nvsqnu[i]);
//   }
// }

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
  for (int i=0; i<rad->num_cross_collisions; ++i) {
    gkyl_array_release(rad->vnu[i]);
    gkyl_array_release(rad->vsqnu[i]);
    if (app->use_gpu) {
      gkyl_array_release(rad->vnu_host[i]);
      gkyl_array_release(rad->vsqnu_host[i]);      
    }
    gkyl_dg_calc_gk_rad_vars_release(rad->calc_gk_rad_vars[i]);

    gk_species_moment_release(app, &rad->moms[i]);
  }
  gkyl_dg_updater_rad_gyrokinetic_release(rad->drag_slvr);
 }
