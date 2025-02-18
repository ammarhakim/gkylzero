#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

// initialize species moment object
void
gk_species_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_species_moment *sm, const char *nm)
{
  assert(is_moment_name_valid(nm));

  sm->is_integrated = strcmp(nm, "Integrated") == 0;

  sm->is_maxwellian_moms = strcmp("MaxwellianMoments", nm) == 0;
  sm->is_bimaxwellian_moms = strcmp("BiMaxwellianMoments", nm) == 0;
  if (sm->is_maxwellian_moms || sm->is_bimaxwellian_moms) {
    struct gkyl_gk_maxwellian_moments_inp inp_mom = {
      .phase_grid = &s->grid,
      .conf_basis = &app->basis,
      .phase_basis = &s->basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .gk_geom = app->gk_geom,
      .vel_map = s->vel_map,
      .divide_jacobgeo = false, 
      .mass = s->info.mass,
      .use_gpu = app->use_gpu,
    };
    sm->gyrokinetic_maxwellian_moms = gkyl_gk_maxwellian_moments_inew(  &inp_mom  );
    if (sm->is_maxwellian_moms) {
      sm->num_mom = 3; // (n, u_par, T/m)
    }
    else {
      sm->num_mom = 4; // (n, u_par, T_par/m, T_perp/m)
    }  
  }
  else {
    sm->mcalc = gkyl_dg_updater_moment_gyrokinetic_new(&s->grid, &app->basis, 
      &s->basis, &app->local, s->info.mass, s->vel_map, app->gk_geom,
      nm, sm->is_integrated, app->use_gpu);    

    sm->num_mom = gkyl_dg_updater_moment_gyrokinetic_num_mom(sm->mcalc);
  }

  if (sm->is_integrated) {
    sm->marr = mkarr(app->use_gpu, sm->num_mom, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu) {
      sm->marr_host = mkarr(false, sm->num_mom, app->local_ext.volume);  
    }
  }
  else {
    sm->marr = mkarr(app->use_gpu, sm->num_mom*app->basis.num_basis, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu) {
      sm->marr_host = mkarr(false, sm->num_mom*app->basis.num_basis, app->local_ext.volume);
    }
    // Bin Op memory for rescaling moment by inverse of Jacobian
    if (app->use_gpu) {
      sm->mem_geo = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis);
    }
    else {
      sm->mem_geo = gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);
    }
  }
}

void
gk_species_moment_calc(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  if (sm->is_maxwellian_moms) {
    gkyl_gk_maxwellian_moments_advance(sm->gyrokinetic_maxwellian_moms, 
      &phase_rng, &conf_rng, fin, sm->marr);
  }
  else if (sm->is_bimaxwellian_moms) {
    gkyl_gk_bimaxwellian_moments_advance(sm->gyrokinetic_maxwellian_moms, 
      &phase_rng, &conf_rng, fin, sm->marr);
  } 
  else {
    gkyl_dg_updater_moment_gyrokinetic_advance(sm->mcalc, 
      &phase_rng, &conf_rng, fin, sm->marr);
  }
}

// release memory for moment data object
void
gk_species_moment_release(const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm)
{
  if (app->use_gpu) {
    gkyl_array_release(sm->marr_host);
  }

  if (sm->is_maxwellian_moms || sm->is_bimaxwellian_moms) {
    gkyl_gk_maxwellian_moments_release(sm->gyrokinetic_maxwellian_moms);
  }
  else {
    gkyl_dg_updater_moment_gyrokinetic_release(sm->mcalc);
  }
  gkyl_array_release(sm->marr);

  // free the weak division memory if not computing integrated moments
  if (!sm->is_integrated) {
    gkyl_dg_bin_op_mem_release(sm->mem_geo);
  }
}
