#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_bgk_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_bgk_collisions *bgk)
{
  int cdim = app->cdim, vdim = app->vdim;
  // allocate nu and initialize it
  bgk->nu_sum = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  bgk->self_nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array *self_nu = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);

  bgk->num_cross_collisions = s->info.collisions.num_cross_collisions;
  
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    app->poly_order+1, 1, s->info.collisions.self_nu, s->info.collisions.ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local, self_nu);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(bgk->self_nu, self_nu);
  gkyl_array_copy(bgk->nu_sum, self_nu);
  gkyl_array_release(self_nu);

  bgk->spitzer_calc = 0;
  bgk->normNu = false;
  if (s->info.collisions.normNu) {
    bgk->normNu = true;
    double nuFrac = s->info.collisions.nuFrac ? s->info.collisions.nuFrac : 1.0;
    double eps0 = 1.0;
    double hbar = 1.0;
    bgk->spitzer_calc = gkyl_spitzer_coll_freq_new(&app->confBasis, app->poly_order+1,
      nuFrac, eps0, hbar, app->use_gpu);
    // Create arrays for scaling collisionality by normalization factor
    // norm_nu is computed from Spitzer calc and is the normalization factor for the local
    // density and thermal velocity, norm_nu_sr = n/(vth_s^2 + vth_r^2)^(3/2)
    // nu_init is the inital collisionality profile, which must be stored so that at every time
    // time step the collisionality profile is properly scaled and the effects are not cumulative
    bgk->norm_nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    bgk->nu_init = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    gkyl_array_copy(bgk->nu_init, bgk->self_nu);
  }

  // Density for if we have cross collisions to compute Greene's factor
  bgk->m0 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  // Host-side copy for I/O
  bgk->nu_sum_host = bgk->nu_sum;
  if (app->use_gpu) 
    bgk->nu_sum_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);

  // allocate moments needed for BGK collisions update
  gk_species_moment_init(app, s, &bgk->moms, "ThreeMoments");

  // Maxwellian correction updater
  struct gkyl_correct_maxwellian_gyrokinetic_inp inp = {
    .phase_grid = &s->grid,
    .conf_grid = &app->grid,
    .phase_basis = &app->basis,
    .conf_basis = &app->confBasis,

    .conf_local = &app->local,
    .conf_local_ext = &app->local_ext,
    .mass = s->info.mass, 
    .gk_geom = app->gk_geom,
    .max_iter = 50, 
    .eps_err = 1.0e-14, 
    .use_gpu = app->use_gpu
  };
  bgk->corr_max = gkyl_correct_maxwellian_gyrokinetic_new(&inp);  
  bgk->proj_max = gkyl_proj_maxwellian_on_basis_new(&s->grid, &app->confBasis, &app->basis, 
    app->poly_order+1, app->use_gpu);
  bgk->fmax = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  bgk->nu_fmax = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  bgk->nu_sum_f = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
}

void 
gk_species_bgk_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_bgk_collisions *bgk)
{  
  bgk->greene_factor_mem = 0;
  if (app->use_gpu)
    bgk->greene_factor_mem = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
  else
    bgk->greene_factor_mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);

  // set pointers to species we cross-collide with
  for (int i=0; i<bgk->num_cross_collisions; ++i) {
    bgk->collide_with[i] = gk_find_species(app, s->info.collisions.collide_with[i]);
    bgk->other_m[i] = bgk->collide_with[i]->info.mass;
    bgk->other_moms[i] = bgk->collide_with[i]->bgk.moms.marr;
    bgk->other_nu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    bgk->cross_nu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    bgk->greene_num[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    bgk->greene_den[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    bgk->greene_factor[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    
    if (bgk->other_m[i] > s->info.mass) {
      gkyl_array_set(bgk->cross_nu[i], sqrt(2), bgk->self_nu);
      gkyl_array_set(bgk->other_nu[i], sqrt(2)*(s->info.mass)/(bgk->other_m[i]), bgk->self_nu);
    } else {
      gkyl_array_set(bgk->cross_nu[i], sqrt(2)*(bgk->other_m[i])/(s->info.mass), bgk->collide_with[i]->bgk.self_nu);
      gkyl_array_set(bgk->other_nu[i], sqrt(2), bgk->collide_with[i]->bgk.self_nu);
    }
    
    gkyl_array_accumulate(bgk->nu_sum, 1.0, bgk->cross_nu[i]);

    bgk->other_mnu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    bgk->other_mnu_m0[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    bgk->self_mnu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    bgk->self_mnu_m0[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    gkyl_array_set(bgk->self_mnu[i], s->info.mass, bgk->cross_nu[i]);
    gkyl_array_set(bgk->other_mnu[i], bgk->other_m[i], bgk->other_nu[i]);

    bgk->cross_moms[i] = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
  }
    
  bgk->betaGreenep1 = 1.0;

  bgk->cross_bgk = gkyl_mom_cross_bgk_gyrokinetic_new(&app->basis, &app->confBasis, app->use_gpu);
}

// computes moments, boundary corrections, and primitive moments
void
gk_species_bgk_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  // compute needed moments
  gk_species_moment_calc(&bgk->moms, species->local, app->local, fin);
  gkyl_array_set_range(bgk->m0, 1.0, bgk->moms.marr, &app->local);
  
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

// computes moments from cross-species collisions
void
gk_species_bgk_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();
  
  wst = gkyl_wall_clock();  
  for (int i=0; i<bgk->num_cross_collisions; ++i) {
    gkyl_dg_mul_op_range(app->confBasis, 0, bgk->self_mnu_m0[i], 0,
      bgk->self_mnu[i], 0, bgk->m0, &app->local);
    gkyl_dg_mul_op_range(app->confBasis, 0, bgk->other_mnu_m0[i], 0,
      bgk->other_mnu[i], 0, bgk->collide_with[i]->bgk.m0, &app->local);

    gkyl_dg_mul_op_range(app->confBasis, 0, bgk->greene_num[i], 0,
      bgk->other_mnu_m0[i], 0, bgk->m0, &app->local);

    gkyl_array_set(bgk->greene_den[i], 1.0, bgk->self_mnu_m0[i]);
    gkyl_array_accumulate(bgk->greene_den[i], 1.0, bgk->other_mnu_m0[i]);

    gkyl_dg_div_op_range(bgk->greene_factor_mem, app->confBasis, 0, bgk->greene_factor[i], 0,
      bgk->greene_num[i], 0, bgk->greene_den[i], &app->local);
    gkyl_array_scale(bgk->greene_factor[i], 2*bgk->betaGreenep1);

    gkyl_mom_cross_bgk_gyrokinetic_advance(bgk->cross_bgk, &app->local, bgk->betaGreenep1, 
      species->info.mass, bgk->moms.marr, bgk->other_m[i], bgk->other_moms[i], 
      bgk->self_nu, bgk->other_nu[i], bgk->cross_moms[i]);
  }
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

// updates the collision terms in the rhs
void
gk_species_bgk_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(bgk->nu_fmax, 0.0);

  // Obtain self-collisions nu*fmax
  gkyl_proj_gkmaxwellian_on_basis_lab_mom(bgk->proj_max, &species->local_ext, &app->local_ext, bgk->moms.marr,
    app->gk_geom->bmag, app->gk_geom->bmag, species->info.mass, bgk->fmax);
  gkyl_correct_maxwellian_gyrokinetic_advance(bgk->corr_max, bgk->fmax, bgk->moms.marr, &app->local, &species->local);
  gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, bgk->fmax, 
    bgk->self_nu, bgk->fmax, &app->local, &species->local);
  gkyl_array_accumulate(bgk->nu_fmax, 1.0, bgk->fmax);

  // Accumulate cross-collisions nu*fmax
  for (int i=0; i<bgk->num_cross_collisions; ++i) {
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(bgk->proj_max, &species->local_ext, &app->local_ext, bgk->cross_moms[i],
      app->gk_geom->bmag, app->gk_geom->bmag, species->info.mass, bgk->fmax);
    gkyl_correct_maxwellian_gyrokinetic_advance(bgk->corr_max, bgk->fmax, bgk->moms.marr, &app->local, &species->local);
    gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, bgk->fmax, 
      bgk->cross_nu[i], bgk->fmax, &app->local, &species->local);
    gkyl_array_accumulate(bgk->nu_fmax, 1.0, bgk->fmax);
  }

  // Accumulate RHS of BGK collisions
  gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, bgk->nu_sum_f, 
    bgk->nu_sum, fin, &app->local, &species->local);  
  gkyl_array_accumulate(rhs, 1.0, bgk->nu_fmax);
  gkyl_array_accumulate(rhs, -1.0, bgk->nu_sum_f);

  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void 
gk_species_bgk_release(const struct gkyl_gyrokinetic_app *app, const struct gk_bgk_collisions *bgk)
{
  gkyl_array_release(bgk->self_nu);
  gkyl_array_release(bgk->nu_sum);
  gkyl_array_release(bgk->m0);

  gkyl_array_release(bgk->fmax);
  gkyl_array_release(bgk->nu_fmax);
  gkyl_array_release(bgk->nu_sum_f);

  if (app->use_gpu) 
    gkyl_array_release(bgk->nu_sum_host);

  gk_species_moment_release(app, &bgk->moms);

  if (bgk->normNu) {
    gkyl_array_release(bgk->norm_nu);
    gkyl_array_release(bgk->nu_init);
    gkyl_spitzer_coll_freq_release(bgk->spitzer_calc);
  }

  if (bgk->num_cross_collisions) {
    gkyl_dg_bin_op_mem_release(bgk->greene_factor_mem);
    for (int i=0; i<bgk->num_cross_collisions; ++i) {
      gkyl_array_release(bgk->other_moms[i]);
      gkyl_array_release(bgk->cross_nu[i]);
      gkyl_array_release(bgk->other_nu[i]);
      gkyl_array_release(bgk->self_mnu[i]);
      gkyl_array_release(bgk->self_mnu_m0[i]);
      gkyl_array_release(bgk->other_mnu[i]);
      gkyl_array_release(bgk->other_mnu_m0[i]);
      gkyl_array_release(bgk->greene_num[i]);
      gkyl_array_release(bgk->greene_den[i]);
      gkyl_array_release(bgk->greene_factor[i]);
    }
    gkyl_mom_cross_bgk_gyrokinetic_release(bgk->cross_bgk);
  }
  gkyl_correct_maxwellian_gyrokinetic_release(bgk->corr_max);
  gkyl_proj_maxwellian_on_basis_release(bgk->proj_max);
 }