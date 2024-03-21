#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_bgk_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_bgk_collisions *bgk)
{
  int cdim = app->cdim, vdim = app->vdim;
  // allocate nu and initialize it
  bgk->nu_sum = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  bgk->self_nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array *self_nu = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  
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

  // Host-side copy for I/O
  bgk->nu_sum_host = bgk->nu_sum;
  if (app->use_gpu) {
    bgk->nu_sum_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  }

  bgk->model_id = s->model_id;
  // allocate moments needed for BGK collisions update
  vm_species_moment_init(app, s, &bgk->moms, "LTEMoments");

  if (bgk->model_id == GKYL_MODEL_SR) {
    bgk->proj_mj = gkyl_proj_mj_on_basis_new(&s->grid, &app->confBasis, &app->basis, 
      app->poly_order+1, app->use_gpu);
  }
  else {
    bgk->proj_max = gkyl_proj_maxwellian_on_basis_new(&s->grid, &app->confBasis, &app->basis, 
      app->poly_order+1, app->use_gpu);
  }

  bgk->correct_all_moms = false; 
  if (s->info.collisions.correct_all_moms) {
    bgk->correct_all_moms = true;
  }
  bgk->corr_lte = gkyl_correct_vlasov_lte_new(&s->grid, &app->confBasis, &app->basis, 
    &app->local, &app->local_ext, &s->local_vel, 
    s->p_over_gamma, s->gamma, s->gamma_inv, 
    bgk->model_id, s->info.mass, app->use_gpu);

  bgk->f_lte = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  bgk->nu_f_lte = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  // BGK updater (also computes stable timestep)
  bgk->up_bgk = gkyl_bgk_collisions_new(&app->confBasis, &app->basis, app->use_gpu);
}

// computes moments
void
vm_species_bgk_moms(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  vm_species_moment_calc(&bgk->moms, species->local, app->local, fin);
  
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

// updates the collision terms in the rhs
void
vm_species_bgk_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(bgk->nu_f_lte, 0.0);

  // Obtain self-collisions nu*f_lte where f_lte is the LTE distribution 
  // e.g., Maxwellian for non-relativistic and Maxwell-Juttner for relativistic
  if (bgk->model_id == GKYL_MODEL_SR) {
    gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(bgk->proj_mj, &species->local, &app->local,
      bgk->moms.marr, bgk->f_lte);
  }
  else { 
    gkyl_proj_maxwellian_on_basis_prim_mom(bgk->proj_max, &species->local, &app->local, 
      bgk->moms.marr, bgk->f_lte);
  }
  // Correct the density of the projected LTE distribution function through rescaling
  gkyl_correct_density_moment_vlasov_lte(bgk->corr_lte, bgk->f_lte, bgk->moms.marr,
    &species->local, &app->local);
  // Correct all the moments of the projected LTE distribution function 
  if (bgk->correct_all_moms) {
    gkyl_correct_all_moments_vlasov_lte(bgk->corr_lte, bgk->f_lte, bgk->moms.marr,
      &species->local, &app->local);
  } 

  gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, bgk->f_lte, 
    bgk->self_nu, bgk->f_lte, &app->local, &species->local);
  gkyl_array_accumulate(bgk->nu_f_lte, 1.0, bgk->f_lte);

  gkyl_bgk_collisions_advance(bgk->up_bgk, &app->local, &species->local, 
    bgk->nu_sum, bgk->nu_f_lte, fin, rhs, species->cflrate);

  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void 
vm_species_bgk_release(const struct gkyl_vlasov_app *app, const struct vm_bgk_collisions *bgk)
{
  gkyl_array_release(bgk->self_nu);
  gkyl_array_release(bgk->nu_sum);

  gkyl_array_release(bgk->f_lte);
  gkyl_array_release(bgk->nu_f_lte);

  if (app->use_gpu) {
    gkyl_array_release(bgk->nu_sum_host);
  }

  if (bgk->normNu) {
    gkyl_array_release(bgk->norm_nu);
    gkyl_array_release(bgk->nu_init);
    gkyl_spitzer_coll_freq_release(bgk->spitzer_calc);
  }

  vm_species_moment_release(app, &bgk->moms);

  if (bgk->model_id == GKYL_MODEL_SR) {
    gkyl_proj_mj_on_basis_release(bgk->proj_mj);
  } 
  else {
    gkyl_proj_maxwellian_on_basis_release(bgk->proj_max);
  }
  gkyl_correct_vlasov_lte_release(bgk->corr_lte);

  gkyl_bgk_collisions_release(bgk->up_bgk);
}
