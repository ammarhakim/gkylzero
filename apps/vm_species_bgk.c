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
  if (app->use_gpu) 
    bgk->nu_sum_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);

  if (s->model_id == GKYL_MODEL_SR) {
    bgk->proj_mj = gkyl_proj_mj_on_basis_new(&s->grid, &app->confBasis, &app->basis,app->poly_order+1);
    bgk->corr_mj = gkyl_correct_mj_new(&s->grid, &app->confBasis,
      &app->basis, &app->local, &app->local_ext, &s->local_vel, s->p_over_gamma, s->gamma, s->gamma_inv, false);
        // Correct the distribution function
    bgk->mj_moms = gkyl_mj_moments_new(&s->grid, &app->confBasis,
      &app->basis, &app->local, &s->local_vel, app->local.volume, app->local_ext.volume, 
      s->p_over_gamma, s->gamma, s->gamma_inv, false);
  } else{
     // allocate moments needed for BGK collisions update
    vm_species_moment_init(app, s, &bgk->moms, "FiveMoments");
    bgk->proj_max = gkyl_proj_maxwellian_on_basis_new(&s->grid, &app->confBasis, &app->basis, 
      app->poly_order+1, app->use_gpu);
  }
  bgk->fmax = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  bgk->nu_fmax = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  // BGK updater (also computes stable timestep)
  bgk->up_bgk = gkyl_bgk_collisions_new(&app->confBasis, &app->basis, app->use_gpu);
}

// computes moments
void
vm_species_bgk_moms(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  // compute needed moments
  vm_species_moment_calc(&bgk->moms, species->local, app->local, fin);
  
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

// updates the collision terms in the rhs
void
vm_species_bgk_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(bgk->nu_fmax, 0.0);

  // Obtain self-collisions nu*fmax
  if (species->model_id == GKYL_MODEL_SR) {

    // Build memory for custom SR moments
    struct gkyl_array *m0, *m1i, *m2;
    m0 = gkyl_array_new(GKYL_DOUBLE, app->confBasis.num_basis,   app->local_ext.volume);
    m1i = gkyl_array_new(GKYL_DOUBLE, app->vdim *  app->confBasis.num_basis,  app->local_ext.volume);
    m2 = gkyl_array_new(GKYL_DOUBLE, app->confBasis.num_basis,   app->local_ext.volume);

    // Compute MJ via the moments
    gkyl_mj_moments_advance(bgk->mj_moms, fin, m0, m1i, m2, &species->local, &app->local);
    gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(bgk->proj_mj, &species->local, &app->local,
      m0, m1i, m2, bgk->fmax);
    gkyl_correct_mj_fix_m0(bgk->corr_mj, bgk->fmax, m0, m1i,
      &species->local, &app->local);
    //gkyl_correct_mj_fix(bgk->corr_mj, bgk->fmax, m0, m1i, m2,
     // &species->local, &app->local, app->poly_order);


    //printf("\n%%------- Exact Moments 1x1v bgk ------\n");
    //struct gkyl_range_iter biter;
    //gkyl_range_iter_init(&biter, &app->local);
    //while (gkyl_range_iter_next(&biter)) {
    //    long midx = gkyl_range_idx(&app->local, biter.idx);
    //    if (midx == 1){
    //      const double *gamma_val = gkyl_array_cfetch(&species->gamma, midx);
    //      const double *num = gkyl_array_cfetch(m0, midx);
    //      const double *vb = gkyl_array_cfetch(m1i, midx);
    //      const double *T = gkyl_array_cfetch(m2, midx);
    //      printf("n = %1.16e;\n",num[0]/sqrt(2));
    //      printf("v = %1.16e;\n",vb[0]/sqrt(2));
    //      printf("T = %1.16e;\n",T[0]/sqrt(2));
    //    }
    //}

    // Release the moment-memory
    gkyl_array_release(m0);
    gkyl_array_release(m1i);
    gkyl_array_release(m2);
  }else{ 
    gkyl_proj_maxwellian_on_basis_lab_mom(bgk->proj_max, &species->local, &app->local, 
      bgk->moms.marr, bgk->fmax);
  }
  gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, bgk->fmax, 
    bgk->self_nu, bgk->fmax, &app->local, &species->local);
  gkyl_array_accumulate(bgk->nu_fmax, 1.0, bgk->fmax);

  gkyl_bgk_collisions_advance(bgk->up_bgk, &app->local, &species->local, 
    bgk->nu_sum, bgk->nu_fmax, fin, rhs, species->cflrate);

  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void 
vm_species_bgk_release(const struct gkyl_vlasov_app *app, const struct vm_bgk_collisions *bgk)
{
  gkyl_array_release(bgk->self_nu);
  gkyl_array_release(bgk->nu_sum);

  gkyl_array_release(bgk->fmax);
  gkyl_array_release(bgk->nu_fmax);

  if (app->use_gpu) {
    gkyl_array_release(bgk->nu_sum_host);
  }

  vm_species_moment_release(app, &bgk->moms);

  if (bgk->normNu) {
    gkyl_array_release(bgk->norm_nu);
    gkyl_array_release(bgk->nu_init);
    gkyl_spitzer_coll_freq_release(bgk->spitzer_calc);
  }

  gkyl_proj_maxwellian_on_basis_release(bgk->proj_max);
  gkyl_proj_mj_on_basis_release(bgk->proj_mj);
  gkyl_correct_mj_release(bgk->corr_mj);
  gkyl_mj_moments_release(bgk->mj_moms);
  gkyl_bgk_collisions_release(bgk->up_bgk);
 }
