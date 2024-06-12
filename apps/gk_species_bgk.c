#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_const.h>

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
    double eps0 = s->info.collisions.eps0 ? s->info.collisions.eps0: GKYL_EPSILON0;
    double hbar = s->info.collisions.hbar ? s->info.collisions.hbar: GKYL_PLANCKS_CONSTANT_H/2/M_PI;
    double eV = s->info.collisions.eV ? s->info.collisions.eV: GKYL_ELEMENTARY_CHARGE;

    struct gkyl_array* bmag_mid_host = app->use_gpu? mkarr(false, 1, 1) : gkyl_array_acquire(app->gk_geom->bmag_mid);
    gkyl_array_copy(bmag_mid_host, app->gk_geom->bmag_mid);
    double *bmag_mid_ptr = gkyl_array_fetch(bmag_mid_host, 0);
    double bmag_mid = s->info.collisions.bmag_mid ? s->info.collisions.bmag_mid : bmag_mid_ptr[0];
    gkyl_array_release(bmag_mid_host);

    // Compute a minimum representable temperature based on the smallest dv in the grid.
    double dv_min[vdim];
    gkyl_velocity_map_reduce_dv_range(s->vel_map, GKYL_MIN, dv_min, s->vel_map->local_vel);

    double tpar_min = (s->info.mass/6.0)*pow(dv_min[0],2);
    double tperp_min = vdim>1 ? (bmag_mid/3.0)*dv_min[1] : tpar_min;
    bgk->vtsq_min = (tpar_min + 2.0*tperp_min)/(3.0*s->info.mass);

    bgk->spitzer_calc = gkyl_spitzer_coll_freq_new(&app->confBasis, app->poly_order+1,
      nuFrac, 1.0, 1.0, app->use_gpu);
    bgk->self_nu_fac = nuFrac*gkyl_calc_norm_nu(s->info.collisions.n_ref, s->info.collisions.n_ref, 
      s->info.mass, s->info.mass, s->info.charge, s->info.charge, 
      s->info.collisions.T_ref, s->info.collisions.T_ref, bmag_mid, eps0, hbar, eV);

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
  // Density and T/m for Spitzer nu
  bgk->m0 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  bgk->vtsq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  // allocate moments needed for BGK collisions update
  gk_species_moment_init(app, s, &bgk->moms, "MaxwellianMoments");

  // Maxwellian correction updater
  bgk->correct_all_moms = false;
  int max_iter = s->info.collisions.max_iter > 0 ? s->info.collisions.max_iter : 50;
  double iter_eps = s->info.collisions.iter_eps > 0 ? s->info.collisions.iter_eps  : 1e-10;
  bool use_last_converged = s->info.collisions.use_last_converged;

  struct gkyl_gyrokinetic_maxwellian_correct_inp inp_corr = {
    .phase_grid = &s->grid,
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .conf_range =  &app->local,
    .conf_range_ext = &app->local_ext,
    .vel_range = &s->local_vel,
    .gk_geom = app->gk_geom,
    .vel_map = s->vel_map,
    .divide_jacobgeo = false, // final Jacobian multiplication will be handled in advance
    .use_last_converged = use_last_converged, // flag for if we utilizing the results of the scheme 
                                              // *even if* it doesn't converge
    .mass = s->info.mass,
    .use_gpu = app->use_gpu,
    .max_iter = max_iter,
    .eps = iter_eps,
  };
  bgk->self_niter = 0;
  bgk->cross_niter = 0;
  bgk->corr_max = gkyl_gyrokinetic_maxwellian_correct_inew( &inp_corr );
  if (s->info.collisions.correct_all_moms) {
    bgk->correct_all_moms = true;
    bgk->corr_stat = gkyl_dynvec_new(GKYL_DOUBLE, 5);
    bgk->is_first_corr_status_write_call = true;
  }

  // Maxwellian projection updater.
  bgk->proj_max = gkyl_proj_maxwellian_on_basis_new(&s->grid, &app->confBasis, &app->basis, 
    app->poly_order+1, s->vel_map, app->use_gpu);
  bgk->fmax = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  bgk->nu_fmax = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  // BGK updater (also computes stable timestep)
  bgk->up_bgk = gkyl_bgk_collisions_new(&app->confBasis, &app->basis, app->use_gpu);
}

void 
gk_species_bgk_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_bgk_collisions *bgk)
{  
  // set pointers to species we cross-collide with
  for (int i=0; i<bgk->num_cross_collisions; ++i) {
    bgk->collide_with[i] = gk_find_species(app, s->info.collisions.collide_with[i]);
    if (s->info.collisions.normNu) {
      double nuFrac = s->info.collisions.nuFrac ? s->info.collisions.nuFrac : 1.0;
      double eps0 = s->info.collisions.eps0 ? s->info.collisions.eps0: GKYL_EPSILON0;
      double hbar = s->info.collisions.hbar ? s->info.collisions.hbar: GKYL_PLANCKS_CONSTANT_H/2/M_PI;
      double eV = s->info.collisions.eV ? s->info.collisions.eV: GKYL_ELEMENTARY_CHARGE;
      struct gkyl_array* bmag_mid_host = app->gk_geom->bmag_mid;
      if (app->use_gpu) {
        bmag_mid_host = mkarr(false, 1, 1);
        gkyl_array_copy(bmag_mid_host, app->gk_geom->bmag_mid);
      }
      double *bmag_mid_ptr = gkyl_array_fetch(bmag_mid_host, 0);
      double bmag_mid = s->info.collisions.bmag_mid ? s->info.collisions.bmag_mid : bmag_mid_ptr[0];
      if (app->use_gpu)
        gkyl_array_release(bmag_mid_host);
      bgk->cross_nu_fac[i] = nuFrac*gkyl_calc_norm_nu(s->info.collisions.n_ref, bgk->collide_with[i]->info.collisions.n_ref, 
        s->info.mass, bgk->collide_with[i]->info.mass, s->info.charge, bgk->collide_with[i]->info.charge, 
        s->info.collisions.T_ref, bgk->collide_with[i]->info.collisions.T_ref, bmag_mid, eps0, hbar, eV);
    }    
    bgk->other_m[i] = bgk->collide_with[i]->info.mass;
    bgk->other_moms[i] = bgk->collide_with[i]->bgk.moms.marr;
    bgk->other_nu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    bgk->cross_nu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    
    if (bgk->other_m[i] > s->info.mass) {
      gkyl_array_set(bgk->cross_nu[i], sqrt(2.), bgk->self_nu);
      gkyl_array_set(bgk->other_nu[i], sqrt(2.)*(s->info.mass)/(bgk->other_m[i]), bgk->self_nu);
    } else {
      gkyl_array_set(bgk->cross_nu[i], sqrt(2.)*(bgk->other_m[i])/(s->info.mass), bgk->collide_with[i]->bgk.self_nu);
      gkyl_array_set(bgk->other_nu[i], sqrt(2.), bgk->collide_with[i]->bgk.self_nu);
    }
    gkyl_array_accumulate(bgk->nu_sum, 1.0, bgk->cross_nu[i]);

    bgk->cross_moms[i] = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
  }
    
  bgk->betaGreenep1 = 1.0;

  bgk->cross_bgk = gkyl_gyrokinetic_cross_prim_moms_bgk_new(&app->basis, &app->confBasis, app->use_gpu);
}

// computes moments, boundary corrections, and primitive moments
void
gk_species_bgk_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  // compute needed Maxwellian moments (n, u_par, T/m) (Jacobian factors already eliminated)
  gk_species_moment_calc(&bgk->moms, species->local, app->local, fin);

  // Calculate self_nu if using spitzer nu
  if (bgk->normNu) {
    gkyl_array_clear(bgk->nu_sum, 0.0);
    // Fetch n and T/m from Maxwellian moments computed for BGK update
    gkyl_array_set_offset(bgk->m0, 1.0, bgk->moms.marr, 0*app->confBasis.num_basis);
    gkyl_array_set_offset(bgk->vtsq, 1.0, bgk->moms.marr, 2*app->confBasis.num_basis);
    gkyl_spitzer_coll_freq_advance_normnu(bgk->spitzer_calc, &app->local, 
      bgk->vtsq, bgk->vtsq_min, bgk->m0, bgk->vtsq, bgk->vtsq_min, bgk->self_nu_fac, bgk->self_nu);
    gkyl_array_accumulate(bgk->nu_sum, 1.0, bgk->self_nu);
  }
  
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
    // Calculate cross_nu if using spitzer nu
    if (bgk->normNu) {
      gkyl_spitzer_coll_freq_advance_normnu(bgk->spitzer_calc, &app->local, 
        bgk->vtsq, bgk->vtsq_min, bgk->collide_with[i]->bgk.m0, 
        bgk->collide_with[i]->bgk.vtsq, bgk->collide_with[i]->bgk.vtsq_min, 
        bgk->cross_nu_fac[i], bgk->cross_nu[i]);
      gkyl_array_set(bgk->other_nu[i], (species->info.mass)/(bgk->other_m[i]), bgk->cross_nu[i]);
      gkyl_array_accumulate(bgk->nu_sum, 1.0, bgk->cross_nu[i]);
    }

    gkyl_gyrokinetic_cross_prim_moms_bgk_advance(bgk->cross_bgk, &app->local, bgk->betaGreenep1, 
      species->info.mass, bgk->moms.marr, bgk->other_m[i], bgk->other_moms[i], 
      bgk->cross_nu[i], bgk->other_nu[i], bgk->cross_moms[i]);

  }
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

// updates the collision terms in the rhs
void
gk_species_bgk_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_bgk_collisions *bgk, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute the self-collisions Maxwellian.
  gkyl_proj_gkmaxwellian_on_basis_prim_mom(bgk->proj_max, &species->local, &app->local, bgk->moms.marr,
    app->gk_geom->bmag, app->gk_geom->bmag, species->info.mass, bgk->fmax);
  gkyl_array_scale_by_cell(bgk->fmax, species->vel_map->jacobvel); // Multiply by the velocity space jacobian.
  // First correct the density
  gkyl_gyrokinetic_maxwellian_correct_density_moment(bgk->corr_max, 
    bgk->fmax, bgk->moms.marr, &species->local, &app->local);
  // Correct all the moments of the projected gyrokinetic maxwellian distribution function.
  if (bgk->correct_all_moms) {
    struct gkyl_gyrokinetic_maxwellian_correct_status status_corr;
    status_corr = gkyl_gyrokinetic_maxwellian_correct_all_moments(bgk->corr_max, 
      bgk->fmax, bgk->moms.marr, &species->local, &app->local);
    double corr_vec[5] = { 0.0 };
    corr_vec[0] = status_corr.num_iter;
    corr_vec[1] = status_corr.iter_converged;
    corr_vec[2] = status_corr.error[0];
    corr_vec[3] = status_corr.error[1];
    corr_vec[4] = status_corr.error[2];

    double corr_vec_global[5] = { 0.0 };
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 5, corr_vec, corr_vec_global);    
    gkyl_dynvec_append(bgk->corr_stat, app->tcurr, corr_vec_global);

    bgk->self_niter += status_corr.num_iter;
  } 
  // Multiply the Maxwellian by the configuration-space Jacobian.
  gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, bgk->fmax, 
    app->gk_geom->jacobgeo, bgk->fmax, &app->local, &species->local);

  // Obtain and accumulate the self-collisions nu*fmax
  gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, bgk->fmax, 
    bgk->self_nu, bgk->fmax, &app->local, &species->local);
  gkyl_array_set(bgk->nu_fmax, 1.0, bgk->fmax);

  // Cross-collisions nu*fmax.
  for (int i=0; i<bgk->num_cross_collisions; ++i) {
    // Compute the Maxwellian.
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(bgk->proj_max, &species->local_ext, &app->local_ext, bgk->cross_moms[i],
      app->gk_geom->bmag, app->gk_geom->bmag, species->info.mass, bgk->fmax);
    gkyl_array_scale_by_cell(bgk->fmax, species->vel_map->jacobvel); // Multiply by the velocity space jacobian.
    // First correct the density
    gkyl_gyrokinetic_maxwellian_correct_density_moment(bgk->corr_max, 
      bgk->fmax, bgk->cross_moms[i], &species->local, &app->local);
    // Correct all the moments of the projected gyrokinetic maxwellian distribution function.
    if (bgk->correct_all_moms) {
      struct gkyl_gyrokinetic_maxwellian_correct_status status_corr;
      status_corr = gkyl_gyrokinetic_maxwellian_correct_all_moments(bgk->corr_max, 
        bgk->fmax, bgk->cross_moms[i], &species->local, &app->local);

      bgk->cross_niter += status_corr.num_iter;
    } 
    // Multiply the Maxwellian by the configuration-space Jacobian.
    gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, bgk->fmax, 
      app->gk_geom->jacobgeo, bgk->fmax, &app->local, &species->local);

    // Compute and accumulate nu*fmax.
    gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, bgk->fmax, 
      bgk->cross_nu[i], bgk->fmax, &app->local, &species->local);
    gkyl_array_accumulate(bgk->nu_fmax, 1.0, bgk->fmax);
  }

  gkyl_bgk_collisions_advance(bgk->up_bgk, &app->local, &species->local, 
    bgk->nu_sum, bgk->nu_fmax, fin, rhs, species->cflrate);

  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void 
gk_species_bgk_release(const struct gkyl_gyrokinetic_app *app, const struct gk_bgk_collisions *bgk)
{
  gkyl_array_release(bgk->self_nu);
  gkyl_array_release(bgk->nu_sum);
  gkyl_array_release(bgk->m0);
  gkyl_array_release(bgk->vtsq);

  gkyl_array_release(bgk->fmax);
  gkyl_array_release(bgk->nu_fmax);

  if (app->use_gpu) 
    gkyl_array_release(bgk->nu_sum_host);

  gk_species_moment_release(app, &bgk->moms);

  if (bgk->normNu) {
    gkyl_array_release(bgk->norm_nu);
    gkyl_array_release(bgk->nu_init);
    gkyl_spitzer_coll_freq_release(bgk->spitzer_calc);
  }

  if (bgk->num_cross_collisions) {
    for (int i=0; i<bgk->num_cross_collisions; ++i) {
      gkyl_array_release(bgk->cross_nu[i]);
      gkyl_array_release(bgk->other_nu[i]);
      gkyl_array_release(bgk->cross_moms[i]);
    }
    gkyl_gyrokinetic_cross_prim_moms_bgk_release(bgk->cross_bgk);
  }
  gkyl_gyrokinetic_maxwellian_correct_release(bgk->corr_max);
  if (bgk->correct_all_moms) {
    gkyl_dynvec_release(bgk->corr_stat);
  }  
  gkyl_proj_maxwellian_on_basis_release(bgk->proj_max);
  gkyl_bgk_collisions_release(bgk->up_bgk);
 }
