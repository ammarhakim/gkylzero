#include <assert.h>
#include <gkyl_vlasov_poisson_priv.h>
#include <gkyl_const.h>

void 
vp_species_lbo_init(struct gkyl_vlasov_poisson_app *app, struct vp_species *s, struct vp_lbo_collisions *lbo)
{
  int cdim = app->cdim, vdim = app->vdim;
  double v_bounds[2*GKYL_MAX_DIM] = { 0.0 };
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = s->info.lower[d];
    v_bounds[d + vdim] = s->info.upper[d];
  }

  // Allocate nu and initialize it.
  lbo->nu_sum = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  lbo->self_nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array *self_nu = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);

  lbo->num_cross_collisions = s->info.collisions.num_cross_collisions;
  
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    app->poly_order+1, 1, s->info.collisions.self_nu, s->info.collisions.ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local, self_nu);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(lbo->self_nu, self_nu);
  gkyl_array_copy(lbo->nu_sum, self_nu);
  gkyl_array_release(self_nu);

  lbo->spitzer_calc = 0;
  lbo->normNu = false;
  if (s->info.collisions.normNu) {
    lbo->normNu = true;
    double nuFrac = s->info.collisions.nuFrac ? s->info.collisions.nuFrac : 1.0;
    double eps0 = s->info.collisions.eps0 ? s->info.collisions.eps0: GKYL_EPSILON0;
    double hbar = s->info.collisions.hbar ? s->info.collisions.hbar: GKYL_PLANCKS_CONSTANT_H/2/M_PI;
    double eV = s->info.collisions.eV ? s->info.collisions.eV: GKYL_ELEMENTARY_CHARGE;

    double T_min = 0.0;
    for (int d=0; d<vdim; d++) T_min += (s->info.mass/6.0)*pow(s->grid.dx[cdim+d],2);
    T_min = T_min/vdim;
    lbo->vtsq_min = T_min/s->info.mass;

    lbo->spitzer_calc = gkyl_spitzer_coll_freq_new(&app->confBasis, app->poly_order+1,
      nuFrac, 1.0, 1.0, app->use_gpu);
    lbo->self_nu_fac = nuFrac * gkyl_calc_norm_nu(s->info.collisions.n_ref, s->info.collisions.n_ref,
      s->info.mass, s->info.mass, s->info.charge, s->info.charge, s->info.collisions.T_ref,
      s->info.collisions.T_ref, 0.0, eps0, hbar, eV);

    // Create arrays for scaling collisionality by normalization factor
    // norm_nu is computed from Spitzer calc and is the normalization factor for the local
    // density and thermal velocity, norm_nu_sr = n/(vth_s^2 + vth_r^2)^(3/2)
    // nu_init is the inital collisionality profile, which must be stored so that at every time
    // time step the collisionality profile is properly scaled and the effects are not cumulative
    lbo->norm_nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->nu_init = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    gkyl_array_copy(lbo->nu_init, lbo->self_nu);
  }

  // Allocate needed arrays (boundary corrections, primitive moments, and nu*primitive moments)
  lbo->boundary_corrections = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);

  lbo->prim_moms = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
  lbo->nu_prim_moms = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
  lbo->m0 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  lbo->vtsq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  // Allocate moments needed for LBO update.
  vp_species_moment_init(app, s, &lbo->moms, "FiveMoments");

  // Edge of velocity space corrections to momentum and energy.
  lbo->bcorr_calc = gkyl_mom_calc_bcorr_lbo_vlasov_new(&s->grid, 
    &app->confBasis, &app->basis, v_bounds, app->use_gpu);
  
  // Primitive moment calculator.
  lbo->coll_pcalc = gkyl_prim_lbo_vlasov_calc_new(&s->grid, 
    &app->confBasis, &app->basis, &app->local, app->use_gpu);

  // LBO updater.
  struct gkyl_dg_lbo_vlasov_drag_auxfields drag_inp = { .nuSum = lbo->nu_sum, .nuPrimMomsSum = lbo->nu_prim_moms };
  struct gkyl_dg_lbo_vlasov_diff_auxfields diff_inp = { .nuSum = lbo->nu_sum, .nuPrimMomsSum = lbo->nu_prim_moms };
  lbo->coll_slvr = gkyl_dg_updater_lbo_vlasov_new(&s->grid, 
    &app->confBasis, &app->basis, &app->local, &drag_inp, &diff_inp, app->use_gpu);
}

void 
vp_species_lbo_cross_init(struct gkyl_vlasov_poisson_app *app, struct vp_species *s, struct vp_lbo_collisions *lbo)
{
  int vdim = app->vdim;

  lbo->cross_calc = gkyl_prim_lbo_vlasov_cross_calc_new(&s->grid, 
    &app->confBasis, &app->basis, &app->local, app->use_gpu);
  
  lbo->cross_nu_prim_moms = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);

  lbo->greene_factor_mem = 0;
  if (app->use_gpu)
    lbo->greene_factor_mem = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
  else
    lbo->greene_factor_mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);

  // Set pointers to species we cross-collide with.
  for (int i=0; i<lbo->num_cross_collisions; ++i) {
    lbo->collide_with[i] = vp_find_species(app, s->info.collisions.collide_with[i]);
    if (s->info.collisions.normNu) {
      double nuFrac = s->info.collisions.nuFrac ? s->info.collisions.nuFrac : 1.0;
      double eps0 = s->info.collisions.eps0 ? s->info.collisions.eps0: GKYL_EPSILON0;
      double hbar = s->info.collisions.hbar ? s->info.collisions.hbar: GKYL_PLANCKS_CONSTANT_H/2/M_PI;
      double eV = s->info.collisions.eV ? s->info.collisions.eV: GKYL_ELEMENTARY_CHARGE;
      lbo->cross_nu_fac[i] = nuFrac*gkyl_calc_norm_nu(s->info.collisions.n_ref, lbo->collide_with[i]->info.collisions.n_ref,
        s->info.mass, lbo->collide_with[i]->info.mass, s->info.charge, lbo->collide_with[i]->info.charge,
       	s->info.collisions.T_ref, lbo->collide_with[i]->info.collisions.T_ref, 0.0, eps0, hbar, eV);
    }

    lbo->other_m[i] = lbo->collide_with[i]->info.mass;
    lbo->other_prim_moms[i] = lbo->collide_with[i]->lbo.prim_moms;
    lbo->other_nu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->cross_prim_moms[i] = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
    lbo->cross_nu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->greene_num[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->greene_den[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->greene_factor[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    
    if (lbo->other_m[i] > s->info.mass) {
      gkyl_array_set(lbo->cross_nu[i], sqrt(2.0), lbo->self_nu);
      gkyl_array_set(lbo->other_nu[i], sqrt(2.0)*(s->info.mass)/(lbo->other_m[i]), lbo->self_nu);
    } else {
      gkyl_array_set(lbo->cross_nu[i], sqrt(2.0)*(lbo->other_m[i])/(s->info.mass), lbo->collide_with[i]->lbo.self_nu);
      gkyl_array_set(lbo->other_nu[i], sqrt(2.0), lbo->collide_with[i]->lbo.self_nu);
    }
    
    gkyl_array_accumulate(lbo->nu_sum, 1.0, lbo->cross_nu[i]);

    lbo->other_mnu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->other_mnu_m0[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    lbo->self_mnu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->self_mnu_m0[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    gkyl_array_set(lbo->self_mnu[i], s->info.mass, lbo->cross_nu[i]);
    gkyl_array_set(lbo->other_mnu[i], lbo->other_m[i], lbo->other_nu[i]);
  }
    
  lbo->betaGreenep1 = 1.0;
}

void
vp_species_lbo_moms(gkyl_vlasov_poisson_app *app, const struct vp_species *species,
  struct vp_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  // Computes moments, boundary corrections, and primitive moments.
  struct timespec wst = gkyl_wall_clock();

  // Compute needed moments.
  vp_species_moment_calc(&lbo->moms, species->local, app->local, fin);
  gkyl_array_set_range(lbo->m0, 1.0, lbo->moms.marr, &app->local);
  
  if (app->use_gpu) {
    // Construct boundary corrections.
    gkyl_mom_calc_bcorr_advance_cu(lbo->bcorr_calc,
      &species->local, &app->local, fin, lbo->boundary_corrections);

    // Construct primitive moments.
    gkyl_prim_lbo_calc_advance_cu(lbo->coll_pcalc, &app->local, 
      lbo->moms.marr, lbo->boundary_corrections,
      lbo->prim_moms);
  } 
  else {
    // Construct boundary corrections.
    gkyl_mom_calc_bcorr_advance(lbo->bcorr_calc,
      &species->local, &app->local, fin, lbo->boundary_corrections);

    // Construct primitive moments.
    gkyl_prim_lbo_calc_advance(lbo->coll_pcalc, &app->local, 
      lbo->moms.marr, lbo->boundary_corrections,
      lbo->prim_moms);
  }

  // Calculate self_nu if using spitzer-nu
  if (lbo->normNu) {
    gkyl_array_set_offset(lbo->vtsq, 1.0, lbo->prim_moms, app->vdim*app->confBasis.num_basis);
    gkyl_spitzer_coll_freq_advance_normnu(lbo->spitzer_calc, &app->local, lbo->vtsq, lbo->vtsq_min,
      lbo->m0, lbo->vtsq, lbo->vtsq_min, lbo->self_nu_fac, lbo->self_nu);
    gkyl_array_set(lbo->nu_sum, 1.0, lbo->self_nu);
  }

  // Scale upar and vtSq by nu.
  for (int d=0; d<app->vdim+1; d++)
    gkyl_dg_mul_op(app->confBasis, d, lbo->nu_prim_moms, d, lbo->prim_moms, 0, lbo->self_nu);
  
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

void
vp_species_lbo_cross_moms(gkyl_vlasov_poisson_app *app, const struct vp_species *species,
  struct vp_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  // Computes moments from cross-species collisions.
  struct timespec wst = gkyl_wall_clock();
  
  wst = gkyl_wall_clock();  
  for (int i=0; i<lbo->num_cross_collisions; ++i) {
    gkyl_dg_mul_op_range(app->confBasis, 0, lbo->self_mnu_m0[i], 0,
      lbo->self_mnu[i], 0, lbo->m0, &app->local);
    gkyl_dg_mul_op_range(app->confBasis, 0, lbo->other_mnu_m0[i], 0,
      lbo->other_mnu[i], 0, lbo->collide_with[i]->lbo.m0, &app->local);

    gkyl_dg_mul_op_range(app->confBasis, 0, lbo->greene_num[i], 0,
      lbo->other_mnu_m0[i], 0, lbo->m0, &app->local);

    gkyl_array_set(lbo->greene_den[i], 1.0, lbo->self_mnu_m0[i]);
    gkyl_array_accumulate(lbo->greene_den[i], 1.0, lbo->other_mnu_m0[i]);

    gkyl_dg_div_op_range(lbo->greene_factor_mem, app->confBasis, 0, lbo->greene_factor[i], 0,
      lbo->greene_num[i], 0, lbo->greene_den[i], &app->local);
    gkyl_array_scale(lbo->greene_factor[i], 2*lbo->betaGreenep1);

    if (app->use_gpu)
      gkyl_prim_lbo_cross_calc_advance_cu(lbo->cross_calc,
        &app->local, 
        lbo->greene_factor[i], 
        species->info.mass, lbo->moms.marr, lbo->prim_moms, 
        lbo->other_m[i], lbo->collide_with[i]->lbo.moms.marr, lbo->other_prim_moms[i],
        lbo->boundary_corrections, 
        lbo->cross_prim_moms[i]);
    else 
      gkyl_prim_lbo_cross_calc_advance(lbo->cross_calc,
        &app->local, 
        lbo->greene_factor[i], 
        species->info.mass, lbo->moms.marr, lbo->prim_moms, 
        lbo->other_m[i], lbo->collide_with[i]->lbo.moms.marr, lbo->other_prim_moms[i],
        lbo->boundary_corrections, 
        lbo->cross_prim_moms[i]);

    // Calculate cross nu if using spitzer nu.
    if (lbo->normNu) {
      gkyl_spitzer_coll_freq_advance_normnu(lbo->spitzer_calc, &app->local, lbo->vtsq, lbo->vtsq_min,
        lbo->collide_with[i]->lbo.m0, lbo->collide_with[i]->lbo.vtsq, lbo->collide_with[i]->lbo.vtsq_min,
        lbo->cross_nu_fac[i], lbo->cross_nu[i]);
      gkyl_array_accumulate(lbo->nu_sum, 1.0, lbo->cross_nu[i]);
    }

    // Scale upar_{rs} and vtSq_{rs} by nu_{rs}.
    for (int d=0; d<app->vdim+1; d++)
      gkyl_dg_mul_op(app->confBasis, d, lbo->cross_nu_prim_moms, d, lbo->cross_prim_moms[i], 0, lbo->cross_nu[i]);

    gkyl_array_accumulate(lbo->nu_prim_moms, 1.0, lbo->cross_nu_prim_moms);

  }
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

void
vp_species_lbo_rhs(gkyl_vlasov_poisson_app *app, const struct vp_species *species,
  struct vp_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  // Updates the collision terms in the rhs.
  struct timespec wst = gkyl_wall_clock();
    
  // Accumulate update due to collisions onto rhs.
  gkyl_dg_updater_lbo_vlasov_advance(lbo->coll_slvr, &species->local,
    fin, species->cflrate, rhs);
  
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void 
vp_species_lbo_release(const struct gkyl_vlasov_poisson_app *app, const struct vp_lbo_collisions *lbo)
{
  gkyl_array_release(lbo->boundary_corrections);
  gkyl_array_release(lbo->prim_moms);
  gkyl_array_release(lbo->self_nu);
  gkyl_array_release(lbo->nu_sum);
  gkyl_array_release(lbo->nu_prim_moms);
  gkyl_array_release(lbo->m0);
  gkyl_array_release(lbo->vtsq);

  vp_species_moment_release(app, &lbo->moms);

  gkyl_mom_calc_bcorr_release(lbo->bcorr_calc);
  gkyl_prim_lbo_calc_release(lbo->coll_pcalc);

  if (lbo->normNu) {
    gkyl_array_release(lbo->norm_nu);
    gkyl_array_release(lbo->nu_init);
    gkyl_spitzer_coll_freq_release(lbo->spitzer_calc);
  }

  if (lbo->num_cross_collisions) {
    gkyl_dg_bin_op_mem_release(lbo->greene_factor_mem);
    gkyl_array_release(lbo->cross_nu_prim_moms);
    for (int i=0; i<lbo->num_cross_collisions; ++i) {
      gkyl_array_release(lbo->cross_prim_moms[i]);
      gkyl_array_release(lbo->cross_nu[i]);
      gkyl_array_release(lbo->other_nu[i]);
      gkyl_array_release(lbo->self_mnu[i]);
      gkyl_array_release(lbo->self_mnu_m0[i]);
      gkyl_array_release(lbo->other_mnu[i]);
      gkyl_array_release(lbo->other_mnu_m0[i]);
      gkyl_array_release(lbo->greene_num[i]);
      gkyl_array_release(lbo->greene_den[i]);
      gkyl_array_release(lbo->greene_factor[i]);
    }
    gkyl_prim_lbo_cross_calc_release(lbo->cross_calc);
  }
  gkyl_dg_updater_lbo_vlasov_release(lbo->coll_slvr);
}
