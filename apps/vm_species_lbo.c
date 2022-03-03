#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_lbo_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_lbo_collisions *lbo)
{
  // TO DO: Expose nu_u and nu_vthsq arrays above species object
  //        for cross-species collisions. Just testing for now JJ 09/24/21
  int cdim = app->cdim, vdim = app->vdim;
  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = s->info.lower[d];
    v_bounds[d + vdim] = s->info.upper[d];
  }

  // allocate nu and initialize it
  lbo->nu_sum = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array *nu_sum = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    app->poly_order+1, 1, s->info.nu, s->info.ctx);

  gkyl_proj_on_basis_advance(proj, 0.0, &app->local, nu_sum);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(lbo->nu_sum, nu_sum);
  gkyl_array_release(nu_sum);

  lbo->boundary_corrections = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume); 

  lbo->u_drift = mkarr(app->use_gpu, vdim*app->confBasis.num_basis, app->local_ext.volume);
  lbo->vth_sq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  lbo->nu_u = mkarr(app->use_gpu, vdim*app->confBasis.num_basis, app->local_ext.volume);
  lbo->nu_vthsq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  // allocate moments needed for LBO update
  vm_species_moment_init(app, s, &lbo->moms, "FiveMoments");
  // create collision equation object and solver
  if (app->use_gpu) {
    // edge of velocity space corrections to momentum and energy 
    lbo->bcorr_type = gkyl_mom_bcorr_lbo_vlasov_cu_dev_new(&app->confBasis, &app->basis, v_bounds);
    lbo->bcorr_calc = gkyl_mom_calc_bcorr_cu_dev_new(&s->grid, lbo->bcorr_type);
    
    // primitive moment calculators
    lbo->coll_prim = gkyl_prim_lbo_vlasov_cu_dev_new(&app->confBasis, &app->basis);
    lbo->coll_pcalc = gkyl_prim_lbo_calc_cu_dev_new(&s->grid, lbo->coll_prim);
    
    // LBO updater
    lbo->coll_slvr = gkyl_dg_updater_lbo_vlasov_cu_dev_new(&s->grid, &app->confBasis, &app->basis, &app->local);
   }
  else {
    // edge of velocity space corrections to momentum and energy 
    lbo->bcorr_type = gkyl_mom_bcorr_lbo_vlasov_new(&app->confBasis, &app->basis, v_bounds);
    lbo->bcorr_calc = gkyl_mom_calc_bcorr_new(&s->grid, lbo->bcorr_type);
    
    // primitive moment calculators
    lbo->coll_prim = gkyl_prim_lbo_vlasov_new(&app->confBasis, &app->basis);
    lbo->coll_pcalc = gkyl_prim_lbo_calc_new(&s->grid, lbo->coll_prim);
    
    // LBO updater
    lbo->coll_slvr = gkyl_dg_updater_lbo_vlasov_new(&s->grid, &app->confBasis, &app->basis, &app->local);
  }
}

double
vm_species_lbo_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  // compute needed moments
  vm_species_moment_calc(&lbo->moms, species->local, app->local, fin);

  if (app->use_gpu) {
    wst = gkyl_wall_clock();
    
    // construct boundary corrections
    gkyl_mom_calc_bcorr_advance_cu(lbo->bcorr_calc,
      &species->local, &app->local, fin, lbo->boundary_corrections);

    // construct primitive moments
    gkyl_prim_lbo_calc_advance_cu(lbo->coll_pcalc, app->confBasis, app->local, 
      lbo->moms.marr, lbo->boundary_corrections,
      lbo->u_drift, lbo->vth_sq);
  
    gkyl_dg_mul_op(app->confBasis, 0, lbo->nu_u, 0, lbo->u_drift, 0, lbo->nu_sum);
    gkyl_dg_mul_op(app->confBasis, 0, lbo->nu_vthsq, 0, lbo->vth_sq, 0, lbo->nu_sum);

    app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);

    wst = gkyl_wall_clock();

    // acccumulate update due to collisions onto rhs
    gkyl_dg_updater_lbo_vlasov_advance_cu(lbo->coll_slvr, &species->local,
      lbo->nu_sum, lbo->nu_u, lbo->nu_vthsq, fin, species->cflrate, rhs);
    
    app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
    
  } else {
    wst = gkyl_wall_clock();    
    // construct boundary corrections
    gkyl_mom_calc_bcorr_advance(lbo->bcorr_calc,
      &species->local, &app->local, fin, lbo->boundary_corrections);

    // construct primitive moments
    gkyl_prim_lbo_calc_advance(lbo->coll_pcalc, app->confBasis, app->local, 
      lbo->moms.marr, lbo->boundary_corrections,
      lbo->u_drift, lbo->vth_sq);
  
    gkyl_dg_mul_op(app->confBasis, 0, lbo->nu_u, 0, lbo->u_drift, 0, lbo->nu_sum);
    gkyl_dg_mul_op(app->confBasis, 0, lbo->nu_vthsq, 0, lbo->vth_sq, 0, lbo->nu_sum);

    app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    

    wst = gkyl_wall_clock();
    
    // acccumulate update due to collisions onto rhs
    gkyl_dg_updater_lbo_vlasov_advance(lbo->coll_slvr, &species->local,
      lbo->nu_sum, lbo->nu_u, lbo->nu_vthsq, fin, species->cflrate, rhs);
    
    app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
  }
  // TODO: This needs to be set properly!
  return 0;
}

void 
vm_species_lbo_release(const struct gkyl_vlasov_app *app, const struct vm_lbo_collisions *lbo)
{
  gkyl_array_release(lbo->boundary_corrections);
  gkyl_array_release(lbo->u_drift);
  gkyl_array_release(lbo->vth_sq);
  gkyl_array_release(lbo->nu_sum);
  gkyl_array_release(lbo->nu_u);
  gkyl_array_release(lbo->nu_vthsq);

  vm_species_moment_release(app, &lbo->moms);

  gkyl_mom_type_release(lbo->bcorr_type);
  gkyl_mom_calc_bcorr_release(lbo->bcorr_calc);
  gkyl_prim_lbo_type_release(lbo->coll_prim);
  gkyl_prim_lbo_calc_release(lbo->coll_pcalc);
  gkyl_dg_updater_lbo_vlasov_release(lbo->coll_slvr);
 }
