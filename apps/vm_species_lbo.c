#include "gkyl_dg_bin_ops.h"
#include <assert.h>
#include <gkyl_vlasov_priv.h>

void
gkyl_lbo_write_mom(struct gkyl_vlasov_app* app, struct vm_lbo_collisions *lbo)
{      
  const char *fmt = "%s.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name);

  if (app->use_gpu) {
    gkyl_array_copy(lbo->boundary_corrections_host, lbo->boundary_corrections);
    gkyl_grid_sub_array_write(&app->grid, &app->local, lbo->boundary_corrections_host, fileNm);
  }
  else {
    gkyl_grid_sub_array_write(&app->grid, &app->local, lbo->boundary_corrections, fileNm);
  }
}

void 
vm_species_lbo_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_lbo_collisions *lbo, bool collides_with_fluid)
{
  int cdim = app->cdim, vdim = app->vdim;
  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = s->info.lower[d];
    v_bounds[d + vdim] = s->info.upper[d];
  }

  // allocate nu and initialize it
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

  lbo->boundary_corrections = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
  lbo->boundary_corrections_host = mkarr(false, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);

  lbo->u_drift = mkarr(app->use_gpu, vdim*app->confBasis.num_basis, app->local_ext.volume);
  lbo->vth_sq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  lbo->u_drift_host = mkarr(false, vdim*app->confBasis.num_basis, app->local_ext.volume);
  lbo->vth_sq_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  lbo->nu_u = mkarr(app->use_gpu, vdim*app->confBasis.num_basis, app->local_ext.volume);
  lbo->nu_vthsq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  lbo->m0 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    
  // allocate moments needed for LBO update
  vm_species_moment_init(app, s, &lbo->moms, "FiveMoments");

  // create collision equation object and solver
  if (app->use_gpu) {
    // edge of velocity space corrections to momentum and energy 
    lbo->bcorr_calc = gkyl_mom_calc_bcorr_lbo_vlasov_cu_dev_new(&s->grid, &app->confBasis, &app->basis, v_bounds);
    
    // primitive moment calculators
    if (collides_with_fluid)
      lbo->coll_pcalc = gkyl_prim_lbo_vlasov_with_fluid_calc_cu_dev_new(&s->grid, &app->confBasis, &app->basis, &app->local);
    else
      lbo->coll_pcalc = gkyl_prim_lbo_vlasov_calc_cu_dev_new(&s->grid, &app->confBasis, &app->basis);
  } else {
    // edge of velocity space corrections to momentum and energy 
    lbo->bcorr_calc = gkyl_mom_calc_bcorr_lbo_vlasov_new(&s->grid, &app->confBasis, &app->basis, v_bounds);
    
    // primitive moment calculators
    if (collides_with_fluid)
      lbo->coll_pcalc = gkyl_prim_lbo_vlasov_with_fluid_calc_new(&s->grid, &app->confBasis, &app->basis, &app->local);
    else
      lbo->coll_pcalc = gkyl_prim_lbo_vlasov_calc_new(&s->grid, &app->confBasis, &app->basis);
  }

  // LBO updater
  lbo->coll_slvr = gkyl_dg_updater_lbo_vlasov_new(&s->grid, &app->confBasis, &app->basis, &app->local, app->use_gpu);
}

void 
vm_species_lbo_cross_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_lbo_collisions *lbo)
{
  int vdim = app->vdim;

  if (app->use_gpu)
    lbo->cross_calc = gkyl_prim_lbo_vlasov_cross_calc_cu_dev_new(&s->grid, &app->confBasis, &app->basis);
  else
    lbo->cross_calc = gkyl_prim_lbo_vlasov_cross_calc_new(&s->grid, &app->confBasis, &app->basis);
  
  lbo->cross_nu_u = mkarr(app->use_gpu, vdim*app->confBasis.num_basis, app->local_ext.volume);
  lbo->cross_nu_vthsq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  lbo->greene_factor_mem = 0;
  if (app->use_gpu)
    lbo->greene_factor_mem = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
  else
    lbo->greene_factor_mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);

  // set pointers to species we cross-collide with
  for (int i=0; i<lbo->num_cross_collisions; ++i) {
    lbo->collide_with[i] = vm_find_species(app, s->info.collisions.collide_with[i]);
    lbo->other_m[i] = lbo->collide_with[i]->info.mass;
    lbo->other_u_drift[i] = lbo->collide_with[i]->lbo.u_drift;
    lbo->other_vth_sq[i] = lbo->collide_with[i]->lbo.vth_sq;
    lbo->other_nu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->cross_u_drift[i] = mkarr(app->use_gpu, vdim*app->confBasis.num_basis, app->local_ext.volume);
    lbo->cross_vth_sq[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->cross_nu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->greene_num[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->greene_den[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->greene_factor[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    
    if (lbo->other_m[i] > s->info.mass) {
      gkyl_array_set(lbo->cross_nu[i], sqrt(2), lbo->self_nu);
      gkyl_array_set(lbo->other_nu[i], (s->info.mass)/(lbo->other_m[i]), lbo->self_nu);
    } else {
      gkyl_array_set(lbo->cross_nu[i], (lbo->other_m[i])/(s->info.mass), lbo->collide_with[i]->lbo.self_nu);
      gkyl_array_set(lbo->other_nu[i], sqrt(2), lbo->collide_with[i]->lbo.self_nu);
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

// computes moments, boundary corrections, and primitive moments
void
vm_species_lbo_moms(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_lbo_collisions *lbo, const struct gkyl_array *fin,
  bool collides_with_fluid, const struct gkyl_array *fluidin[])
{
  struct timespec wst = gkyl_wall_clock();
  // compute needed moments
  vm_species_moment_calc(&lbo->moms, species->local, app->local, fin);
  gkyl_array_set_range(lbo->m0, 1.0, lbo->moms.marr, app->local);

  // get pointer to fluid species if kinetic species is colliding with fluid species
  if (collides_with_fluid)
    gkyl_prim_lbo_vlasov_with_fluid_set_auxfields(gkyl_prim_lbo_calc_get_prim(lbo->coll_pcalc),
      (struct gkyl_prim_lbo_vlasov_with_fluid_auxfields) { .fluid = fluidin[species->fluid_index]  });
  
  if (app->use_gpu) {
    wst = gkyl_wall_clock();

    // construct boundary corrections
    gkyl_mom_calc_bcorr_advance_cu(lbo->bcorr_calc,
      &species->local, &app->local, fin, lbo->boundary_corrections);

    gkyl_lbo_write_mom(app, lbo);

    // construct primitive moments    
    gkyl_prim_lbo_calc_advance_cu(lbo->coll_pcalc, app->confBasis, &app->local, 
      lbo->moms.marr, lbo->boundary_corrections,
      lbo->u_drift, lbo->vth_sq);
  
    gkyl_dg_mul_op(app->confBasis, 0, lbo->nu_u, 0, lbo->u_drift, 0, lbo->self_nu);
    gkyl_dg_mul_op(app->confBasis, 0, lbo->nu_vthsq, 0, lbo->vth_sq, 0, lbo->self_nu);

    app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
  } else {
    wst = gkyl_wall_clock();
    
    // construct boundary corrections
    gkyl_mom_calc_bcorr_advance(lbo->bcorr_calc,
      &species->local, &app->local, fin, lbo->boundary_corrections);

    gkyl_lbo_write_mom(app, lbo);

    // construct primitive moments    
    gkyl_prim_lbo_calc_advance(lbo->coll_pcalc, app->confBasis, &app->local, 
      lbo->moms.marr, lbo->boundary_corrections,
      lbo->u_drift, lbo->vth_sq);
    
    gkyl_dg_mul_op(app->confBasis, 0, lbo->nu_u, 0, lbo->u_drift, 0, lbo->self_nu);
    gkyl_dg_mul_op(app->confBasis, 0, lbo->nu_vthsq, 0, lbo->vth_sq, 0, lbo->self_nu);

    app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
  }
}

// computes moments from cross-species collisions
void
vm_species_lbo_cross_moms(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_lbo_collisions *lbo, const struct gkyl_array *fin,
  bool collides_with_fluid, const struct gkyl_array *fluidin[])
{
  struct timespec wst = gkyl_wall_clock();

  // get pointer to fluid species if kinetic species is colliding with fluid species
  if (collides_with_fluid)
    gkyl_prim_lbo_vlasov_with_fluid_set_auxfields(gkyl_prim_lbo_cross_calc_get_prim(lbo->cross_calc),
      (struct gkyl_prim_lbo_vlasov_with_fluid_auxfields) { .fluid = fluidin[species->fluid_index]  });
  
  wst = gkyl_wall_clock();  
  for (int i=0; i<lbo->num_cross_collisions; ++i) {
    gkyl_dg_mul_op_range(app->confBasis, 0, lbo->self_mnu_m0[i], 0,
      lbo->self_mnu[i], 0, lbo->m0, app->local);
    gkyl_dg_mul_op_range(app->confBasis, 0, lbo->other_mnu_m0[i], 0,
      lbo->other_mnu[i], 0, lbo->collide_with[i]->lbo.m0, app->local);

    gkyl_dg_mul_op_range(app->confBasis, 0, lbo->greene_num[i], 0,
      lbo->other_mnu_m0[i], 0, lbo->m0, app->local);

    gkyl_array_set(lbo->greene_den[i], 1.0, lbo->self_mnu_m0[i]);
    gkyl_array_accumulate(lbo->greene_den[i], 1.0, lbo->other_mnu_m0[i]);

    gkyl_dg_div_op_range(lbo->greene_factor_mem, app->confBasis, 0, lbo->greene_factor[i], 0,
      lbo->greene_num[i], 0, lbo->greene_den[i], app->local);
    gkyl_array_scale(lbo->greene_factor[i], 2*lbo->betaGreenep1);

    if (app->use_gpu)
      gkyl_prim_lbo_cross_calc_advance_cu(lbo->cross_calc,
        app->confBasis, &app->local, 
        lbo->greene_factor[i], 
        species->info.mass, lbo->u_drift, lbo->vth_sq, 
        lbo->other_m[i], lbo->other_u_drift[i], lbo->other_vth_sq[i],
        lbo->moms.marr, lbo->boundary_corrections, 
        lbo->cross_u_drift[i], lbo->cross_vth_sq[i]);
    else 
      gkyl_prim_lbo_cross_calc_advance(lbo->cross_calc,
        app->confBasis, &app->local, 
        lbo->greene_factor[i], 
        species->info.mass, lbo->u_drift, lbo->vth_sq, 
        lbo->other_m[i], lbo->other_u_drift[i], lbo->other_vth_sq[i],
        lbo->moms.marr, lbo->boundary_corrections, 
        lbo->cross_u_drift[i], lbo->cross_vth_sq[i]);

    gkyl_dg_mul_op(app->confBasis, 0, lbo->cross_nu_u, 0, lbo->cross_u_drift[i], 0, lbo->cross_nu[i]);
    gkyl_dg_mul_op(app->confBasis, 0, lbo->cross_nu_vthsq, 0, lbo->cross_vth_sq[i], 0, lbo->cross_nu[i]);
    gkyl_array_accumulate(lbo->nu_u, 1.0, lbo->cross_nu_u);
    gkyl_array_accumulate(lbo->nu_vthsq, 1.0, lbo->cross_nu_vthsq);
  }
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

// computes cross-primitive moments and updates the collision terms in the rhs
double
vm_species_lbo_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  if (app->use_gpu) {
    wst = gkyl_wall_clock();
    
    // accumulate update due to collisions onto rhs
    gkyl_dg_updater_lbo_vlasov_advance_cu(lbo->coll_slvr, &species->local,
      lbo->nu_sum, lbo->nu_u, lbo->nu_vthsq, fin, species->cflrate, rhs);
    
    app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
    
  } else {
    wst = gkyl_wall_clock();
    
    // accumulate update due to collisions onto rhs
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
  gkyl_array_release(lbo->boundary_corrections_host);
  gkyl_array_release(lbo->u_drift_host);
  gkyl_array_release(lbo->vth_sq_host);
  gkyl_array_release(lbo->self_nu);
  gkyl_array_release(lbo->nu_sum);
  gkyl_array_release(lbo->nu_u);
  gkyl_array_release(lbo->nu_vthsq);
  gkyl_array_release(lbo->m0);

  vm_species_moment_release(app, &lbo->moms);

  gkyl_mom_calc_bcorr_release(lbo->bcorr_calc);
  gkyl_prim_lbo_calc_release(lbo->coll_pcalc);

  if (lbo->num_cross_collisions) {
    gkyl_dg_bin_op_mem_release(lbo->greene_factor_mem);
    gkyl_array_release(lbo->cross_nu_u);
    gkyl_array_release(lbo->cross_nu_vthsq);
    for (int i=0; i<lbo->num_cross_collisions; ++i) {
      gkyl_array_release(lbo->cross_u_drift[i]);
      gkyl_array_release(lbo->cross_vth_sq[i]);
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
