#include <gkyl_flash_priv.h>

void
flash_species_lbo_init(gkyl_flash_app *app, struct flash_species *s, struct flash_lbo_collisions *lbo) 
{
  int cdim = app->cdim;
  int vdim = app->vdim;
  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = s->info.lower[d];
    v_bounds[d+vdim] = s->info.upper[d];
  }

  // allocate nu and initialize it
  lbo->nu_sum = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_kinetic_ext.volume);
  lbo->self_nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_kinetic_ext.volume);
  struct gkyl_array *self_nu = mkarr(false, app->confBasis.num_basis, app->local_kinetic_ext.volume);

  lbo->num_cross_collisions = s->info.collisions.num_cross_collisions;
  
  // Verify that user specified collision frequency profile
  if (!s->info.collisions.self_nu) {
    gkyl_exit("No collision frequency profile. Specify `self_nu` collision frequency function in the species collision object.\n");
  }
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid_kinetic, &app->confBasis,
    app->poly_order+1, 1, s->info.collisions.self_nu, s->info.collisions.ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local_kinetic, self_nu);
  gkyl_proj_on_basis_release(proj);

  gkyl_array_copy(lbo->self_nu, self_nu);
  gkyl_array_copy(lbo->nu_sum, self_nu);
  gkyl_array_release(self_nu);

  // Allocate needed arrays (boundary corrections, primitive moments, and nu*primitive moments)
  lbo->boundary_corrections = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_kinetic_ext.volume);

  lbo->prim_moms = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_kinetic_ext.volume);
  lbo->nu_prim_moms = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_kinetic_ext.volume);
  lbo->m0 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_kinetic_ext.volume);

  // edge of velocity space corrections to momentum and energy 
  // Note: PKPM model still has a momentum correction to insure M1 = 0 from collisions
  lbo->bcorr_calc = gkyl_mom_calc_bcorr_lbo_pkpm_new(&s->grid, 
    &app->confBasis, &app->basis, v_bounds, s->info.mass, app->use_gpu);
  // primitive moment calculators
  lbo->coll_pcalc = gkyl_prim_lbo_pkpm_calc_new(&s->grid, 
    &app->confBasis, &app->basis, &app->local_kinetic, app->use_gpu);

      
  // LBO updater
  struct gkyl_dg_lbo_pkpm_drag_auxfields drag_inp = { .nuSum = lbo->nu_sum, .nuPrimMomsSum = lbo->nu_prim_moms };
  struct gkyl_dg_lbo_pkpm_diff_auxfields diff_inp = { .nuSum = lbo->nu_sum, .nuPrimMomsSum = lbo->nu_prim_moms };
  lbo->coll_slvr = gkyl_dg_updater_lbo_pkpm_new(&s->grid, 
    &app->confBasis, &app->basis, &app->local_kinetic, &drag_inp, &diff_inp, app->use_gpu);
}

// Compute required moments, primitive moments, and boundary corrections.
void
flash_species_lbo_moms(gkyl_flash_app *app, const struct flash_species *species,
  struct flash_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_set_range(lbo->m0, 1.0, species->pkpm_moms, &app->local_kinetic);

  // construct boundary corrections
  gkyl_mom_calc_bcorr_advance(lbo->bcorr_calc,
    &species->local, &app->local_kinetic, fin, lbo->boundary_corrections);

  // Compute primitive moments.
  gkyl_prim_lbo_calc_advance(lbo->coll_pcalc, &app->local_kinetic, 
    species->pkpm_moms, lbo->boundary_corrections, lbo->self_nu, lbo->prim_moms);

  for (int d=0; d<app->vdim; d++)
    gkyl_dg_mul_op(&app->confBasis, d, lbo->nu_prim_moms, d, lbo->prim_moms, 0, lbo->self_nu);
  gkyl_dg_mul_op(&app->confBasis, app->vdim, lbo->nu_prim_moms, app->vdim, lbo->prim_moms, 0, lbo->self_nu);
  
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);   
}

void
flash_species_lbo_rhs(gkyl_flash_app *app, const struct flash_species *species, 
  struct flash_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
    
  // accumulate update due to collisions onto rhs
  gkyl_dg_updater_lbo_pkpm_advance(lbo->coll_slvr, &species->local,
    fin, species->cflrate_f, rhs);
  
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void
flash_species_lbo_compute_cflrate(gkyl_flash_app *app, const struct flash_species *species,
  struct flash_lbo_collisions *lbo)
{
  gkyl_dg_updater_flash_compute_cflrate_lbo(species->updater, &species->grid, 
    &species->local, &app->local_kinetic, lbo->nu_sum, lbo->nu_prim_moms, species->cflrate_f);
}

void 
flash_species_lbo_release(const struct flash_lbo_collisions *lbo)
{
  gkyl_array_release(lbo->boundary_corrections);
  gkyl_array_release(lbo->prim_moms);
  gkyl_array_release(lbo->self_nu);
  gkyl_array_release(lbo->nu_sum);
  gkyl_array_release(lbo->nu_prim_moms);
  gkyl_array_release(lbo->m0);

  gkyl_mom_calc_bcorr_release(lbo->bcorr_calc);
  gkyl_prim_lbo_calc_release(lbo->coll_pcalc);

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
  gkyl_dg_updater_lbo_pkpm_release(lbo->coll_slvr);
 }