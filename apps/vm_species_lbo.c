#include "gkyl_dg_bin_ops.h"
#include "gkyl_eqn_type.h"
#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_lbo_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_lbo_collisions *lbo)
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

  lbo->model_id = GKYL_MODEL_DEFAULT;
  if (s->model_id == GKYL_MODEL_PKPM) {
    lbo->model_id = GKYL_MODEL_PKPM;
    // Only energy corrections for pkpm model
    lbo->boundary_corrections = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);   
    lbo->prim_moms = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->nu_prim_moms = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    // edge of velocity space corrections to *only* energy (in PKPM model)
    lbo->bcorr_calc = gkyl_mom_calc_bcorr_lbo_vlasov_pkpm_new(&s->grid, 
      &app->confBasis, &app->basis, v_bounds, s->info.mass, app->use_gpu);
    // primitive moment calculators
    lbo->coll_pcalc = gkyl_prim_lbo_vlasov_pkpm_calc_new(&s->grid, 
      &app->confBasis, &app->basis, &app->local, app->use_gpu);

    // Since PKPM used vth^2 to penalize fluxes, need to apply BCs to vth^2
    // allocate buffer for applying BCs
    long buff_sz = 0;
    // compute buffer size needed
    for (int d=0; d<app->cdim; ++d) {
      long vol = app->skin_ghost.lower_skin[d].volume;
      buff_sz = buff_sz > vol ? buff_sz : vol;
    }
    lbo->bc_buffer = mkarr(app->use_gpu, app->confBasis.num_basis, buff_sz);

    // determine which directions are not periodic
    int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
    for (int d=0; d<num_periodic_dir; ++d)
      is_np[app->periodic_dirs[d]] = 0;

    for (int dir=0; dir<app->cdim; ++dir) {
      lbo->lower_bc[dir] = lbo->upper_bc[dir] = GKYL_SPECIES_COPY;
      if (is_np[dir]) {
        const enum gkyl_species_bc_type *bc;
        // Use vm_species to get BC type if not periodic
        if (dir == 0)
          bc = s->info.bcx;
        else if (dir == 1)
          bc = s->info.bcy;
        else
          bc = s->info.bcz;

        lbo->lower_bc[dir] = bc[0];
        lbo->upper_bc[dir] = bc[1];
      }
    }

    int ghost[GKYL_MAX_DIM] = {0.0};
    for (int d=0; d<app->cdim; ++d)
      ghost[d] = 1;

    for (int d=0; d<app->cdim; ++d) {
      // Lower BC updater. Copy BCs by default.
      enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
      if (lbo->lower_bc[d] == GKYL_SPECIES_COPY)
        bctype = GKYL_BC_COPY;
      else if (lbo->lower_bc[d] == GKYL_SPECIES_ABSORB)
        bctype = GKYL_BC_ABSORB;

      lbo->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, &app->local_ext, ghost, bctype,
                                      app->basis_on_dev.confBasis, lbo->prim_moms->ncomp, app->cdim, app->use_gpu);
      // Upper BC updater. Copy BCs by default.
      if (lbo->upper_bc[d] == GKYL_SPECIES_COPY)
        bctype = GKYL_BC_COPY;
      else if (lbo->upper_bc[d] == GKYL_SPECIES_ABSORB)
        bctype = GKYL_BC_ABSORB;

      lbo->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, &app->local_ext, ghost, bctype,
                                      app->basis_on_dev.confBasis, lbo->prim_moms->ncomp, app->cdim, app->use_gpu);
    }
  }
  else {
    lbo->boundary_corrections = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);

    lbo->prim_moms = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
    lbo->nu_prim_moms = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
    lbo->m0 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
      
    // allocate moments needed for LBO update
    vm_species_moment_init(app, s, &lbo->moms, "FiveMoments");

    // edge of velocity space corrections to momentum and energy 
    lbo->bcorr_calc = gkyl_mom_calc_bcorr_lbo_vlasov_new(&s->grid, 
      &app->confBasis, &app->basis, v_bounds, app->use_gpu);
    
    // primitive moment calculator
    lbo->coll_pcalc = gkyl_prim_lbo_vlasov_calc_new(&s->grid, 
      &app->confBasis, &app->basis, &app->local, app->use_gpu);
  }

  // LBO updater
  lbo->coll_slvr = gkyl_dg_updater_lbo_vlasov_new(&s->grid, 
    &app->confBasis, &app->basis, &app->local, s->model_id, app->use_gpu);
}

void 
vm_species_lbo_cross_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_lbo_collisions *lbo)
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

  // set pointers to species we cross-collide with
  for (int i=0; i<lbo->num_cross_collisions; ++i) {
    lbo->collide_with[i] = vm_find_species(app, s->info.collisions.collide_with[i]);
    lbo->other_m[i] = lbo->collide_with[i]->info.mass;
    lbo->other_prim_moms[i] = lbo->collide_with[i]->lbo.prim_moms;
    lbo->other_nu[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    lbo->cross_prim_moms[i] = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
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
  struct vm_lbo_collisions *lbo, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();
  
  if (species->model_id == GKYL_MODEL_PKPM) {
    // Set pointer to pressure tensor for use in boundary corrections  
    gkyl_prim_lbo_vlasov_pkpm_set_auxfields(gkyl_prim_lbo_calc_get_prim(lbo->coll_pcalc),
      (struct gkyl_prim_lbo_vlasov_pkpm_auxfields) { .pvar = species->pkpm_fluid_species->p });
  }
  else {
    // compute needed moments
    vm_species_moment_calc(&lbo->moms, species->local, app->local, fin);
    gkyl_array_set_range(lbo->m0, 1.0, lbo->moms.marr, app->local);
  }
  
  if (app->use_gpu) {
    wst = gkyl_wall_clock();

    // construct boundary corrections
    gkyl_mom_calc_bcorr_advance_cu(lbo->bcorr_calc,
      &species->local, &app->local, fin, lbo->boundary_corrections);

    // construct primitive moments  
    if (species->model_id == GKYL_MODEL_PKPM) {
      // PKPM moments already computed before this, so just fetch results
      gkyl_prim_lbo_calc_advance_cu(lbo->coll_pcalc, &app->local, 
        species->pkpm_moms.marr, lbo->boundary_corrections,
        lbo->prim_moms);
      // Apply BCs to vth^2 for use in penalization
      vm_species_lbo_apply_bc(app, lbo);
      // PKPM model only has vth^2
      gkyl_dg_mul_op(app->confBasis, 0, lbo->nu_prim_moms, 0, lbo->prim_moms, 0, lbo->self_nu);
    }
    else {
      gkyl_prim_lbo_calc_advance_cu(lbo->coll_pcalc, &app->local, 
        lbo->moms.marr, lbo->boundary_corrections,
        lbo->prim_moms);

      for (int d=0; d<app->vdim; d++)
        gkyl_dg_mul_op(app->confBasis, d, lbo->nu_prim_moms, d, lbo->prim_moms, 0, lbo->self_nu);
      gkyl_dg_mul_op(app->confBasis, app->vdim, lbo->nu_prim_moms, app->vdim, lbo->prim_moms, 0, lbo->self_nu);
    }  

    app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
  } 
  else {
    wst = gkyl_wall_clock();
    
    // construct boundary corrections
    gkyl_mom_calc_bcorr_advance(lbo->bcorr_calc,
      &species->local, &app->local, fin, lbo->boundary_corrections);

    // construct primitive moments  
    if (species->model_id == GKYL_MODEL_PKPM) {
      // PKPM moments already computed before this, so just fetch results
      gkyl_prim_lbo_calc_advance(lbo->coll_pcalc, &app->local, 
        species->pkpm_moms.marr, lbo->boundary_corrections,
        lbo->prim_moms);
      // Apply BCs to vth^2 for use in penalization
      vm_species_lbo_apply_bc(app, lbo);
      // PKPM model only has vth^2
      gkyl_dg_mul_op(app->confBasis, 0, lbo->nu_prim_moms, 0, lbo->prim_moms, 0, lbo->self_nu);
    }
    else {
      gkyl_prim_lbo_calc_advance(lbo->coll_pcalc, &app->local, 
        lbo->moms.marr, lbo->boundary_corrections,
        lbo->prim_moms);

      for (int d=0; d<app->vdim; d++)
        gkyl_dg_mul_op(app->confBasis, d, lbo->nu_prim_moms, d, lbo->prim_moms, 0, lbo->self_nu);
      gkyl_dg_mul_op(app->confBasis, app->vdim, lbo->nu_prim_moms, app->vdim, lbo->prim_moms, 0, lbo->self_nu);
    }  

    app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
  }
}

// computes moments from cross-species collisions
void
vm_species_lbo_cross_moms(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_lbo_collisions *lbo, const struct gkyl_array *fin)
{
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


    for (int d=0; d<app->vdim; d++)
      gkyl_dg_mul_op(app->confBasis, d, lbo->cross_nu_prim_moms, d, lbo->cross_prim_moms[i], 0, lbo->cross_nu[i]);
    gkyl_dg_mul_op(app->confBasis, app->vdim, lbo->cross_nu_prim_moms, app->vdim, lbo->cross_prim_moms[i], 0, lbo->cross_nu[i]);

    gkyl_array_accumulate(lbo->nu_prim_moms, 1.0, lbo->cross_nu_prim_moms);
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
      lbo->nu_sum, lbo->nu_prim_moms, fin, species->cflrate, rhs);
    
    app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
    
  } else {
    wst = gkyl_wall_clock();
    
    // accumulate update due to collisions onto rhs
    gkyl_dg_updater_lbo_vlasov_advance(lbo->coll_slvr, &species->local,
      lbo->nu_sum, lbo->nu_prim_moms, fin, species->cflrate, rhs);
    
    app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
  }
  // TODO: This needs to be set properly!
  return 0;
}

// Determine which directions are periodic and which directions are copy,
// and then apply boundary conditions to primitive moments (only used by PKPM)
void
vm_species_lbo_apply_bc(gkyl_vlasov_app *app, const struct vm_lbo_collisions *lbo)
{
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d) {
    gkyl_array_copy_to_buffer(lbo->bc_buffer->data, lbo->prim_moms, app->skin_ghost.lower_skin[d]);
    gkyl_array_copy_from_buffer(lbo->prim_moms, lbo->bc_buffer->data, app->skin_ghost.upper_ghost[d]);

    gkyl_array_copy_to_buffer(lbo->bc_buffer->data, lbo->prim_moms, app->skin_ghost.upper_skin[d]);
    gkyl_array_copy_from_buffer(lbo->prim_moms, lbo->bc_buffer->data, app->skin_ghost.lower_ghost[d]);
    is_np_bc[app->periodic_dirs[d]] = 0;
  }
  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {

      switch (lbo->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(lbo->bc_lo[d], lbo->bc_buffer, lbo->prim_moms);
          break;
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
        case GKYL_SPECIES_FIXED_FUNC:
          assert(false);
          break;
        default:
          break;
      }

      switch (lbo->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(lbo->bc_up[d], lbo->bc_buffer, lbo->prim_moms);
          break;
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
        case GKYL_SPECIES_FIXED_FUNC:
          assert(false);
          break;
        default:
          break;
      }
    }
  }
}

void 
vm_species_lbo_release(const struct gkyl_vlasov_app *app, const struct vm_lbo_collisions *lbo)
{
  gkyl_array_release(lbo->boundary_corrections);
  gkyl_array_release(lbo->prim_moms);
  gkyl_array_release(lbo->self_nu);
  gkyl_array_release(lbo->nu_sum);
  gkyl_array_release(lbo->nu_prim_moms);

  if (lbo->model_id != GKYL_MODEL_PKPM) {
    gkyl_array_release(lbo->m0);
    vm_species_moment_release(app, &lbo->moms);
  }
  else {
    // Copy BCs are allocated by default for PKPM. Need to free.
    for (int d=0; d<app->cdim; ++d) {
      gkyl_bc_basic_release(lbo->bc_lo[d]);
      gkyl_bc_basic_release(lbo->bc_up[d]);
    }
    gkyl_array_release(lbo->bc_buffer);
  }

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
  gkyl_dg_updater_lbo_vlasov_release(lbo->coll_slvr);
 }
