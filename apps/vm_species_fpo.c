#include <assert.h>
#include <gkyl_vlasov_priv.h>

#include <gkyl_array_rio.h>

void 
vm_species_fpo_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_fpo_collisions *fpo)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;
  struct gkyl_basis surf_basis;

  // initialize surface basis for potentials on velocity space edges
  gkyl_cart_modal_serendip(&surf_basis, pdim-1, app->poly_order);

  // allocate gamma and initialize it
  fpo->gamma = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array *gamma_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    app->poly_order+1, 1, s->info.collisions.self_nu, s->info.collisions.ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local, gamma_host);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(fpo->gamma, gamma_host);
  gkyl_array_release(gamma_host);

  // initialize the potentials and solver for potentials
  fpo->h = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  fpo->g = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
 
  fpo->h_surf = mkarr(app->use_gpu, 3*surf_basis.num_basis, s->local_ext.volume);
  fpo->g_surf = mkarr(app->use_gpu, 3*surf_basis.num_basis, s->local_ext.volume); 
  fpo->dhdv_surf = mkarr(app->use_gpu, 3*surf_basis.num_basis, s->local_ext.volume);
  fpo->dgdv_surf = mkarr(app->use_gpu, 9*surf_basis.num_basis, s->local_ext.volume); 
  fpo->d2gdv2_surf = mkarr(app->use_gpu, 9*surf_basis.num_basis, s->local_ext.volume); 

  fpo->m0 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  // fpo->boundary_corrections = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
  fpo->prim_moms = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);

  fpo->pot_slvr = gkyl_proj_maxwellian_pots_on_basis_new(&s->grid, &app->confBasis, &app->basis, app->poly_order+1);

  // allocate moments needed for FPO update
  vm_species_moment_init(app, s, &fpo->moms, "FiveMoments");
  
  fpo->prim_calc = gkyl_dg_prim_vars_vlasov_new(&app->confBasis, &app->basis, "prim", app->use_gpu);

  // initialize drag and diffusion coefficients
  fpo->drag_coeff = mkarr(app->use_gpu, 3*app->basis.num_basis, s->local_ext.volume);
  fpo->diff_coeff = mkarr(app->use_gpu, 9*app->basis.num_basis, s->local_ext.volume);

  // initialize FPO updater
  fpo->coll_slvr = gkyl_dg_updater_fpo_vlasov_new(&s->grid, &app->basis, &s->local, app->use_gpu);
}

// computes drag coefficient and diffusion tensor
void
vm_species_fpo_drag_diff_coeffs(gkyl_vlasov_app *app, const struct vm_species *s,
  struct vm_fpo_collisions *fpo, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  // Calculate five moments
  vm_species_moment_calc(&fpo->moms, s->local, app->local, fin);
  gkyl_array_set_range(fpo->m0, 1.0, fpo->moms.marr, &app->local);

  // Calculate primitive moments vtsq and u_i
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &app->local);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linc = gkyl_range_idx(&app->local, conf_iter.idx);

    const double *moms_d = gkyl_array_cfetch(fpo->moms.marr, linc);
    double *prim_moms_d = gkyl_array_fetch(fpo->prim_moms, linc);

    fpo->prim_calc->kernel(fpo->prim_calc, conf_iter.idx, moms_d, prim_moms_d);
  }

  // calculate maxwellian potentials
  gkyl_proj_maxwellian_pots_on_basis_lab_mom(fpo->pot_slvr, &s->local, &app->local,
    fpo->m0, fpo->prim_moms, 
    fpo->h, fpo->g, fpo->h_surf, fpo->g_surf, fpo->dhdv_surf, fpo->dgdv_surf, fpo->d2gdv2_surf);

  // Calculate drag and diffusion coefficients
  gkyl_calc_fpo_drag_coeff_recovery(&s->grid, app->basis, &s->local, &app->local, fpo->gamma,
    fpo->h, fpo->dhdv_surf, fpo->drag_coeff); 
  gkyl_calc_fpo_diff_coeff_recovery(&s->grid, app->basis, &s->local, &app->local, fpo->gamma,
    fpo->g, fpo->g_surf, fpo->dgdv_surf, fpo->d2gdv2_surf, fpo->diff_coeff); 

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

// updates the collision terms in the rhs
void
vm_species_fpo_rhs(gkyl_vlasov_app *app, const struct vm_species *s,
  struct vm_fpo_collisions *fpo, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  if (app->use_gpu) {
    wst = gkyl_wall_clock();
    
    // accumulate update due to collisions onto rhs
    gkyl_dg_updater_fpo_vlasov_advance_cu(fpo->coll_slvr, &s->local,
      fpo->drag_coeff, fpo->diff_coeff, fin, s->cflrate, rhs);
    
    app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
    
  } 
  else {
    wst = gkyl_wall_clock();
    
    // accumulate update due to collisions onto rhs
    gkyl_dg_updater_fpo_vlasov_advance(fpo->coll_slvr, &s->local,
      fpo->drag_coeff, fpo->diff_coeff, fin, s->cflrate, rhs);

    app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
  }
}

void 
vm_species_fpo_release(const struct gkyl_vlasov_app *app, const struct vm_fpo_collisions *fpo)
{
  gkyl_array_release(fpo->gamma);
  gkyl_array_release(fpo->h);
  gkyl_array_release(fpo->g);
  gkyl_array_release(fpo->h_surf);
  gkyl_array_release(fpo->g_surf);
  gkyl_array_release(fpo->dhdv_surf);
  gkyl_array_release(fpo->dgdv_surf);
  gkyl_array_release(fpo->d2gdv2_surf);

  gkyl_proj_maxwellian_pots_on_basis_release(fpo->pot_slvr);

  vm_species_moment_release(app, &fpo->moms);

  gkyl_array_release(fpo->drag_coeff);
  gkyl_array_release(fpo->diff_coeff);

  gkyl_dg_updater_fpo_vlasov_release(fpo->coll_slvr);
}
