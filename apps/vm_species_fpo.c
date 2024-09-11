#include <assert.h>
#include <gkyl_vlasov_priv.h>

#include <gkyl_array_rio.h>

static void
create_offsets(const int num_up_dirs, const int update_dirs[2], 
  const struct gkyl_range *range, const int idxc[GKYL_MAX_DIM], long offsets[9])
{
  
  // Check if we're at an upper or lower edge in each direction
  bool is_edge_upper[2], is_edge_lower[2];
  for (int i=0; i<num_up_dirs; ++i) {
    is_edge_lower[i] = idxc[update_dirs[i]] == range->lower[update_dirs[i]];
    is_edge_upper[i] = idxc[update_dirs[i]] == range->upper[update_dirs[i]];
  }

  // Construct the offsets *only* in the directions being updated.
  // No need to load the neighbors that are not needed for the update.
  int lower_offset[GKYL_MAX_DIM] = {0};
  int upper_offset[GKYL_MAX_DIM] = {0};
  for (int d=0; d<num_up_dirs; ++d) {
    int dir = update_dirs[d];
    lower_offset[dir] = -1 + is_edge_lower[d];
    upper_offset[dir] = 1 - is_edge_upper[d];
  }  

  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, lower_offset, upper_offset);
  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);
  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

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
 
  fpo->h_surf = mkarr(app->use_gpu, vdim*surf_basis.num_basis, s->local_ext.volume);
  fpo->g_surf = mkarr(app->use_gpu, vdim*surf_basis.num_basis, s->local_ext.volume); 
  fpo->dhdv_surf = mkarr(app->use_gpu, vdim*surf_basis.num_basis, s->local_ext.volume);
  fpo->dgdv_surf = mkarr(app->use_gpu, 2*vdim*surf_basis.num_basis, s->local_ext.volume); 
  fpo->d2gdv2_surf = mkarr(app->use_gpu, vdim*surf_basis.num_basis, s->local_ext.volume); 

  fpo->pot_slvr = gkyl_proj_maxwellian_pots_on_basis_new(&s->grid, &app->confBasis, &app->basis, app->poly_order+1, app->use_gpu);

  // allocate moments needed for FPO update
  vm_species_moment_init(app, s, &fpo->lte_moms, "LTEMoments");
  vm_species_moment_init(app, s, &fpo->moms, "FiveMoments");

  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = s->info.lower[d];
    v_bounds[d + vdim] = s->info.upper[d];
  }

  // initialize drag and diffusion coefficient arrays
  int num_surf_quad_nodes = (int)(pow(app->poly_order+1, pdim-1));
  fpo->drag_coeff = mkarr(app->use_gpu, vdim*app->basis.num_basis, s->local_ext.volume);
  fpo->drag_coeff_surf = mkarr(app->use_gpu, vdim*surf_basis.num_basis, s->local_ext.volume);
  fpo->sgn_drag_coeff_surf = mkarr(app->use_gpu, vdim*num_surf_quad_nodes, s->local_ext.volume); 
  fpo->const_sgn_drag_coeff_surf = mkarr(app->use_gpu, vdim, s->local_ext.volume);

  fpo->diff_coeff = mkarr(app->use_gpu, vdim*vdim*app->basis.num_basis, s->local_ext.volume);
  fpo->diff_coeff_surf = mkarr(app->use_gpu, 2*vdim*vdim*surf_basis.num_basis, s->local_ext.volume);

  // velocity space boundary corrections for momentum and energy
  fpo->fpo_moms = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
  fpo->boundary_corrections = mkarr(app->use_gpu, 2*(vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
  fpo->drag_diff_coeff_corrs = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);

  // Array of offsets for each phase space cell. Stored as:
  // 3-cell: (vx), (vy), (vx)
  // 9-cell: (vx-vy), (vx-vz), (vy-vz) 
  struct gkyl_array *offsets_ho = mk_int_arr(0, (vdim*3 + vdim*9), s->local_ext.volume);
  fpo->offsets = mk_int_arr(app->use_gpu, (vdim*3 + vdim*9), s->local_ext.volume);

  // Initializing the offset array here, need to find a better place to put this.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s->local);
  int idxc[GKYL_MAX_DIM];
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idxc);
    long linp = gkyl_range_idx(&s->local_ext, iter.idx);

    int *offsets_d = (int *)gkyl_array_fetch(offsets_ho, linp);

    long offsets_arr[9];
    for (int d1=0; d1<vdim; ++d1) {
      int dir1 = d1 + cdim;
      int update_dir[] = {dir1};
      create_offsets(1, update_dir, &s->local, idxc, offsets_arr);

      int idx[3][GKYL_MAX_DIM];
      for (int i=0; i<3; ++i) {
        offsets_d[d1*3 + i] = (int)offsets_arr[i];
        gkyl_range_inv_idx(&s->local, linp+offsets_arr[i], idx[i]);
      }

      long offsets_arr[9];
      for (int d2=d1+1; d2<vdim; ++d2) {
        int dir2 = d2 + cdim;
        int update_dirs[] = {dir1, dir2};
        create_offsets(2, update_dirs, &s->local, idxc, offsets_arr);

        for (int i=0; i<9; ++i) {
          offsets_d[3*vdim + (d1+d2-1) + i] = (int)offsets_arr[i];
        }
      }
    }
  }
  gkyl_array_copy(fpo->offsets, offsets_ho);

  const struct gkyl_mom_type* fpo_mom_type = gkyl_mom_fpo_vlasov_new(&app->confBasis, 
    &app->basis, &s->local, app->use_gpu);

  struct gkyl_mom_fpo_vlasov_auxfields fpo_mom_auxfields = { 
    .a = fpo->drag_coeff, .D = fpo->diff_coeff};
  gkyl_mom_fpo_vlasov_set_auxfields(fpo_mom_type, fpo_mom_auxfields);
  
  fpo->fpo_mom_calc = gkyl_mom_calc_new(&s->grid, fpo_mom_type, app->use_gpu);
  fpo->bcorr_calc = gkyl_mom_calc_bcorr_fpo_vlasov_new(&s->grid,
    &app->confBasis, &app->basis, &s->local, v_bounds, fpo->diff_coeff, app->use_gpu);

  fpo->coeffs_correct_calc = gkyl_fpo_coeffs_correct_new(app->use_gpu, &s->grid, &app->confBasis, &app->local);

  // initialize FPO updater
  fpo->coll_slvr = gkyl_dg_updater_fpo_vlasov_new(&s->grid, &app->basis, &s->local, app->use_gpu);
}

// computes drag coefficient and diffusion tensor
void
vm_species_fpo_drag_diff_coeffs(gkyl_vlasov_app *app, const struct vm_species *s,
  struct vm_fpo_collisions *fpo, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  // calculate needed moments
  vm_species_moment_calc(&fpo->lte_moms, s->local, app->local, fin);
  vm_species_moment_calc(&fpo->moms, s->local, app->local, fin);

  // calculate maxwellian potentials
  gkyl_proj_maxwellian_pots_on_basis_advance(fpo->pot_slvr, &s->local, &app->local, 
    fpo->lte_moms.marr, fpo->h, fpo->g, fpo->h_surf, fpo->g_surf,
    fpo->dhdv_surf, fpo->dgdv_surf, fpo->d2gdv2_surf);

  // calculate drag and diffusion coefficients
  gkyl_calc_fpo_drag_coeff_recovery(app->use_gpu, &s->grid, app->basis, &s->local, &app->local, 
    fpo->offsets, fpo->gamma,
    fpo->h, fpo->dhdv_surf, fpo->drag_coeff, fpo->drag_coeff_surf); 

  gkyl_calc_fpo_diff_coeff_recovery(&s->grid, app->basis, &s->local, &app->local, fpo->gamma,
    fpo->g, fpo->g_surf, fpo->dgdv_surf, fpo->d2gdv2_surf, 
    fpo->diff_coeff, fpo->diff_coeff_surf); 

  // Calculate corrections for momentum and energy conservation
  if (app->use_gpu) {
    // calculate volume corrections
    gkyl_mom_calc_advance_cu(fpo->fpo_mom_calc, &s->local, &app->local, fin, fpo->fpo_moms);

    // calculate boundary corrections
    gkyl_mom_calc_bcorr_advance_cu(fpo->bcorr_calc,
      &s->local, &app->local, fin, fpo->boundary_corrections);
  }
  else {
    // calculate volume corrections
    gkyl_mom_calc_advance(fpo->fpo_mom_calc, &s->local, &app->local, fin, fpo->fpo_moms);

    // calculate boundary corrections
    gkyl_mom_calc_bcorr_advance(fpo->bcorr_calc,
      &s->local, &app->local, fin, fpo->boundary_corrections);
  }

  // solve linear system for corrections and accumulate onto drag/diff coefficients
  gkyl_fpo_coeffs_correct_advance(fpo->coeffs_correct_calc, 
    &app->local, &s->local, fpo->fpo_moms, fpo->boundary_corrections, 
    fpo->moms.marr, fpo->drag_diff_coeff_corrs, 
    fpo->drag_coeff, fpo->drag_coeff_surf, 
    fpo->diff_coeff, fpo->diff_coeff_surf);

  gkyl_mom_calc_advance(fpo->fpo_mom_calc, &s->local, &app->local, fin, fpo->fpo_moms);
  gkyl_mom_calc_bcorr_advance(fpo->bcorr_calc,&s->local, &app->local, fin, fpo->boundary_corrections);

  // compute sign information at interfaces with corrected drag coefficient
  gkyl_calc_fpo_sgn_drag_coeff(app->basis, &s->local, fpo->drag_coeff_surf, 
    fpo->sgn_drag_coeff_surf, fpo->const_sgn_drag_coeff_surf);

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

// updates the collision terms in the rhs
void
vm_species_fpo_rhs(gkyl_vlasov_app *app, const struct vm_species *s,
  struct vm_fpo_collisions *fpo, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  wst = gkyl_wall_clock();
 
  // accumulate update due to collisions onto rhs
  gkyl_dg_updater_fpo_vlasov_advance(fpo->coll_slvr, &s->local,
    fpo->drag_coeff, fpo->drag_coeff_surf, 
    fpo->sgn_drag_coeff_surf, fpo->const_sgn_drag_coeff_surf,
    fpo->diff_coeff, fpo->diff_coeff_surf, fin, s->cflrate, rhs);

  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
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
