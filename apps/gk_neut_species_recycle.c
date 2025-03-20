#include <assert.h>
#include <gkyl_mom_canonical_pb.h>
#include <gkyl_gyrokinetic_priv.h>

void
gk_neut_species_recycle_init(struct gkyl_gyrokinetic_app *app, struct gk_recycle_wall *recyc,
  int dir, enum gkyl_edge_loc edge, struct gkyl_gyrokinetic_emission_inp *params, struct gkyl_array *f0,
  struct gk_neut_species *s, bool use_gpu)
{
  recyc->params = params;
  recyc->num_species = params->num_species;
  recyc->edge = edge;
  recyc->dir = dir;

  int bdir = (recyc->edge == GKYL_LOWER_EDGE) ? 2*recyc->dir : 2*recyc->dir+1;
  int d = recyc->dir;
  int e = (recyc->edge == GKYL_LOWER_EDGE) ? 0 : 1;
  struct gkyl_range *skin_r = e==0? &s->lower_skin[d] : &s->upper_skin[d];
  struct gkyl_range *ghost_r = e==0? &s->lower_ghost[d] : &s->upper_ghost[d];

  recyc->emit_buff_r = &s->bflux.flux_r[bdir];
  recyc->emit_ghost_r = (recyc->edge == GKYL_LOWER_EDGE) ? &s->lower_ghost[d] : &s->upper_ghost[d];

  // Calculate flux associated with unit Maxwellian projected in f0.
  // Store in init_bflux_arr.
  recyc->f0_flux_slvr[bdir] = gkyl_boundary_flux_new(recyc->dir, e, &s->grid,
    skin_r, ghost_r, s->eqn_vlasov, true, app->use_gpu);  
  recyc->init_bflux_arr = mkarr(app->use_gpu, s->basis.num_basis, recyc->emit_buff_r->volume);
  gkyl_boundary_flux_advance(recyc->f0_flux_slvr[bdir], f0, f0);
  gkyl_array_copy_range_to_range(recyc->init_bflux_arr, f0, recyc->emit_buff_r,
    recyc->emit_ghost_r);
}

void
gk_neut_species_recycle_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_recycle_wall *recyc)
{
  int cdim = app->cdim;
  int vdim = app->vdim+1; // from gk_neut_species
  int bdir = (recyc->edge == GKYL_LOWER_EDGE) ? 2*recyc->dir : 2*recyc->dir+1;
 
  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    ghost[cdim+d] = 0;
  }

  // This buffer contains the projection of the unit Maxwellian.
  recyc->buffer = (recyc->edge == GKYL_LOWER_EDGE) ? s->bc_buffer_lo_recyc : s->bc_buffer_up_recyc;

  // Define necessary grid, ranges, and array for calculating the desired Maxwellian for ghost.
  recyc->emit_grid = &s->bflux.boundary_grid[bdir];
  recyc->emit_cbuff_r = &s->bflux.conf_r[bdir];
  recyc->emit_skin_r = (recyc->edge == GKYL_LOWER_EDGE) ? &s->lower_skin[recyc->dir] : &s->upper_skin[recyc->dir];
  recyc->f_emit = mkarr(app->use_gpu, s->basis.num_basis, recyc->emit_buff_r->volume);

  // For writing diagnostics, if needed.
  recyc->f_diag = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
  recyc->diag_out = mkarr(app->use_gpu, app->basis.num_basis, app->local.volume);
  recyc->diag_out_ho = recyc->diag_out;
  if (app->use_gpu) {
    recyc->diag_out_ho = mkarr(false, app->basis.num_basis, app->local.volume);
  }
   
  // Calculate the flux associated with unit-density Maxwellian.
  gkyl_bc_emission_flux_ranges(&recyc->emit_normal_r, recyc->dir + cdim, recyc->emit_buff_r,
    ghost, recyc->edge);
  
  recyc->init_flux = mkarr(app->use_gpu, app->basis.num_basis, recyc->emit_cbuff_r->volume);
  recyc->emit_flux = mkarr(app->use_gpu, app->basis.num_basis, recyc->emit_cbuff_r->volume);

  struct gkyl_mom_canonical_pb_auxfields can_pb_inp = {.hamil = s->hamil};
  recyc->init_flux_slvr = gkyl_dg_updater_moment_new(recyc->emit_grid, &app->basis,
    &s->basis, recyc->emit_cbuff_r, &s->local_vel, recyc->emit_buff_r, s->model_id,
    &can_pb_inp, "M0", false, app->use_gpu);
  
  gkyl_dg_updater_moment_advance(recyc->init_flux_slvr, &recyc->emit_normal_r, recyc->emit_cbuff_r, recyc->init_bflux_arr,
				 recyc->init_flux);

  // Define memory for div bin op for calculating correct scaling factor.
  if (app->use_gpu) {
    recyc->mem_geo = gkyl_dg_bin_op_mem_cu_dev_new(recyc->emit_cbuff_r->volume, app->basis.num_basis);
  }
  else {
    recyc->mem_geo = gkyl_dg_bin_op_mem_new(recyc->emit_cbuff_r->volume, app->basis.num_basis);
  }
  
  // Initialize inelastic emission spectrums
  for (int i=0; i<recyc->num_species; ++i) {
    
    recyc->impact_species[i] = gk_find_species(app, recyc->params->in_species[i]);
    struct gk_species *gks = recyc->impact_species[i];

    recyc->frac = recyc->params->rec_frac;
    recyc->spectrum[i] = mkarr(app->use_gpu, s->basis.num_basis, recyc->emit_buff_r->volume);

    recyc->impact_grid[i] = &gks->bflux_solver.boundary_grid[bdir];
    recyc->impact_buff_r[i] = &gks->bflux_solver.flux_r[bdir];
    recyc->impact_cbuff_r[i] = &gks->bflux_solver.conf_r[bdir];

    recyc->flux_slvr[i] = gkyl_dg_updater_moment_gyrokinetic_new(recyc->impact_grid[i], &app->basis,
      &gks->basis, recyc->emit_cbuff_r, gks->info.mass, gks->info.charge, gks->vel_map,
      app->gk_geom, app->field->phi_smooth, "M0", 0, app->use_gpu);
    
    recyc->flux[i] = mkarr(app->use_gpu, app->basis.num_basis, recyc->impact_cbuff_r[i]->volume);    
    recyc->bflux_arr[i] = gks->bflux_solver.flux_arr[bdir];

    gkyl_bc_emission_flux_ranges(&recyc->impact_normal_r[i], recyc->dir + cdim, recyc->impact_buff_r[i],
      ghost, recyc->edge);
  }
}

void
gk_neut_species_recycle_apply_bc(struct gkyl_gyrokinetic_app *app, const struct gk_recycle_wall *recyc,
  const struct gk_neut_species *s, struct gkyl_array *fout)
{
  gkyl_array_clear(recyc->f_emit, 0.0); // Zero emitted distribution before beginning accumulate

  // Inelastic emission contribution
  // This relies on the calculation of the ion flux (bflux_arr)
  // from gk_species_bflux_rhs_solver.
  for (int i=0; i<recyc->num_species; ++i) {
    // Clear array from previous step and copy unit-density Maxwellian from buffer.
    gkyl_array_clear(recyc->spectrum[i], 0.0);
    gkyl_array_accumulate(recyc->spectrum[i], recyc->frac, recyc->buffer);

    int species_idx;
    species_idx = gk_find_species_idx(app, recyc->impact_species[i]->info.name);

    // Calculate M0 moment of ion flux.
    gkyl_dg_updater_moment_gyrokinetic_advance(recyc->flux_slvr[i], &recyc->impact_normal_r[i],
      recyc->emit_cbuff_r, recyc->bflux_arr[i], recyc->flux[i]);

    // Calculate scaling factor from ratio of ion flux to unit-density flux.
    gkyl_dg_div_op_range(recyc->mem_geo, app->basis, 0, recyc->flux[i], 0, recyc->flux[i],
    			 0, recyc->init_flux, recyc->emit_cbuff_r);

    gkyl_dg_mul_conf_phase_op_accumulate_range(&app->basis, &s->basis, recyc->f_emit, 1.0,
     recyc->flux[i], recyc->spectrum[i], recyc->impact_cbuff_r[i], recyc->emit_buff_r);    
  }

  gkyl_array_set_range_to_range(fout, 1.0, recyc->f_emit, recyc->emit_ghost_r,
    recyc->emit_buff_r);
}

void
gk_neut_species_recycle_write_flux(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
				   struct gk_recycle_wall *recyc, double tm, int frame)
{
  // Output boundary flux from ions and neutral ghost cells
  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }
  );

  const char *edge = (recyc->edge == GKYL_LOWER_EDGE)? "lower" : "upper";
  struct gkyl_range *cskin_r = (recyc->edge == GKYL_LOWER_EDGE) ? &app->lower_skin[recyc->dir] : &app->upper_skin[recyc->dir];
  int bdir = (recyc->edge == GKYL_LOWER_EDGE) ? 2*recyc->dir : 2*recyc->dir+1;
  
  for (int i=0; i<recyc->num_species; ++i) {

    const char *fmt = "%s-%s_recyc_flux_%s_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, recyc->impact_species[i]->info.name,
      edge, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, recyc->impact_species[i]->info.name,
      edge, frame);
    
    gkyl_dg_updater_moment_gyrokinetic_advance(recyc->flux_slvr[i], &recyc->impact_normal_r[i],
      recyc->emit_cbuff_r, recyc->bflux_arr[i], recyc->flux[i]);
    gkyl_array_clear(recyc->diag_out, 0.0);
    gkyl_array_copy_range_to_range(recyc->diag_out, recyc->flux[i], cskin_r,
      recyc->emit_cbuff_r);

    if (app->use_gpu) {
      gkyl_array_clear(recyc->diag_out_ho, 0.0);
      gkyl_array_copy(recyc->diag_out_ho, recyc->diag_out);
    }

    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(app->comm, &app->grid, &app->local,
      mt, recyc->diag_out_ho, fileNm);
    app->stat.diag_io_tm += gkyl_time_diff_now_sec(wtm);
  }

  const char *fmt = "%s-%s_recyc_flux_%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->info.name,
    edge, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, s->info.name,
    edge, frame);
  
  // calc neut moment
  gkyl_array_clear(recyc->f_diag, 0.0);
  gkyl_array_copy_range_to_range(recyc->f_diag, recyc->f_emit, recyc->emit_skin_r,
    recyc->emit_buff_r);
  gkyl_boundary_flux_advance(recyc->f0_flux_slvr[bdir], recyc->f_diag, recyc->f_diag);
  gkyl_array_copy_range_to_range(recyc->init_bflux_arr, recyc->f_diag, recyc->emit_buff_r,
    recyc->emit_ghost_r);
  gkyl_dg_updater_moment_advance(recyc->init_flux_slvr, &recyc->emit_normal_r, recyc->emit_cbuff_r, recyc->init_bflux_arr,
  				 recyc->emit_flux);

  gkyl_array_clear(recyc->diag_out, 0.0);
  gkyl_array_copy_range_to_range(recyc->diag_out, recyc->emit_flux, cskin_r,
    recyc->emit_cbuff_r);

  if (app->use_gpu) {
    gkyl_array_clear(recyc->diag_out_ho, 0.0);
    gkyl_array_copy(recyc->diag_out_ho, recyc->diag_out);
  }
  
  struct timespec wtm = gkyl_wall_clock();
  gkyl_comm_array_write(app->comm, &app->grid, &app->local,
      mt, recyc->diag_out_ho, fileNm);
  app->stat.diag_io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.n_diag_io += 1;

  gk_array_meta_release(mt); 
}

void
gk_neut_species_recycle_release(const struct gkyl_gyrokinetic_app *app, const struct gk_recycle_wall *recyc)
{
  int bdir = (recyc->edge == GKYL_LOWER_EDGE) ? 2*recyc->dir : 2*recyc->dir+1;
  gkyl_array_release(recyc->f_emit);
  gkyl_array_release(recyc->init_flux);
  gkyl_array_release(recyc->init_bflux_arr);
  gkyl_array_release(recyc->f_diag);
  gkyl_array_release(recyc->emit_flux);
  gkyl_array_release(recyc->diag_out);
  gkyl_dg_updater_moment_release(recyc->init_flux_slvr);
  gkyl_dg_bin_op_mem_release(recyc->mem_geo);
  gkyl_boundary_flux_release(recyc->f0_flux_slvr[bdir]);

  if (app->use_gpu) {
    gkyl_array_release(recyc->diag_out_ho);
  }
  for (int i=0; i<recyc->num_species; ++i) {
    gkyl_array_release(recyc->spectrum[i]);
    gkyl_array_release(recyc->flux[i]);
    gkyl_dg_updater_moment_release(recyc->flux_slvr[i]);
  }
}
