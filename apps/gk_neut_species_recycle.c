#include <assert.h>
#include <gkyl_mom_canonical_pb.h>
#include <gkyl_gyrokinetic_priv.h>

void
gk_neut_species_recycle_init(struct gkyl_gyrokinetic_app *app, struct gk_recycle_wall *recyc,
  int dir, enum gkyl_edge_loc edge, void *ctx, struct gkyl_array *f0,
  struct gk_neut_species *s, bool use_gpu)
{
  struct gkyl_bc_emission_ctx *params = ctx;
  recyc->params = params;
  recyc->num_species = params->num_species;
  recyc->edge = edge;
  recyc->dir = dir;
  recyc->elastic = params->elastic;
  recyc->t_bound = params->t_bound;
  recyc->f0 = f0;

  int bdir = (recyc->edge == GKYL_LOWER_EDGE) ? 2*recyc->dir : 2*recyc->dir+1;

  recyc->emit_buff_r = &s->bflux.flux_r[bdir];
  recyc->emit_ghost_r = (recyc->edge == GKYL_LOWER_EDGE) ? &s->lower_ghost[recyc->dir] : &s->upper_ghost[recyc->dir];

  recyc->f0_flux_slvr = gkyl_ghost_surf_calc_new(&s->grid, s->eqn_vlasov, app->cdim, app->use_gpu);
  recyc->init_bflux_arr = mkarr(app->use_gpu, s->basis.num_basis, recyc->emit_buff_r->volume);

  if (app->use_gpu) {
    gkyl_ghost_surf_calc_advance_cu(recyc->f0_flux_slvr, &s->local_ext, recyc->f0, recyc->f0);
  } else {
    gkyl_ghost_surf_calc_advance(recyc->f0_flux_slvr, &s->local_ext, recyc->f0, recyc->f0);
  }
  gkyl_array_copy_range_to_range(recyc->init_bflux_arr, recyc->f0, recyc->emit_buff_r,
    recyc->emit_ghost_r);

  gkyl_ghost_surf_calc_release(recyc->f0_flux_slvr);
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

  recyc->emit_grid = &s->bflux.boundary_grid[bdir];
  recyc->emit_cbuff_r = &s->bflux.conf_r[bdir];
  recyc->emit_skin_r = (recyc->edge == GKYL_LOWER_EDGE) ? &s->lower_skin[recyc->dir] : &s->upper_skin[recyc->dir];
  recyc->buffer = (recyc->edge == GKYL_LOWER_EDGE) ? s->bc_buffer_lo_recyc : s->bc_buffer_up_recyc;
  
  recyc->f_emit = mkarr(app->use_gpu, s->basis.num_basis, recyc->emit_buff_r->volume);

  struct gkyl_array *proj_buffer = mkarr(false, s->basis.num_basis, recyc->emit_buff_r->volume);

  // Calculate the flux
  gkyl_bc_emission_flux_ranges(&recyc->emit_normal_r, recyc->dir + cdim, recyc->emit_buff_r,
    ghost, recyc->edge);

  recyc->init_flux = mkarr(app->use_gpu, app->basis.num_basis, recyc->emit_cbuff_r->volume);
  recyc->init_conf_grid = &s->bflux.conf_boundary_grid[bdir];
  
  struct gkyl_mom_canonical_pb_auxfields can_pb_inp = {.hamil = s->hamil};
  recyc->init_flux_slvr = gkyl_dg_updater_moment_new(recyc->emit_grid, &app->basis,
    &s->basis, recyc->emit_cbuff_r, &s->local_vel, recyc->emit_buff_r, s->model_id,
    &can_pb_inp, "M0", false, app->use_gpu);
  
  gkyl_dg_updater_moment_advance(recyc->init_flux_slvr, &recyc->emit_normal_r, recyc->emit_cbuff_r, recyc->init_bflux_arr,
				 recyc->init_flux);

  if (app->use_gpu) {
    recyc->mem_geo = gkyl_dg_bin_op_mem_cu_dev_new(recyc->emit_cbuff_r->volume, app->basis.num_basis);
  }
  else {
    recyc->mem_geo = gkyl_dg_bin_op_mem_new(recyc->emit_cbuff_r->volume, app->basis.num_basis);
  }
  
  // Initialize elastic component of emission
  if (recyc->elastic) {
    recyc->elastic_yield = mkarr(app->use_gpu, s->basis.num_basis, recyc->emit_buff_r->volume);
    recyc->elastic_update = gkyl_bc_emission_elastic_new(recyc->params->elastic_model,
      recyc->elastic_yield, recyc->dir, recyc->edge, cdim, vdim, s->info.mass, s->f->ncomp, recyc->emit_grid,
      recyc->emit_buff_r, app->poly_order, s->basis_on_dev, &s->basis, proj_buffer,
      app->use_gpu);
  }

  // Initialize inelastic emission spectrums
  for (int i=0; i<recyc->num_species; ++i) {

    const struct gkyl_emission_yield_constant *model = container_of(recyc->params->yield_model[i],
      struct gkyl_emission_yield_constant, yield);
    recyc->delta[i] = model->delta;
    recyc->spectrum[i] = mkarr(app->use_gpu, s->basis.num_basis, recyc->emit_buff_r->volume);
    
    recyc->impact_species[i] = gk_find_species(app, recyc->params->in_species[i]);
    recyc->impact_grid[i] = &recyc->impact_species[i]->bflux.boundary_grid[bdir];
    recyc->impact_conf_grid[i] = &recyc->impact_species[i]->bflux.conf_boundary_grid[bdir];
    recyc->impact_skin_r[i] = (recyc->edge == GKYL_LOWER_EDGE) ? &recyc->impact_species[i]->lower_skin[recyc->dir] : &recyc->impact_species[i]->upper_skin[recyc->dir];
    recyc->impact_ghost_r[i] = (recyc->edge == GKYL_LOWER_EDGE) ? &recyc->impact_species[i]->lower_ghost[recyc->dir] : &recyc->impact_species[i]->upper_ghost[recyc->dir];
    recyc->impact_buff_r[i] = &recyc->impact_species[i]->bflux.flux_r[bdir];
    recyc->impact_cbuff_r[i] = &recyc->impact_species[i]->bflux.conf_r[bdir];

    recyc->flux_slvr[i] = gkyl_dg_updater_moment_gyrokinetic_new(recyc->impact_grid[i], &app->basis,
      &recyc->impact_species[i]->basis, recyc->emit_cbuff_r, recyc->impact_species[i]->info.mass, recyc->impact_species[i]->vel_map,
      app->gk_geom, "M0", 0, app->use_gpu);
    
    recyc->flux[i] = mkarr(app->use_gpu, app->basis.num_basis, recyc->impact_cbuff_r[i]->volume);
    recyc->bflux_arr[i] = recyc->impact_species[i]->bflux.flux_arr[bdir];

    gkyl_bc_emission_flux_ranges(&recyc->impact_normal_r[i], recyc->dir + cdim, recyc->impact_buff_r[i],
      ghost, recyc->edge);
  }
  gkyl_array_release(proj_buffer);
}

void
gk_neut_species_recycle_apply_bc(struct gkyl_gyrokinetic_app *app, const struct gk_recycle_wall *recyc,
  const struct gk_neut_species *s, struct gkyl_array *fout)
{
  // Optional scaling of emission with time
  double t_scale = 1.0;
  /* if (recyc->t_bound) */
  /*   t_scale = sin(M_PI*tcurr/(2.0*recyc->t_bound)); */

  gkyl_array_clear(recyc->f_emit, 0.0); // Zero emitted distribution before beginning accumulate

  // Elastic emission contribution
  if (recyc->elastic) {
    gkyl_bc_emission_elastic_advance(recyc->elastic_update, recyc->emit_skin_r, recyc->buffer, fout,
      recyc->f_emit, recyc->elastic_yield, &s->basis);
  }
  // Inelastic emission contribution
  for (int i=0; i<recyc->num_species; ++i) {
    gkyl_array_clear(recyc->spectrum[i], 0.0);
    gkyl_array_accumulate(recyc->spectrum[i], recyc->delta[i], recyc->buffer);

    int species_idx;
    species_idx = gk_find_species_idx(app, recyc->impact_species[i]->info.name);
    
    gkyl_dg_updater_moment_gyrokinetic_advance(recyc->flux_slvr[i], &recyc->impact_normal_r[i],
      recyc->emit_cbuff_r, recyc->bflux_arr[i], recyc->flux[i]);

    /* const char *fmt = "recyc_flux_edge_%d.gkyl"; */
    /* int sz = gkyl_calc_strlen(fmt, recyc->edge); */
    /* char fileNm[sz+1]; // ensures no buffer overflow */
    /* snprintf(fileNm, sizeof fileNm, fmt, recyc->edge); */
    /* gkyl_grid_sub_array_write(recyc->impact_conf_grid[i], recyc->emit_cbuff_r, 0, recyc->flux[i], fileNm); */

    gkyl_array_set_range_to_range(fout, t_scale, recyc->f_emit, recyc->emit_ghost_r,
      recyc->emit_buff_r);

    gkyl_dg_div_op_range(recyc->mem_geo, app->basis, 0, recyc->flux[i], 0, recyc->flux[i],
    			 0, recyc->init_flux, recyc->emit_cbuff_r);

    gkyl_dg_mul_conf_phase_op_accumulate_range(&app->basis, &s->basis, recyc->f_emit, 1.0,
     recyc->flux[i], recyc->spectrum[i], recyc->impact_cbuff_r[i], recyc->emit_buff_r);
    
  }
  const char *fmt = "recyc_f_emit_edge_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, recyc->edge);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, recyc->edge);
  gkyl_grid_sub_array_write(recyc->emit_grid, recyc->emit_buff_r, 0, recyc->f_emit, fileNm);

  gkyl_array_set_range_to_range(fout, t_scale, recyc->f_emit, recyc->emit_ghost_r,
    recyc->emit_buff_r);
}

void
gk_neut_species_recycle_release(const struct gk_recycle_wall *recyc)
{
  gkyl_array_release(recyc->f_emit);
  gkyl_array_release(recyc->init_flux);
  gkyl_array_release(recyc->init_bflux_arr);
  gkyl_dg_updater_moment_release(recyc->init_flux_slvr);
  gkyl_dg_bin_op_mem_release(recyc->mem_geo);
  if (recyc->elastic) {
    gkyl_array_release(recyc->elastic_yield);
    gkyl_bc_emission_elastic_release(recyc->elastic_update);
  }
  for (int i=0; i<recyc->num_species; ++i) {
    gkyl_array_release(recyc->spectrum[i]);
    gkyl_array_release(recyc->flux[i]);
    gkyl_dg_updater_moment_release(recyc->flux_slvr[i]);
  }
}
