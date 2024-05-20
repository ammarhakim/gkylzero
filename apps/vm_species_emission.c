#include <assert.h>
#include <gkyl_vlasov_priv.h>

void
vm_species_emission_init(struct gkyl_vlasov_app *app, struct vm_emitting_wall *emit,
  int dir, enum gkyl_edge_loc edge, void *ctx, bool use_gpu)
{
  struct vm_emission_ctx *params = ctx;
  emit->params = params;
  emit->num_species = params->num_species;
  emit->edge = edge;
  emit->dir = dir;
}

void
vm_species_emission_cross_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_emitting_wall *emit)
{
  int cdim = app->cdim;
  int vdim = app->vdim;
  int bdir = (emit->edge == GKYL_LOWER_EDGE) ? 2*emit->dir : 2*emit->dir+1;

  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    ghost[cdim+d] = 0;
  }

  emit->emit_grid = &s->bflux.boundary_grid[bdir];
  emit->emit_buff_r = &s->bflux.flux_r[bdir];
  emit->emit_ghost_r = (emit->edge == GKYL_LOWER_EDGE) ? &s->lower_ghost[emit->dir] : &s->upper_ghost[emit->dir];

  for (int i=0; i<emit->num_species; ++i) {
    emit->impact_species[i] = vm_find_species(app, emit->params->in_species[i]);
    emit->impact_grid[i] = &emit->impact_species[i]->bflux.boundary_grid[bdir];

    emit->flux_slvr[i] = gkyl_dg_updater_moment_new(emit->impact_grid[i], &app->confBasis,
      &app->basis, NULL, NULL, emit->impact_species[i]->model_id, 0, "Integrated", 1,
      emit->impact_species[i]->info.mass, app->use_gpu);

    emit->impact_skin_r[i] = (emit->edge == GKYL_LOWER_EDGE) ? &emit->impact_species[i]->lower_skin[emit->dir] : &emit->impact_species[i]->upper_skin[emit->dir];
    emit->impact_ghost_r[i] = (emit->edge == GKYL_LOWER_EDGE) ? &emit->impact_species[i]->lower_ghost[emit->dir] : &emit->impact_species[i]->upper_ghost[emit->dir];
    emit->impact_buff_r[i] = &emit->impact_species[i]->bflux.flux_r[bdir];
    emit->impact_cbuff_r[i] = &emit->impact_species[i]->bflux.conf_r[bdir];

    emit->yield[i] = mkarr(app->use_gpu, app->basis.num_basis, emit->impact_buff_r[i]->volume);
    emit->spectrum[i] = mkarr(app->use_gpu, app->basis.num_basis, emit->emit_buff_r->volume);
    emit->weight[i] = mkarr(app->use_gpu, app->confBasis.num_basis,
      emit->impact_cbuff_r[i]->volume);
    emit->flux[i] = mkarr(app->use_gpu, app->confBasis.num_basis, emit->impact_cbuff_r[i]->volume);
    emit->bflux_arr[i] = emit->impact_species[i]->bflux.flux_arr[bdir];
    emit->k[i] = mkarr(app->use_gpu, app->confBasis.num_basis, emit->impact_cbuff_r[i]->volume);

    gkyl_bc_emission_flux_ranges(&emit->impact_normal_r[i], emit->dir, emit->impact_ghost_r[i],
      ghost, emit->edge);
    
    emit->update[i] = gkyl_bc_emission_spectrum_new(emit->params->norm_type[i],
      emit->params->yield_type[i], emit->params->norm_params[i], emit->params->yield_params[i],
      emit->yield[i], emit->spectrum[i], emit->dir, emit->edge, cdim, vdim, 
      emit->impact_buff_r[i], emit->impact_ghost_r[i], emit->impact_grid[i], app->poly_order,
      &app->basis, app->use_gpu);
  }
}

void
vm_species_emission_rhs(struct gkyl_vlasov_app *app, struct vm_emitting_wall *emit, struct gkyl_array *rhs[])
{
  for (int i=0; i<emit->num_species; ++i) {
    int species_idx;
    species_idx = vm_find_species_idx(app, emit->impact_species[i]->info.name);
    gkyl_dg_updater_moment_advance(emit->flux_slvr[i], emit->impact_buff_r[i],
      emit->impact_cbuff_r[i], emit->bflux_arr[i], emit->flux[i]);
    
    gkyl_bc_emission_spectrum_advance(emit->update[i], emit->impact_skin_r[i], emit->impact_buff_r[i], emit->impact_cbuff_r[i], &emit->emit_ghost_r[i], rhs[species_idx], emit->yield[i], emit->spectrum[i], emit->weight[i], emit->flux[i], emit->k[i]);
  }
}

void
vm_species_emission_apply_bc()
{
  int i = 1;
  /* for (int i=0; i<emit->num_species; ++i) { */
  /*   gkyl_bc_emission_spectrum_advance(emit->update[i], ); */
  /* } */
}

/* void */
/* vm_species_emission_release() */
/* { */
  
/* } */
