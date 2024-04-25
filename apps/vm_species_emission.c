#include <assert.h>
#include <gkyl_vlasov_priv.h>

void
vm_species_emission_init(struct gkyl_vlasov_app *app, struct vm_emitting_wall *emit,
  int dir, enum gkyl_edge_loc edge, struct gkyl_range *ghost_r, void *ctx)
{
  struct vm_emission_ctx *params = ctx;
  emit->params = params;
  emit->num_species = params->num_species;
  emit->f_emit = mkarr(app->use_gpu, app->basis.num_basis, ghost_r->volume);
  emit->edge = edge;
  emit->dir = dir;
}

void
vm_species_emission_cross_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_emitting_wall *emit)
{
  int cdim = app->cdim;
  int vdim = app->vdim;

  emit->emit_skin_r = (emit->edge == GKYL_LOWER_EDGE) ? &s->lower_skin[emit->dir] : &s->upper_skin[emit->dir];
  emit->emit_ghost_r = (emit->edge == GKYL_LOWER_EDGE) ? &s->lower_ghost[emit->dir] : &s->upper_ghost[emit->dir];

  for (int i=0; i<emit->num_species; ++i) {
    emit->impact_species[i] = vm_find_species(app, emit->params->in_species[i]);
    
    emit->impact_skin_r[i] = (emit->edge == GKYL_LOWER_EDGE) ? &emit->impact_species[i]->lower_skin[emit->dir] : &emit->impact_species[i]->upper_skin[emit->dir];
    emit->impact_ghost_r[i] = (emit->edge == GKYL_LOWER_EDGE) ? &emit->impact_species[i]->lower_ghost[emit->dir] : &emit->impact_species[i]->upper_ghost[emit->dir];
    
    emit->yield[i] = mkarr(app->use_gpu, app->confBasis.num_basis, emit->impact_ghost_r[i]->volume);
    emit->spectrum[i] = mkarr(app->use_gpu, app->basis.num_basis, emit->emit_ghost_r->volume);

    emit->update[i] = gkyl_bc_emission_spectrum_new(emit->params->norm_type[i],
      emit->params->yield_type[i], emit->params->norm_params[i], emit->params->yield_params[i],
      emit->yield[i], emit->spectrum[i], emit->dir, emit->edge, cdim, vdim, 
      emit->impact_skin_r[i], emit->impact_ghost_r[i], &app->grid, app->use_gpu);
  }
}

/* void */
/* vm_species_emission_apply_bc() */
/* { */
/*   for (int i=0; i<emit->num_species; ++i) { */
/*     gkyl_bc_emission_spectrum_advance(); */
/*   } */
/* } */

/* void */
/* vm_species_emission_release() */
/* { */
  
/* } */
