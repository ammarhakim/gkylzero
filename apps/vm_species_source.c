#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_source_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_source *src)
{
  src->calc_bflux = false;
  if (s->source_id == GKYL_BFLUX_SOURCE) {
    src->calc_bflux = true;
    assert(s->info.source.source_length);
    assert(s->info.source.source_species);
    src->source_length = s->info.source.source_length;
    src->source_species = vm_find_species(app, s->info.source.source_species);
    src->source_species_idx = vm_find_species_idx(app, s->info.source.source_species);
    vm_species_bflux_init(app, src->source_species, &src->bflux); // boundary flux updater
    if (app->use_gpu) {
      src->scale_ptr = gkyl_cu_malloc((2 + app->vdim)*sizeof(double));
    }
    else {
      src->scale_ptr = gkyl_malloc((2 + app->vdim)*sizeof(double));
    }
  }

  // we need to ensure source has same shape as distribution function
  src->source = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  src->scale_factor = 1.0;
  
  vm_species_projection_init(app, s, s->info.source.projection, &src->proj_source);
}

void
vm_species_source_calc(gkyl_vlasov_app *app, struct vm_species *s, 
  struct vm_source *src, double tm)
{
  if (s->source_id) {
    vm_species_projection_calc(app, s, &src->proj_source, src->source, tm);
  }
}

// computes rhs of the boundary flux
void
vm_species_source_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_source *src, const struct gkyl_array *fin[], struct gkyl_array *rhs[])
{
  int species_idx;
  species_idx = vm_find_species_idx(app, species->info.name);
  // use boundary fluxes to scale source profile
  if (species->source_id == GKYL_BFLUX_SOURCE) {
    src->scale_factor = 0;
    double z[app->confBasis.num_basis];
    double red_mom[1];

    for (int d=0; d<app->cdim; ++d) {
      gkyl_array_reduce(src->scale_ptr, src->bflux.integ_moms[2*d].marr, GKYL_SUM);
      if (app->use_gpu) {
        gkyl_cu_memcpy(red_mom, src->scale_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
      }
      else {
        red_mom[0] = src->scale_ptr[0];
      }
      src->scale_factor += red_mom[0];
      gkyl_array_reduce(src->scale_ptr, src->bflux.integ_moms[2*d+1].marr, GKYL_SUM);
      if (app->use_gpu) {
        gkyl_cu_memcpy(red_mom, src->scale_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
      }
      else {
        red_mom[0] = src->scale_ptr[0];
      }
      src->scale_factor += red_mom[0];
    }
    src->scale_factor = src->scale_factor/src->source_length;
  }

  gkyl_array_accumulate(rhs[species_idx], src->scale_factor, src->source);

  // bflux calculation needs to be after source. The source uses bflux from the previous stage.
  if (src->calc_bflux) {
    vm_species_bflux_rhs(app, src->source_species, &src->bflux, fin[src->source_species_idx],
      rhs[src->source_species_idx]);
  }
}

void
vm_species_source_release(const struct gkyl_vlasov_app *app, const struct vm_source *src)
{
  gkyl_array_release(src->source);
  vm_species_projection_release(app, &src->proj_source);
  if (src->calc_bflux) {
    if (app->use_gpu) {
      gkyl_cu_free(src->scale_ptr);
    } 
    else {
      gkyl_free(src->scale_ptr);
    }    
    vm_species_bflux_release(app, &src->bflux);
  }
}
