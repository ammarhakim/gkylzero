#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_source_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_source *src)
{
  if (s->source_id == GKYL_BFLUX_SOURCE) {
    assert(s->info.source.source_length);
    assert(s->info.source.source_species);
    src->source_length = s->info.source.source_length;
    src->source_species = vm_find_species(app, s->info.source.source_species);
    src->source_species->calc_bflux = true; // ensure that boundary fluxes will be calculated
    if (app->use_gpu)
      src->scale_ptr = gkyl_cu_malloc(sizeof(double));
    else
      src->scale_ptr = gkyl_malloc(sizeof(double));
  }
  
  // we need to ensure source has same shape as distribution function
  src->source = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  src->scale_factor = 1.0;
  
  src->source_host = src->source;
  if (app->use_gpu)
    src->source_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);

    src->source_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp)    {
      .grid = &s->grid,
      .basis = &app->basis,
      .qtype = GKYL_GAUSS_QUAD,
      .num_quad = app->basis.poly_order+1,
      .num_ret_vals = 1,
      .eval = s->info.source.profile,
      .ctx = s->info.source.ctx
    }
  );
}

// computes rhs of the boundary flux
void
vm_species_source_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  double red_mom[app->vdim+2];
  // use boundary fluxes to scale source profile
  if (species->source_id == GKYL_BFLUX_SOURCE) {
    src->scale_factor = 0;
    double z[app->confBasis.num_basis];
    double red_mom[1];
    
    for (int d=0; d<app->cdim; ++d) {
      gkyl_array_reduce(src->scale_ptr, src->source_species->bflux.integ_moms[2*d].marr, GKYL_SUM);
      if (app->use_gpu)
        gkyl_cu_memcpy(red_mom, src->scale_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
      else
        red_mom[0] = src->scale_ptr[0];
      src->scale_factor += red_mom[0];
      gkyl_array_reduce(src->scale_ptr, src->source_species->bflux.integ_moms[2*d+1].marr, GKYL_SUM);
      if (app->use_gpu)
        gkyl_cu_memcpy(red_mom, src->scale_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
      else
        red_mom[0] = src->scale_ptr[0];
      src->scale_factor += red_mom[0];
    }
    src->scale_factor = src->scale_factor/src->source_length;
  }
  gkyl_array_accumulate(rhs, src->scale_factor, src->source);
}

void
vm_species_source_release(const struct gkyl_vlasov_app *app, const struct vm_source *src)
{
  gkyl_array_release(src->source);
  if (app->use_gpu) {
    gkyl_array_release(src->source_host);
    gkyl_cu_free(src->scale_ptr);
  } else {
    gkyl_free(src->scale_ptr);
  }
  gkyl_proj_on_basis_release(src->source_proj);
}
