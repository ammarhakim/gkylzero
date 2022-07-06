#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_source_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_source *src)
{
  // ensure that boundary fluxes will be calculated if needed for source
  if (s->source_id == GKYL_BFLUX_SOURCE) {
    s->boundary_fluxes = true;
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
double
vm_species_source_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  // GKYL_BFLUX_SOURCE currently assumes 1X
  if (s->source_id == GKYL_BFLUX_SOURCE) {
    src->scale_factor = 1.0;
  }
  gkyl_array_accumulate(rhs, src->scale_factor, src->source);
  return 0;
}

void
vm_species_source_release(const struct gkyl_vlasov_app *app, const struct vm_source *src)
{
  gkyl_array_release(src->source);
  if (app->use_gpu)
    gkyl_array_release(src->source_host);

  gkyl_proj_on_basis_release(src->source_proj);
}
