#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_source *src)
{
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

void
gk_species_source_calc(gkyl_gyrokinetic_app *app, struct gk_species *species, double tm)
{
  if (species->source_id) {
    gkyl_proj_on_basis_advance(species->src.source_proj, tm, &species->local_ext, species->src.source_host);
    if (app->use_gpu) // note: source_host is same as source when not on GPUs
      gkyl_array_copy(species->src.source, species->src.source_host);
  }
}

// computes rhs of the boundary flux
void
gk_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_source *src, const struct gkyl_array *fin[], struct gkyl_array *rhs[])
{
  int species_idx;
  species_idx = gk_find_species_idx(app, species->info.name);

  gkyl_array_accumulate(rhs[species_idx], src->scale_factor, src->source);
}

void
gk_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src)
{
  gkyl_array_release(src->source);
  if (app->use_gpu) 
    gkyl_array_release(src->source_host);

  gkyl_proj_on_basis_release(src->source_proj);
}
