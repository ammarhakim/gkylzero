#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_neut_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
  struct gk_neut_source *src)
{
  // we need to ensure source has same shape as distribution function
  src->source = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  src->source_id = s->source_id;
  src->source_host = src->source;
  if (app->use_gpu)
    src->source_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);

  src->write_source = s->info.source.write_source; // optional flag to write out source

  gk_neut_species_projection_init(app, s, s->info.source.projection, &src->proj_source);
}

void
gk_neut_species_source_calc(gkyl_gyrokinetic_app *app, const struct gk_neut_species *s, 
  struct gk_neut_source *src, double tm)
{
  gk_neut_species_projection_calc(app, s, &src->proj_source, src->source, tm);
}

// computes rhs of the source
void
gk_neut_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_neut_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkyl_array_accumulate(rhs, 1.0, src->source);
}

void
gk_neut_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_neut_source *src)
{
  gk_neut_species_projection_release(app, &src->proj_source);
  gkyl_array_release(src->source);
  if (app->use_gpu) 
    gkyl_array_release(src->source_host);
}
