#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_neut_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
  struct gk_source *src)
{
  // we need to ensure source has same shape as distribution function
  src->source = mkarr(app->use_gpu, app->neut_basis.num_basis, s->local_ext.volume);
  src->source_id = s->source_id;
  src->source_host = src->source;
  if (app->use_gpu)
    src->source_host = mkarr(false, app->neut_basis.num_basis, s->local_ext.volume);

  src->write_source = s->info.source.write_source; // optional flag to write out source

  src->num_sources = s->info.source.num_sources;
  for (int k=0; k<s->info.source.num_sources; k++)
    gk_neut_species_projection_init(app, s, s->info.source.projection[k], &src->proj_source[k]);
}

void
gk_neut_species_source_calc(gkyl_gyrokinetic_app *app, const struct gk_neut_species *s, 
  struct gk_source *src, double tm)
{
  struct gkyl_array *source_tmp = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  for (int k=0; k<s->info.source.num_sources; k++) {
    gk_neut_species_projection_calc(app, s, &src->proj_source[k], source_tmp, tm);
    gkyl_array_accumulate(src->source, 1., source_tmp);
  }
  gkyl_array_release(source_tmp);
}

// computes rhs of the source
void
gk_neut_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkyl_array_accumulate(rhs, 1.0, src->source);
}

void
gk_neut_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src)
{
  for (int k=0; k<src->num_sources; k++)
    gk_neut_species_projection_release(app, &src->proj_source[k]);
  gkyl_array_release(src->source);
  if (app->use_gpu) 
    gkyl_array_release(src->source_host);
}
