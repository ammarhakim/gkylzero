#include <string.h>
#include <gkyl_bc_emission.h>

struct gkyl_bc_emission_ctx*
gkyl_bc_emission_new(int num_species, double t_bound, bool elastic, struct gkyl_spectrum_model *spectrum_model[], struct gkyl_yield_model *yield_model[], struct gkyl_elastic_model *elastic_model, char in_species[][128])
{
  struct gkyl_bc_emission_ctx *ctx = gkyl_malloc(sizeof(struct gkyl_bc_emission_ctx));
  
  ctx->num_species = num_species;
  ctx->t_bound = t_bound;
  ctx->elastic = elastic;
  for (int i=0; i<num_species; ++i) {
    ctx->spectrum_model[i] = gkyl_spectrum_model_acquire(spectrum_model[i]);
    ctx->yield_model[i] = gkyl_yield_model_acquire(yield_model[i]);
    strcpy(ctx->in_species[i], in_species[i]);
  }
  ctx->elastic_model = gkyl_elastic_model_acquire(elastic_model);

  return ctx;
}

void gkyl_bc_emission_release(struct gkyl_bc_emission_ctx *ctx)
{
  for (int i=0; i<ctx->num_species; ++i) {
    gkyl_spectrum_model_release(ctx->spectrum_model[i]);
    gkyl_yield_model_release(ctx->yield_model[i]);
  }
  gkyl_elastic_model_release(ctx->elastic_model);
  // Release updater memory.
  gkyl_free(ctx);
}
