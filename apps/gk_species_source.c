#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src)
{
  int vdim = app->vdim;
  // we need to ensure source has same shape as distribution function
  src->source = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  src->source_id = s->source_id;
  src->source_host = src->source;
  if (app->use_gpu)
    src->source_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);

  src->write_source = s->info.source.write_source; // optional flag to write out source

  src->num_sources = s->info.source.num_sources;
  for (int k=0; k<s->info.source.num_sources; k++)
    gk_species_projection_init(app, s, s->info.source.projection[k], &src->proj_source[k]);

  // Initialize moments
  // allocate data for integrated moments
  gk_species_moment_init(app, s, &s->src.integ_moms, "Integrated");

  // allocate data for diagnostic moments
  int ndm = s->info.num_diag_moments;
  s->src.moms = gkyl_malloc(sizeof(struct gk_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    gk_species_moment_init(app, s, &s->src.moms[m], s->info.diag_moments[m]);

  if (app->use_gpu) 
    s->src.red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
  // allocate dynamic-vector to store all-reduced integrated moments 
  s->src.integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  s->src.is_first_integ_write_call = true;

}

void
gk_species_source_calc(gkyl_gyrokinetic_app *app, const struct gk_species *s, 
  struct gk_source *src, double tm)
{
  struct gkyl_array *source_tmp = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  for (int k=0; k<s->info.source.num_sources; k++) {
    gk_species_projection_calc(app, s, &src->proj_source[k], source_tmp, tm);
    gkyl_array_accumulate(src->source, 1., source_tmp);
  }
  gkyl_array_release(source_tmp);
}

// computes rhs of the source
void
gk_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkyl_array_accumulate(rhs, 1.0, src->source);
}

void
gk_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src)
{
  for (int k=0; k<src->num_sources; k++)
    gk_species_projection_release(app, &src->proj_source[k]);
  gkyl_array_release(src->source);
  if (app->use_gpu) 
    gkyl_array_release(src->source_host);
}

void
gkyl_gyrokinetic_app_calc_integrated_source_mom(gkyl_gyrokinetic_app* app, double tm)
{
  int vdim = app->vdim;
  double avals[2+vdim], avals_global[2+vdim];

  struct timespec wst = gkyl_wall_clock();

  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];

    if (s->src.source_id) {
      struct timespec wst = gkyl_wall_clock();

      gk_species_moment_calc(&s->src.integ_moms, s->local, app->local, s->src.source); 

      // reduce to compute sum over whole domain, append to diagnostics
      if (app->use_gpu) {
        gkyl_array_reduce_range(s->src.red_integ_diag, s->src.integ_moms.marr, GKYL_SUM, &(app->local));
        gkyl_cu_memcpy(avals, s->src.red_integ_diag, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
      }
      else {
        gkyl_array_reduce_range(avals, s->src.integ_moms.marr_host, GKYL_SUM, &(app->local));
      }

      gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, avals, avals_global);
      gkyl_dynvec_append(s->src.integ_diag, tm, avals_global);

      app->stat.mom_tm += gkyl_time_diff_now_sec(wst);
      app->stat.nmom += 1;
    }
  }

  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}
