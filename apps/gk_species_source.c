#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void
gk_species_source_write_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
}

void
gk_species_source_write_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = gks->basis.id
    }
  );

  // Write out the source distribution function
  const char *fmt = "%s-%s_source_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(gks->src.source_host, gks->src.source);
  }

  struct timespec wtm = gkyl_wall_clock();
  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->src.source_host, fileNm);
  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.n_io += 1;

  gk_array_meta_release(mt); 
}

void
gk_species_source_write_init_only(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gk_species_source_write_enabled(app, gks, tm, frame);
  gks->src.write_func = gk_species_source_write_disabled;
}

void
gk_species_source_write_mom_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
}

void
gk_species_source_write_mom_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }
  );

  for (int m=0; m<gks->src.num_diag_mom; ++m) {
    gk_species_moment_calc(&gks->src.moms[m], gks->local, app->local, gks->src.source);
    app->stat.n_mom += 1;

    const char *fmt = "%s-%s_source_%s_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name,
      gks->info.source.diagnostics.diag_moments[m], frame);
    char fileNm[sz+1]; // Ensures no buffer overflow.
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name,
      gks->info.source.diagnostics.diag_moments[m], frame);

    // Rescale moment by inverse of Jacobian. 
    // For Maxwellian and bi-Maxwellian moments, we only need to re-scale
    // the density (the 0th component).
    gkyl_dg_div_op_range(gks->moms[m].mem_geo, app->basis, 
      0, gks->src.moms[m].marr, 0, gks->src.moms[m].marr, 0, 
      app->gk_geom->jacobgeo, &app->local);      

    if (app->use_gpu) {
      gkyl_array_copy(gks->src.moms[m].marr_host, gks->src.moms[m].marr);
    }

    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
      gks->src.moms[m].marr_host, fileNm);
    app->stat.diag_io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.n_diag_io += 1;
  }
  gk_array_meta_release(mt); 

  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += 1;
}

void
gk_species_source_write_mom_init_only(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gk_species_source_write_mom_enabled(app, gks, tm, frame);
  gks->src.write_mom_func = gk_species_source_write_mom_disabled;
}

void
gk_species_source_calc_integrated_mom_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
}

void
gk_species_source_calc_integrated_mom_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  double tm_prev = gkyl_dynvec_getlast_tm(gks->src.integ_diag);

  int num_mom = gks->src.integ_moms.num_mom;
  double avals_global[num_mom];

  gk_species_moment_calc(&gks->src.integ_moms, gks->local, app->local, gks->src.source); 
  app->stat.n_mom += 1;

  // Reduce to compute sum over whole domain, append to diagnostics
  gkyl_array_reduce_range(gks->src.red_integ_diag, gks->src.integ_moms.marr, GKYL_SUM, &app->local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom, 
    gks->src.red_integ_diag, gks->src.red_integ_diag_global);
  if (app->use_gpu) {
    gkyl_cu_memcpy(avals_global, gks->src.red_integ_diag_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(avals_global, gks->src.red_integ_diag_global, sizeof(double[num_mom]));
  }

  if (gks->info.source.diagnostics.time_integrated) {
    // This assumes time-independent sources. For time dependent ones
    // step the source contributions in RK3 like we do for boundary fluxes.
    double avals_global_prev[num_mom];
    for (int k=0; k<num_mom; k++)
      avals_global_prev[k] = 0.0;
    gkyl_dynvec_getlast(gks->src.integ_diag, avals_global_prev);
  
    double tau = tm - tm_prev;
    for (int k=0; k<num_mom; k++)
      avals_global[k] = avals_global_prev[k] + tau*avals_global[k];
  }

  gkyl_dynvec_append(gks->src.integ_diag, tm, avals_global);

  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += 1;
}

void
gk_species_source_write_integrated_mom_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
}

void
gk_species_source_write_integrated_mom_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Write out integrated diagnostic moments.
    const char *fmt = "%s-%s_source_%s.gkyl";

    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "integrated_moms");
    char fileNm[sz+1]; // Ensures no buffer overflow.
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "integrated_moms");

    if (gks->src.is_first_integ_write_call) {
      gkyl_dynvec_write(gks->src.integ_diag, fileNm);
      gks->src.is_first_integ_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(gks->src.integ_diag, fileNm);
    }
  }
  gkyl_dynvec_clear(gks->src.integ_diag);

  app->stat.diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag_io += 1;
}

void 
gk_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src)
{
  src->source_id = s->info.source.source_id;

  // Default function pointers.
  src->write_func = gk_species_source_write_disabled;
  src->write_mom_func = gk_species_source_write_mom_disabled;
  src->calc_integrated_mom_func = gk_species_source_calc_integrated_mom_disabled;
  src->write_integrated_mom_func = gk_species_source_write_integrated_mom_disabled;

  if (src->source_id) {
    // Allocate source array.
    src->source = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
    src->source_host = src->source;
    if (app->use_gpu) {
      src->source_host = mkarr(false, src->source->ncomp, src->source->size); 
    }

    src->evolve = s->info.source.evolve; // Whether the source is time dependent.

    src->num_sources = s->info.source.num_sources;
    for (int k=0; k<s->info.source.num_sources; k++) {
      gk_species_projection_init(app, s, s->info.source.projection[k], &src->proj_source[k]);
    }

    // Allocate data and updaters for diagnostic moments.
    src->num_diag_mom = s->info.source.diagnostics.num_diag_moments;
    if (src->num_diag_mom == 0) {
      src->num_diag_mom = s->info.num_diag_moments;
      for (int m=0; m<src->num_diag_mom; ++m) {
        strcpy(s->info.source.diagnostics.diag_moments[m], s->info.diag_moments[m]);
      }
    }

    src->moms = gkyl_malloc(sizeof(struct gk_species_moment[src->num_diag_mom]));
    for (int m=0; m<src->num_diag_mom; ++m) {
      gk_species_moment_init(app, s, &src->moms[m], s->info.source.diagnostics.diag_moments[m], false);
    }

    // Allocate data and updaters for integrated moments.
    src->num_diag_int_mom = s->info.source.diagnostics.num_integrated_diag_moments;
    assert(src->num_diag_int_mom < 2); // 1 int moment allowed now.
    if (src->evolve || src->num_diag_int_mom > 0) {
      gk_species_moment_init(app, s, &src->integ_moms,
        src->num_diag_int_mom == 0? "FourMoments" : s->info.source.diagnostics.integrated_diag_moments[0], true);
      int num_mom = src->integ_moms.num_mom;
      if (app->use_gpu) {
        src->red_integ_diag = gkyl_cu_malloc(sizeof(double[num_mom]));
        src->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[num_mom]));
      } 
      else {
        src->red_integ_diag = gkyl_malloc(sizeof(double[num_mom]));
        src->red_integ_diag_global = gkyl_malloc(sizeof(double[num_mom]));
      }
      // Allocate dynamic-vector to store all-reduced integrated moments.
      src->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, num_mom);
      src->is_first_integ_write_call = true;
    }
    
    // Set function pointers chosen at runtime.
    if (src->evolve) {
      src->write_func = gk_species_source_write_enabled;
      src->write_mom_func = gk_species_source_write_mom_enabled;
      src->calc_integrated_mom_func = gk_species_source_calc_integrated_mom_enabled;
      src->write_integrated_mom_func = gk_species_source_write_integrated_mom_enabled;
    }
    else {
      src->write_func = gk_species_source_write_init_only;
      src->write_mom_func = gk_species_source_write_mom_init_only;
      if (src->num_diag_int_mom > 0) {
        // User requested integrated diagnostics.
        src->calc_integrated_mom_func = gk_species_source_calc_integrated_mom_enabled;
        src->write_integrated_mom_func = gk_species_source_write_integrated_mom_enabled;
      } else {
        src->calc_integrated_mom_func = gk_species_source_calc_integrated_mom_disabled;
        src->write_integrated_mom_func = gk_species_source_write_integrated_mom_disabled;
      }
    }
  }
}

void
gk_species_source_calc(gkyl_gyrokinetic_app *app, const struct gk_species *s, 
  struct gk_source *src, struct gkyl_array *f_buffer, double tm)
{
  if (src->source_id) {
    for (int k=0; k<s->info.source.num_sources; k++) {
      gk_species_projection_calc(app, s, &src->proj_source[k], f_buffer, tm);
      gkyl_array_accumulate(src->source, 1., f_buffer);
    }
  }
}

void
gk_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  if (src->source_id) {
    gkyl_array_accumulate(rhs, 1.0, src->source);
  }
}

void
gk_species_source_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gks->src.write_func(app, gks, tm, frame);
}

void
gk_species_source_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gks->src.write_mom_func(app, gks, tm, frame);
}

void
gk_species_source_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  gks->src.calc_integrated_mom_func(app, gks, tm);
}

void
gk_species_source_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  gks->src.write_integrated_mom_func(app, gks);
}

void
gk_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src)
{
  if (src->source_id) {
    gkyl_array_release(src->source);
    if (app->use_gpu) {
      gkyl_array_release(src->source_host);
    }
    for (int k=0; k<src->num_sources; k++) {
      gk_species_projection_release(app, &src->proj_source[k]);
    }

    // Release moment data.
    for (int i=0; i<src->num_diag_mom; ++i) {
      gk_species_moment_release(app, &src->moms[i]);
    }
    gkyl_free(src->moms);

    if (src->evolve || src->num_diag_int_mom > 0) {
      gk_species_moment_release(app, &src->integ_moms); 
      if (app->use_gpu) {
        gkyl_cu_free(src->red_integ_diag);
        gkyl_cu_free(src->red_integ_diag_global);
      }
      else {
        gkyl_free(src->red_integ_diag);
        gkyl_free(src->red_integ_diag_global);
      }  
      gkyl_dynvec_release(src->integ_diag);
    }
  }
}
