#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src)
{
  src->source_id = s->info.source.source_id;

  if (src->source_id) {
    int vdim = app->vdim;
    // we need to ensure source has same shape as distribution function
    src->source = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
    src->source_host = src->source;
    if (app->use_gpu) {
      src->source_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);
    }

    src->evolve = s->info.source.evolve; // Whether the source is time dependent.

    src->num_sources = s->info.source.num_sources;
    for (int k=0; k<s->info.source.num_sources; k++) {
      gk_species_projection_init(app, s, s->info.source.projection[k], &src->proj_source[k]);
    }

    // Allocate data and updaters for diagnostic moments.
    src->num_diag_moments = s->info.num_diag_moments;
    s->src.moms = gkyl_malloc(sizeof(struct gk_species_moment[src->num_diag_moments]));
    for (int m=0; m<src->num_diag_moments; ++m) {
      gk_species_moment_init(app, s, &s->src.moms[m], s->info.diag_moments[m]);
    }

    // Allocate data and updaters for integrated moments.
    gk_species_moment_init(app, s, &s->src.integ_moms, "Integrated");
    if (app->use_gpu) {
      s->src.red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
      s->src.red_integ_diag_global = gkyl_cu_malloc(sizeof(double[vdim+2]));
    } 
    else {
      s->src.red_integ_diag = gkyl_malloc(sizeof(double[vdim+2]));
      s->src.red_integ_diag_global = gkyl_malloc(sizeof(double[vdim+2]));
    }
    // allocate dynamic-vector to store all-reduced integrated moments 
    s->src.integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
    s->src.is_first_integ_write_call = true;
  }
}

void
gk_species_source_calc(gkyl_gyrokinetic_app *app, const struct gk_species *s, 
  struct gk_source *src, double tm)
{
  if (src->source_id) {
    struct gkyl_array *source_tmp = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
    for (int k=0; k<s->info.source.num_sources; k++) {
      gk_species_projection_calc(app, s, &src->proj_source[k], source_tmp, tm);
      gkyl_array_accumulate(src->source, 1., source_tmp);
    }
    gkyl_array_release(source_tmp);
  }
}

// Compute rhs of the source.
void
gk_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  if (src->source_id) {
    gkyl_array_accumulate(rhs, 1.0, src->source);
  }
}

// Source write funcs
void
gk_species_source_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  if (gks->src.source_id && (gks->src.evolve || frame == 0)) {
    struct timespec wst = gkyl_wall_clock();
    struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->basis.id
      }
    );

    // Write out the source distribution function
    const char *fmt = "%s-%s_source_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

    // copy data from device to host before writing it out
    if (app->use_gpu) {
      gkyl_array_copy(gks->src.source_host, gks->src.source);
    }

    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->src.source_host, fileNm);
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;

    gk_array_meta_release(mt);   
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gk_species_source_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  if (gks->src.source_id && (gks->src.evolve || frame == 0)) {
    struct timespec wst = gkyl_wall_clock();
    struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->confBasis.id
      }
    );

    for (int m=0; m<gks->info.num_diag_moments; ++m) {
      gk_species_moment_calc(&gks->src.moms[m], gks->local, app->local, gks->src.source);
      app->stat.nmom += 1;

      const char *fmt = "%s-%s_source_%s_%d.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name,
        gks->info.diag_moments[m], frame);
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name,
        gks->info.diag_moments[m], frame);

      if (!gks->src.moms[m].is_bimaxwellian_moms && !gks->src.moms[m].is_maxwellian_moms) {
        // Rescale moment by inverse of Jacobian
        gkyl_dg_div_op_range(gks->moms[m].mem_geo, app->confBasis, 
          0, gks->src.moms[m].marr, 0, gks->src.moms[m].marr, 0, 
          app->gk_geom->jacobgeo, &app->local);      
      }

      if (app->use_gpu) {
        gkyl_array_copy(gks->src.moms[m].marr_host, gks->src.moms[m].marr);
      }

      struct timespec wtm = gkyl_wall_clock();
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
        gks->src.moms[m].marr_host, fileNm);
      app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.nio += 1;
    }

    gk_array_meta_release(mt);   
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gk_species_source_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  if (gks->src.source_id && gks->src.evolve) {
    struct timespec wst = gkyl_wall_clock();

    int vdim = app->vdim;
    double avals_global[2+vdim];

    gk_species_moment_calc(&gks->src.integ_moms, gks->local, app->local, gks->src.source); 
    // reduce to compute sum over whole domain, append to diagnostics
    gkyl_array_reduce_range(gks->src.red_integ_diag, gks->src.integ_moms.marr, GKYL_SUM, &app->local);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, 
      gks->src.red_integ_diag, gks->src.red_integ_diag_global);
    if (app->use_gpu) {
      gkyl_cu_memcpy(avals_global, gks->src.red_integ_diag_global, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      memcpy(avals_global, gks->src.red_integ_diag_global, sizeof(double[2+vdim]));
    }
    gkyl_dynvec_append(gks->src.integ_diag, tm, avals_global);
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gk_species_source_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  if (gks->src.source_id && gks->src.evolve) {
    struct timespec wst = gkyl_wall_clock();

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out integrated diagnostic moments
      const char *fmt = "%s-%s_source_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "integrated_moms");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "integrated_moms");

      struct timespec wtm = gkyl_wall_clock();
      if (gks->src.is_first_integ_write_call) {
        gkyl_dynvec_write(gks->src.integ_diag, fileNm);
        gks->src.is_first_integ_write_call = false;
      }
      else {
        gkyl_dynvec_awrite(gks->src.integ_diag, fileNm);
      }
      app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.nio += 1;
    }
    gkyl_dynvec_clear(gks->src.integ_diag);
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
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
    for (int i=0; i<src->num_diag_moments; ++i) {
      gk_species_moment_release(app, &src->moms[i]);
    }
    gkyl_free(src->moms);
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
