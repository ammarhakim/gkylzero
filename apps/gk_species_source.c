#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src)
{
  src->source_id = s->info.source.source_id;

  if (src->source_id) {
    int vdim = app->vdim;  
    src->calc_bflux = false;
    if (src->source_id == GKYL_BFLUX_SOURCE) {
      src->calc_bflux = true;
      assert(s->info.source.source_length);
      assert(s->info.source.source_species);
      src->source_length = s->info.source.source_length;
      src->source_species = gk_find_species(app, s->info.source.source_species);
      src->source_species_idx = gk_find_species_idx(app, s->info.source.source_species);
      if (app->use_gpu) {
        src->scale_ptr = gkyl_cu_malloc((vdim+2)*sizeof(double));
      }
      else {
        src->scale_ptr = gkyl_malloc((vdim+2)*sizeof(double));
      }
    }

    // Allocate source array.
    src->source = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
    src->source_host = src->source;
    if (app->use_gpu) {
      src->source_host = mkarr(false, src->source->ncomp, src->source->size); 
    }
    src->scale_factor = 1.0;

    src->evolve = s->info.source.evolve; // Whether the source is time dependent.

    src->num_sources = s->info.source.num_sources;
    for (int k=0; k<s->info.source.num_sources; k++) {
      gk_species_projection_init(app, s, s->info.source.projection[k], &src->proj_source[k]);
    }

    // Allocate data and updaters for diagnostic moments.
    src->num_diag_moments = s->info.source.diagnostics.num_diag_moments;
    if (src->num_diag_moments == 0) {
      src->num_diag_moments = s->info.num_diag_moments;
      for (int m=0; m<src->num_diag_moments; ++m) {
        strcpy(s->info.source.diagnostics.diag_moments[m], s->info.diag_moments[m]);
      }
    }

    s->src.moms = gkyl_malloc(sizeof(struct gk_species_moment[src->num_diag_moments]));
    for (int m=0; m<src->num_diag_moments; ++m) {
      gk_species_moment_init(app, s, &s->src.moms[m], s->info.source.diagnostics.diag_moments[m], false);
    }

    // Allocate data and updaters for integrated moments.
    int src_num_diag_int_moms = s->info.source.diagnostics.num_integrated_diag_moments;
    assert(src_num_diag_int_moms < 2); // 1 int moment allowed now.
    gk_species_moment_init(app, s, &s->src.integ_moms,
      src_num_diag_int_moms == 0? "FourMoments" : s->info.source.diagnostics.integrated_diag_moments[0], true);
    int num_mom = s->src.integ_moms.num_mom;
    if (app->use_gpu) {
      s->src.red_integ_diag = gkyl_cu_malloc(sizeof(double[num_mom]));
      s->src.red_integ_diag_global = gkyl_cu_malloc(sizeof(double[num_mom]));
    } 
    else {
      s->src.red_integ_diag = gkyl_malloc(sizeof(double[num_mom]));
      s->src.red_integ_diag_global = gkyl_malloc(sizeof(double[num_mom]));
    }
    // Allocate dynamic-vector to store all-reduced integrated moments.
    s->src.integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, num_mom);
    s->src.is_first_integ_write_call = true;
  }
}

void
gk_species_source_calc(gkyl_gyrokinetic_app *app, const struct gk_species *s, 
  struct gk_source *src, double tm)
{
  if (src->source_id) {
    struct gkyl_array *source_tmp = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
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
  struct gk_source *src, const struct gkyl_array *fin[], struct gkyl_array *rhs[])
{
  if (src->source_id) {
    int species_idx;
    species_idx = gk_find_species_idx(app, s->info.name);
    // use boundary fluxes to scale source profile
    if (src->calc_bflux) {
      src->scale_factor = 0.0;
      double z[app->basis.num_basis];
      double red_mom[1] = { 0.0 };

      for (int d=0; d<app->cdim; ++d) {
        gkyl_array_reduce(src->scale_ptr, src->source_species->bflux_solver.gammai[2*d].marr, GKYL_SUM);
        if (app->use_gpu) {
          gkyl_cu_memcpy(red_mom, src->scale_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
        }
        else {
          red_mom[0] = src->scale_ptr[0];
        }
        double red_mom_global[1] = { 0.0 };
        gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, red_mom, red_mom_global);
        src->scale_factor += red_mom_global[0];
        gkyl_array_reduce(src->scale_ptr, src->source_species->bflux_solver.gammai[2*d+1].marr, GKYL_SUM);
        if (app->use_gpu) {
          gkyl_cu_memcpy(red_mom, src->scale_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
        }
        else {
          red_mom[0] = src->scale_ptr[0];
        }
        gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, red_mom, red_mom_global);
        src->scale_factor += red_mom_global[0];
      }
      src->scale_factor = src->scale_factor/src->source_length;
    }
    gkyl_array_accumulate(rhs[species_idx], src->scale_factor, src->source);
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
        .basis_type = app->basis.id
      }
    );

    for (int m=0; m<gks->src.num_diag_moments; ++m) {
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
}

void
gk_species_source_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  if (gks->src.source_id && gks->src.evolve) {
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
}

void
gk_species_source_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  if (gks->src.source_id && gks->src.evolve) {
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
}

void
gk_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src)
{
  if (src->source_id) {
    gkyl_array_release(src->source);
    if (app->use_gpu) {
      gkyl_array_release(src->source_host);
    }

    if (src->calc_bflux) {
      if (app->use_gpu) {
        gkyl_cu_free(src->scale_ptr);
      } 
      else {
        gkyl_free(src->scale_ptr);
      }
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
