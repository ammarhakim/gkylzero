#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src)
{
  src->source_id = s->info.source.source_id;

  if (src->source_id) {
    // Allocate source array.
    src->source = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
    src->source_host = src->source;
    if (app->use_gpu) {
      src->source_host = mkarr(false, src->source->ncomp, src->source->size); 
    }

    src->evolve = s->info.source.evolve || (s->info.source.num_adapt_sources > 0); // Whether the source is time dependent.

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

    // Set up the adaptive source.
    src->num_adapt_sources = s->info.source.num_adapt_sources;
    assert(src->num_adapt_sources <= src->num_sources); // Adaptive source should be a subset of the sources.
    if(src->num_adapt_sources > 0){
      for (int k = 0; k < src->num_adapt_sources; ++k) {
        struct gk_adapt_source *adapt_src = &src->adapt[k];

        adapt_src->adapt_particle = s->info.source.adapt[k].adapt_particle;
        adapt_src->adapt_energy = s->info.source.adapt[k].adapt_energy;

        adapt_src->adapt_species = gk_find_species(app, s->info.source.adapt[k].adapt_species_name);

        adapt_src->particle_src_curr = s->info.source.projection[k].particle;
        adapt_src->energy_src_curr = s->info.source.projection[k].energy;
        // The temperature computation makes sense only if we inject particles
        adapt_src->temperature_curr = s->info.source.projection[k].particle > 0?
          2./3. * adapt_src->energy_src_curr/adapt_src->particle_src_curr : 1.0;

        gk_species_moment_init(app, s, &adapt_src->integ_mom, "ThreeMoments", true);

        int num_mom = adapt_src->integ_mom.num_mom;
        if (app->use_gpu){
          adapt_src->red_integ_mom = gkyl_cu_malloc(sizeof(double[num_mom]));
          adapt_src->red_integ_mom_global = gkyl_cu_malloc(sizeof(double[num_mom]));
        }
        else {
          adapt_src->red_integ_mom = gkyl_malloc(sizeof(double[num_mom]));
          adapt_src->red_integ_mom_global = gkyl_malloc(sizeof(double[num_mom]));
        }

        adapt_src->num_boundaries = s->info.source.adapt[k].num_boundaries;
        adapt_src->range_bflux = gkyl_malloc(sizeof(struct gkyl_range) * adapt_src->num_boundaries);
        adapt_src->range_mom = gkyl_malloc(sizeof(struct gkyl_range) * adapt_src->num_boundaries);
        adapt_src->range_conf = gkyl_malloc(sizeof(struct gkyl_range) * adapt_src->num_boundaries);
        for (int j=0; j < adapt_src->num_boundaries; ++j) {
          int dir = s->info.source.adapt[k].dir[j];
          int edge = s->info.source.adapt[k].edge[j];
          adapt_src->range_bflux[j] = edge == GKYL_LOWER_EDGE ? s->lower_ghost[dir] : s->upper_ghost[dir];
          adapt_src->range_mom[j] = edge == GKYL_LOWER_EDGE ? s->lower_ghost[dir] : s->upper_ghost[dir];
          adapt_src->range_conf[j] = edge == GKYL_LOWER_EDGE ? app->lower_ghost[dir] : app->upper_ghost[dir];
          adapt_src->dir[j] = dir;
          adapt_src->edge[j] = edge;

          if (edge == GKYL_LOWER_EDGE) {
            // Specific scenario if we are in a inner wall limited case.
            if (s->lower_bc[dir].type == GKYL_SPECIES_GK_IWL) { 
              adapt_src->range_mom[j] = s->lower_ghost_par_sol;
              // need to create a configuration space ghost range that encompasses the SOL only.
              double xLCFS = s->lower_bc[j].aux_parameter;
              int idxLCFS_m = (xLCFS-1e-8 - app->grid.lower[0])/app->grid.dx[0]+1;
              gkyl_range_shorten_from_below(&adapt_src->range_conf[j], &app->lower_ghost[j], 0, app->grid.cells[0]-idxLCFS_m+1);
            }
            // We check that the boundaries are not set to zero flux (prevents the call of bflux).
            assert(s->lower_bc[dir].type != GKYL_SPECIES_ZERO_FLUX);
          }
          else {
            if (s->upper_bc[dir].type == GKYL_SPECIES_GK_IWL) { 
              adapt_src->range_mom[j] = s->upper_ghost_par_sol;
              double xLCFS = s->lower_bc[j].aux_parameter;
              int idxLCFS_m = (xLCFS-1e-8 - app->grid.lower[0])/app->grid.dx[0]+1;
              gkyl_range_shorten_from_below(&adapt_src->range_conf[j], &app->upper_ghost[j], 0, app->grid.cells[0]-idxLCFS_m+1);
            }
            assert(s->upper_bc[dir].type != GKYL_SPECIES_ZERO_FLUX);
          }
        }
      }
    }
  }
}

void
gk_species_source_calc(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src, double tm)
{
  if (src->source_id) {
    gkyl_array_clear(src->source, 0.0);
    struct gkyl_array *source_tmp = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
    for (int k=0; k<s->info.source.num_sources; k++) {
      gk_species_projection_calc(app, s, &src->proj_source[k], source_tmp, tm);
      gkyl_array_accumulate(src->source, 1., source_tmp);
    }
    gkyl_array_release(source_tmp);
  }
}

void
gk_species_source_adapt(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src, double tm) {  

  for (int k=0; k < s->info.source.num_adapt_sources; ++k) {
    struct gk_adapt_source *adapt_src = &src->adapt[k];
    struct gk_species *s_adapt = adapt_src->adapt_species;

    adapt_src->particle_rate_loss = 0.0;
    adapt_src->energy_rate_loss = 0.0;
    int num_mom = adapt_src->integ_mom.num_mom;
    // Accumulate energy and particle losses through the boundaries considered by the current adaptive source.
    for (int j=0; j < adapt_src->num_boundaries; ++j) {

      // Compute the moment of the bflux to get the integrated loss.
      gk_species_bflux_get_flux(&s_adapt->bflux, adapt_src->dir[j], adapt_src->edge[j], src->source, &adapt_src->range_bflux[j]);
      gk_species_moment_calc(&adapt_src->integ_mom, adapt_src->range_mom[j], adapt_src->range_conf[j], src->source);
      app->stat.n_mom += 1;
      
      // Reduce the moment over the specified range and store it in the global array.
      double red_int_mom_global[num_mom];
      gkyl_array_reduce_range(adapt_src->red_integ_mom, adapt_src->integ_mom.marr, GKYL_SUM, &adapt_src->range_conf[j]);
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom, adapt_src->red_integ_mom, adapt_src->red_integ_mom_global);
      if (app->use_gpu) {
        gkyl_cu_memcpy(red_int_mom_global, adapt_src->red_integ_mom_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
      }
      else {
        memcpy(red_int_mom_global, adapt_src->red_integ_mom_global, sizeof(double[num_mom]));
      }
      adapt_src->particle_rate_loss += red_int_mom_global[0]; // n
      adapt_src->energy_rate_loss += 0.5 * s_adapt->info.mass * red_int_mom_global[num_mom-1]; // 1/2 * m * v^2
    }

    // Particle and energy rate update.
    // balance = user target + loss
    double particle_src_new = adapt_src->adapt_particle? 
      s->info.source.projection[k].particle + adapt_src->particle_rate_loss
      : s->info.source.projection[k].particle;

    double energy_src_new = adapt_src->adapt_energy?
      s->info.source.projection[k].energy + adapt_src->energy_rate_loss
      : s->info.source.projection[k].energy;

    // Compute the target temperature of the source following the rule:
    // T = 2/3 * Q/G (T: src temperature, Q: src energy rate, G: total particle rate)
    double temperature_new = particle_src_new == 0.0 ? s->info.source.projection[k].temp_max/2 : 2./3. * energy_src_new/particle_src_new;

    // Update the density and temperature moments of the source
    gkyl_array_set(src->proj_source[k].dens, particle_src_new, src->proj_source[k].shape_conf);
    gkyl_array_set(src->proj_source[k].vtsq, temperature_new, src->proj_source[k].one_conf);
    // Upar=0 at initialization and is never changed.

    // Refresh the current values of particle, energy and temperature (can be used for control).
    adapt_src->particle_src_curr = particle_src_new;
    adapt_src->energy_src_curr = energy_src_new;
    adapt_src->temperature_curr = temperature_new;
  }

  // Reproject the source
  gk_species_source_calc(app, s, src, tm);
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

// Source write funcs.
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

    // Write out the source distribution function.
    const char *fmt = "%s-%s_source_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
    char fileNm[sz+1]; // Ensures no buffer overflow.
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
    for (int k = 0; k < src->num_adapt_sources; ++k) {
      gk_species_moment_release(app, &src->adapt[k].integ_mom);
      if (app->use_gpu){
        gkyl_cu_free(src->adapt[k].red_integ_mom);
        gkyl_cu_free(src->adapt[k].red_integ_mom_global);
        gkyl_cu_free(src->adapt[k].range_bflux);
        gkyl_cu_free(src->adapt[k].range_mom);
        gkyl_cu_free(src->adapt[k].range_conf);
      }
      else {
        gkyl_free(src->adapt[k].red_integ_mom);
        gkyl_free(src->adapt[k].red_integ_mom_global);
        gkyl_free(src->adapt[k].range_bflux);
        gkyl_free(src->adapt[k].range_mom);
        gkyl_free(src->adapt[k].range_conf);
      }
    }
  }
}
