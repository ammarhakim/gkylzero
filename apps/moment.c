#include <gkyl_moment_priv.h>

gkyl_moment_app*
gkyl_moment_app_new(struct gkyl_moment *mom)
{
  disable_denorm_float();
  
  struct gkyl_moment_app *app = gkyl_malloc(sizeof(gkyl_moment_app));

  int ndim = app->ndim = mom->ndim;
  strcpy(app->name, mom->name);
  app->tcurr = 0.0; // reset on init

  int ghost[3] = { 2, 2, 2 }; // 2 ghost-cells for wave
  if (mom->scheme_type == GKYL_MOMENT_MP)
    for (int d=0; d<3; ++d) ghost[d] = 3; // 3 for MP scheme
  
  gkyl_rect_grid_init(&app->grid, ndim, mom->lower, mom->upper, mom->cells);
  gkyl_create_grid_ranges(&app->grid, ghost, &app->local_ext, &app->local);
  
  skin_ghost_ranges_init(&app->skin_ghost, &app->local_ext, ghost);

  app->c2p_ctx = app->mapc2p = 0;  
  app->has_mapc2p = mom->mapc2p ? true : false;

  if (app->has_mapc2p) {
    // initialize computational to physical space mapping
    app->c2p_ctx = mom->c2p_ctx;
    app->mapc2p = mom->mapc2p;

    // we project mapc2p on p=1 basis functions
    struct gkyl_basis basis;
    gkyl_cart_modal_tensor(&basis, ndim, 1);

    // initialize DG field representing mapping
    struct gkyl_array *c2p = mkarr(false, ndim*basis.num_basis, app->local_ext.volume);
    gkyl_eval_on_nodes *ev_c2p = gkyl_eval_on_nodes_new(&app->grid, &basis, ndim, mom->mapc2p, mom->c2p_ctx);
    gkyl_eval_on_nodes_advance(ev_c2p, 0.0, &app->local_ext, c2p);

    // write DG projection of mapc2p to file
    const char *fmt = "%s-mapc2p.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name);
    gkyl_grid_sub_array_write(&app->grid, &app->local, c2p, fileNm);

    gkyl_array_release(c2p);
    gkyl_eval_on_nodes_release(ev_c2p);
  }

  // create geometry object
  app->geom = gkyl_wave_geom_new(&app->grid, &app->local_ext,
    app->mapc2p, app->c2p_ctx);

  double cfl_frac = mom->cfl_frac == 0 ? 0.95 : mom->cfl_frac;
  app->cfl = 1.0*cfl_frac;

  app->num_periodic_dir = mom->num_periodic_dir;
  for (int d=0; d<ndim; ++d)
    app->periodic_dirs[d] = mom->periodic_dirs[d];

  // construct list of directions to skip
  for (int d=0; d<3; ++d)
    app->is_dir_skipped[d] = 0;
  for (int i=0; i<mom->num_skip_dirs; ++i)
    app->is_dir_skipped[mom->skip_dirs[i]] = 1;

  app->has_field = 0;
  // initialize field if we have one
  if (mom->field.init) {
    app->has_field = 1;
    moment_field_init(mom, &mom->field, app, &app->field);
  }

  int ns = app->num_species = mom->num_species;
  // allocate space to store species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct moment_species[ns])) : 0;
  // create species grid & ranges
  for (int i=0; i<ns; ++i)
    moment_species_init(mom, &mom->species[i], app, &app->species[i]);

  // check if we should update sources
  app->update_sources = 0;
  if (app->has_field && ns>0) {
    app->update_sources = 1; // only update if field and species are present
    moment_coupling_init(app, &app->sources);
  }

  // initialize stat object to all zeros
  app->stat = (struct gkyl_moment_stat) {
  };

  return app;
}

double
gkyl_moment_app_max_dt(gkyl_moment_app* app)
{
  double max_dt = DBL_MAX;
  for (int i=0;  i<app->num_species; ++i) 
    max_dt = fmin(max_dt, moment_species_max_dt(app, &app->species[i]));

  if (app->has_field)
    max_dt = fmin(max_dt, moment_field_max_dt(app, &app->field));

  return max_dt;
}

void
gkyl_moment_app_apply_ic(gkyl_moment_app* app, double t0)
{
  app->tcurr = t0;
  gkyl_moment_app_apply_ic_field(app, t0);
  for (int i=0;  i<app->num_species; ++i)
    gkyl_moment_app_apply_ic_species(app, i, t0);
}

void
gkyl_moment_app_apply_ic_field(gkyl_moment_app* app, double t0)
{
  if (app->has_field != 1) return;

  app->tcurr = t0;
  gkyl_fv_proj *proj = gkyl_fv_proj_new(&app->grid, 2, 8, app->field.init, app->field.ctx);
  
  gkyl_fv_proj_advance(proj, t0, &app->local, app->field.f[0]);
  gkyl_fv_proj_release(proj);

  moment_field_apply_bc(app, t0, &app->field, app->field.f[0]);
}

void
gkyl_moment_app_apply_ic_species(gkyl_moment_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  gkyl_fv_proj *proj = gkyl_fv_proj_new(&app->grid, 2, app->species[sidx].num_equations,
    app->species[sidx].init, app->species[sidx].ctx);
  
  gkyl_fv_proj_advance(proj, t0, &app->local, app->species[sidx].f[0]);
  gkyl_fv_proj_release(proj);

  moment_species_apply_bc(app, t0, &app->species[sidx], app->species[sidx].f[0]);
}

void
gkyl_moment_app_write(const gkyl_moment_app* app, double tm, int frame)
{
  gkyl_moment_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i)
    gkyl_moment_app_write_species(app, i, tm, frame);
}

void
gkyl_moment_app_write_field(const gkyl_moment_app* app, double tm, int frame)
{
  if (app->has_field != 1) return;

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "field", frame);
  gkyl_grid_sub_array_write(&app->grid, &app->local, app->field.f[0], fileNm.str);
  cstr_drop(&fileNm);

  // write external EM field if it is present
  if (app->field.ext_em) {
    cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "ext_em_field", frame);
    gkyl_grid_sub_array_write(&app->grid, &app->local, app->field.ext_em, fileNm.str);
    cstr_drop(&fileNm);
  }
}

void
gkyl_moment_app_write_field_energy(gkyl_moment_app *app)
{
  if (app->has_field) {
    // write out field energy
    cstr fileNm = cstr_from_fmt("%s-field-energy.gkyl", app->name);

    if (app->field.is_first_energy_write_call) {
      // write to a new file (this ensure previous output is removed)
      gkyl_dynvec_write(app->field.integ_energy, fileNm.str);
      app->field.is_first_energy_write_call = false;
    }
    else {
      // append to existing file
      gkyl_dynvec_awrite(app->field.integ_energy, fileNm.str);
    }
    gkyl_dynvec_clear(app->field.integ_energy);
    
    cstr_drop(&fileNm);
  }
}

void
gkyl_moment_app_write_integrated_mom(gkyl_moment_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    // write out diagnostic moments
    cstr fileNm = cstr_from_fmt("%s-%s-%s.gkyl", app->name, app->species[i].name,
      "imom");
    
    if (app->species[i].is_first_q_write_call) {
      gkyl_dynvec_write(app->species[i].integ_q, fileNm.str);
      app->species[i].is_first_q_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(app->species[i].integ_q, fileNm.str);
    }
    gkyl_dynvec_clear(app->species[i].integ_q);

    cstr_drop(&fileNm);
  }
}

void
gkyl_moment_app_write_species(const gkyl_moment_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, app->species[sidx].name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[sidx].name, frame);
  
  gkyl_grid_sub_array_write(&app->grid, &app->local, app->species[sidx].f[0], fileNm);
}

// internal function that takes a single time-step
static struct gkyl_update_status
moment_update(gkyl_moment_app* app, double dt0)
{
  int ns = app->num_species, ndim = app->ndim;
  bool have_nans_occured = false;
  
  double dt_suggested = DBL_MAX;
  
  // time-stepper states
  enum {
    UPDATE_DONE = 0,
    PRE_UPDATE,
    POST_UPDATE,
    FIRST_COUPLING_UPDATE,
    FIELD_UPDATE,
    SPECIES_UPDATE,
    SECOND_COUPLING_UPDATE,
    UPDATE_REDO,
  } state = PRE_UPDATE;

  double tcurr = app->tcurr, dt = dt0;
  while (state != UPDATE_DONE) {
    switch (state) {
      case PRE_UPDATE:
        state = FIRST_COUPLING_UPDATE; // next state
          
        // copy old solution in case we need to redo this step
        for (int i=0; i<ns; ++i)
          gkyl_array_copy(app->species[i].fdup, app->species[i].f[0]);
        if (app->has_field)
          gkyl_array_copy(app->field.fdup, app->field.f[0]);

        break;
          
      
      case FIRST_COUPLING_UPDATE:
        state = FIELD_UPDATE; // next state

        if (app->update_sources) {
          struct timespec src1_tm = gkyl_wall_clock();
          moment_coupling_update(app, &app->sources, 0, tcurr, dt/2);
          app->stat.sources_tm += gkyl_time_diff_now_sec(src1_tm);
        }
            
        break;

      case FIELD_UPDATE:
        state = SPECIES_UPDATE; // next state

        if (app->has_field) {
          struct timespec fl_tm = gkyl_wall_clock();
          struct gkyl_update_status s = moment_field_update(app, &app->field, tcurr, dt);
          if (!s.success) {
            app->stat.nfail += 1;
            dt = s.dt_suggested;
            state = UPDATE_REDO;
            break;
          }
            
          dt_suggested = fmin(dt_suggested, s.dt_suggested);
          app->stat.field_tm += gkyl_time_diff_now_sec(fl_tm);
        }
          
        break;

      case SPECIES_UPDATE:
        state = SECOND_COUPLING_UPDATE; // next state

        struct timespec sp_tm = gkyl_wall_clock();
        for (int i=0; i<ns; ++i) {         
          struct gkyl_update_status s =
            moment_species_update(app, &app->species[i], tcurr, dt);

          if (!s.success) {
            app->stat.nfail += 1;
            dt = s.dt_suggested;
            state = UPDATE_REDO;
            break;
          }
          dt_suggested = fmin(dt_suggested, s.dt_suggested);
        }
        app->stat.species_tm += gkyl_time_diff_now_sec(sp_tm);
         
        break;

      case SECOND_COUPLING_UPDATE:
        state = POST_UPDATE; // next state

        if (app->update_sources) {
          struct timespec src2_tm = gkyl_wall_clock();
          moment_coupling_update(app, &app->sources, 1, tcurr, dt/2);
          app->stat.sources_tm += gkyl_time_diff_now_sec(src2_tm);
        }

        break;

      case POST_UPDATE:
        state = UPDATE_DONE;

        // copy solution in prep for next time-step
        for (int i=0; i<ns; ++i) {
          // check for nans before copying
          if (check_for_nans(app->species[i].f[ndim], app->local))
            have_nans_occured = true;
          else // only copy in case no nans, so old solution can be written out
            gkyl_array_copy(app->species[i].f[0], app->species[i].f[ndim]);
        }
        
        if (app->has_field)
          gkyl_array_copy(app->field.f[0], app->field.f[ndim]);
          
        break;

      case UPDATE_REDO:
        state = PRE_UPDATE; // start all-over again
          
        // restore solution and retake step
        for (int i=0; i<ns; ++i)
          gkyl_array_copy(app->species[i].f[0], app->species[i].fdup);
        if (app->has_field)
          gkyl_array_copy(app->field.f[0], app->field.fdup);
          
        break;

      case UPDATE_DONE: // unreachable code! (suppresses warning)
        break;
    }
  }

  return (struct gkyl_update_status) {
    .success = have_nans_occured ? false : true,
    .dt_actual = dt,
    .dt_suggested = dt_suggested,
  };
}

struct gkyl_update_status
gkyl_moment_update(gkyl_moment_app* app, double dt)
{
  app->stat.nup += 1;
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_update_status status = moment_update(app, dt);
  app->tcurr += status.dt_actual;
  
  app->stat.total_tm += gkyl_time_diff_now_sec(wst);
  
  return status;
}

void
gkyl_moment_app_calc_field_energy(gkyl_moment_app* app, double tm)
{
  if (app->has_field) {
    double energy[6] = { 0.0 };
    calc_integ_quant(6, app->grid.cellVolume, app->field.f[0], app->geom,
      app->local, integ_sq, energy);
    gkyl_dynvec_append(app->field.integ_energy, tm, energy);
  }
}

void
gkyl_moment_app_calc_integrated_mom(gkyl_moment_app *app, double tm)
{
  for (int sidx=0; sidx<app->num_species; ++sidx) {
    int meqn = app->species[sidx].num_equations;
    double q_integ[meqn];
    calc_integ_quant(meqn, app->grid.cellVolume, app->species[sidx].f[0], app->geom,
      app->local, integ_unit, q_integ);
    gkyl_dynvec_append(app->species[sidx].integ_q, tm, q_integ);
  }
}

struct gkyl_moment_stat
gkyl_moment_app_stat(gkyl_moment_app* app)
{
  return app->stat;
}

void
gkyl_moment_app_stat_write(const gkyl_moment_app* app)
{
  const char *fmt = "%s-%s";
  int sz = gkyl_calc_strlen(fmt, app->name, "stat.json");
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, "stat.json");

  char buff[70];
  time_t t = time(NULL);
  struct tm curr_tm = *localtime(&t);

  // compute total number of cells updated in simulation
  long tot_cells_up = app->local.volume*app->num_species*app->ndim*app->stat.nup;

  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  with_file (fp, fileNm, "a") {
    fprintf(fp, "{\n");

    if (strftime(buff, sizeof buff, "%c", &curr_tm))
      fprintf(fp, " date : %s\n", buff);

    fprintf(fp, " nup : %ld,\n", app->stat.nup);
    fprintf(fp, " nfail : %ld,\n", app->stat.nfail);
    fprintf(fp, " total_tm : %lg,\n", app->stat.total_tm);
    fprintf(fp, " species_tm : %lg,\n", app->stat.species_tm);
    fprintf(fp, " field_tm : %lg,\n", app->stat.field_tm);
    fprintf(fp, " sources_tm : %lg\n", app->stat.sources_tm);

    for (int i=0; i<app->num_species; ++i) {
      long tot_bad_cells = 0L;      
      for (int d=0; d<app->ndim; ++d) {
        struct gkyl_wave_prop_stats wvs = gkyl_wave_prop_stats(app->species[i].slvr[d]);
        fprintf(fp, " %s_n_bad_1D_sweeps[%d] = %ld\n", app->species[i].name, d, wvs.n_bad_advance_calls);
        fprintf(fp, " %s_n_bad_cells[%d] = %ld\n", app->species[i].name, d, wvs.n_bad_cells);
        fprintf(fp, " %s_n_max_bad_cells[%d] = %ld\n", app->species[i].name, d, wvs.n_max_bad_cells);

        tot_bad_cells += wvs.n_bad_cells;
      }
      fprintf(fp, " %s_bad_cell_frac = %lg\n", app->species[i].name, (double) tot_bad_cells/tot_cells_up );
    }

  
    fprintf(fp, "}\n");
  }
}

void
gkyl_moment_app_release(gkyl_moment_app* app)
{
  for (int i=0; i<app->num_species; ++i)
    moment_species_release(&app->species[i]);
  gkyl_free(app->species);

  if (app->has_field)
    moment_field_release(&app->field);

  if (app->update_sources)
    moment_coupling_release(&app->sources);

  gkyl_wave_geom_release(app->geom);

  gkyl_free(app);
}
