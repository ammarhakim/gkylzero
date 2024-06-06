#include <gkyl_alloc.h>

#include <gkyl_gyrokinetic_mb.h>
#include <gkyl_gyrokinetic_mb_priv.h>
#include <gkyl_dflt.h>

gkyl_gyrokinetic_mb_app*
gkyl_gyrokinetic_mb_app_new(struct gkyl_gk_mb *inp)
{
  disable_denorm_float();

  gkyl_gyrokinetic_mb_app *app = gkyl_malloc(sizeof(gkyl_gyrokinetic_mb_app));

  strcpy(app->name, inp->name);

  app->cdim = inp->cdim;
  app->vdim = inp->vdim;
  app->poly_order = inp->poly_order;

  app->cfl = inp->cfl_frac == 0 ? 1.0 : inp->cfl_frac;

#ifdef GKYL_HAVE_CUDA
  app->use_gpu = inp->use_gpu;
#else
  app->use_gpu = false; // can't use GPUs if we don't have them!
#endif

  app->num_blocks = inp->num_blocks;
  app->num_blocks_local = inp->num_blocks;

  app->blocks = gkyl_malloc(app->num_blocks * sizeof(struct gkyl_gyrokinetic_app));

  app->btopo = gkyl_block_topo_new(app->cdim, app->num_blocks);

  for (int bidx=0; bidx<app->num_blocks; bidx++) {
    struct gkyl_gk *blinp = inp->blocks[bidx];

    // Block name = <the name of the app>_b#.
    const char *fmt = "%s_b%d";
    int sz = gkyl_calc_strlen(fmt, inp->name, bidx);
    char appNm[sz+1]; // ensures no buffer overflow
    snprintf(appNm, sizeof appNm, fmt, inp->name, bidx);
    strcpy(blinp->name, appNm);

    // Fill the input for each block with common values.
    blinp->cdim = inp->cdim;
    blinp->vdim = inp->vdim;
    blinp->poly_order = inp->poly_order;
    blinp->basis_type = inp->basis_type;

    blinp->use_gpu = inp->use_gpu;

    blinp->num_species = inp->num_species;
    for (int s=0; s<inp->num_species; s++)
      blinp->species[s] = inp->species[s];

    blinp->num_neut_species = inp->num_neut_species;
    for (int s=0; s<inp->num_neut_species; s++)
      blinp->neut_species[s] = inp->neut_species[s];

    blinp->skip_field = inp->skip_field;
    blinp->field = inp->field;

    blinp->num_periodic_dir = inp->num_periodic_dir;
    memcpy(blinp->periodic_dirs, inp->periodic_dirs, inp->num_periodic_dir*sizeof(int));

    app->btopo->conn[bidx] = blinp->block_connections;

    // Create a new app for each block.
    app->blocks[bidx] = gkyl_gyrokinetic_app_new(blinp);
  }

  return app;
}

void
gkyl_gyrokinetic_mb_app_apply_ic(gkyl_gyrokinetic_mb_app* app, double t0)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_apply_ic(app->blocks[bidx], t0);
  }
}

static void
calc_field(gkyl_gyrokinetic_mb_app* app, double tcurr, int fidx)
{
  // Compute fields.


  // NYI.
}

static void
calc_field_and_apply_bc(gkyl_gyrokinetic_mb_app* app, double tcurr, int fidx)
{
  // Compute fields and apply BCs.

  // Compute the field.
  calc_field(app, tcurr, fidx);

  // Apply boundary conditions.

  // NYI.

  // We need to 
  //   a) apply the BCs on physical boundaries.
  //   b) sync within each block
  //   c) sync between blocks.
  const int nblocks = app->num_blocks_local;
  for (int bidx=0; bidx<nblocks; bidx++) {
    struct gkyl_gyrokinetic_app *gkbl = app->blocks[bidx];

    struct gkyl_array *f_charged[gkbl->num_species];
    struct gkyl_array *f_neut[gkbl->num_neut_species];
    for (int i=0; i<gkbl->num_species; ++i) {
      f_charged[i] = gkbl->species[i].distfs[fidx];
    }
    for (int i=0; i<gkbl->num_neut_species; ++i) {
      if (!gkbl->neut_species[i].info.is_static) {
        f_neut[i] = gkbl->neut_species[i].distfs[fidx];
      }
    }

    gkyl_gyrokinetic_apply_bc(gkbl, f_charged, f_neut);
  }
}

static void
rk3_fill_fin_fout(const struct gkyl_gyrokinetic_app *gkbl, int fidx_in, int fidx_out,
  const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[],
  struct gkyl_array *fout[], struct gkyl_array *fout_neut[])
{
  // Fill arrays of distributions functions with the appropriate distributions
  // specified by the fidx_in and fidx_out indices.

  for (int i=0; i<gkbl->num_species; ++i) {
    fin[i] = gkbl->species[i].distfs[fidx_in];
    fout[i] = gkbl->species[i].distfs[fidx_out];
  }
  for (int i=0; i<gkbl->num_neut_species; ++i) {
    fin_neut[i] = gkbl->neut_species[i].distfs[0];
    if (!gkbl->neut_species[i].info.is_static) {
      fin_neut[i] = gkbl->neut_species[i].distfs[fidx_in];
      fout_neut[i] = gkbl->neut_species[i].distfs[fidx_out];
    }
  }
}

// Take time-step using the RK3 method. Also sets the status object
// which has the actual and suggested dts used. These can be different
// from the actual time-step.
static struct gkyl_update_status
rk3(gkyl_gyrokinetic_mb_app* app, double dt0)
{
  const int nblocks = app->num_blocks_local;

  struct gkyl_update_status st[nblocks];
  for (int bidx=0; bidx<nblocks; bidx++)
    st[bidx].success = true;

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  int fidx_in = 0, fidx_out = 1;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
      {
        // Indices in gk_species' distfs of input and output distributions.
        const int fidx_in = 0, fidx_out = 1;

        // Compute df/dt = ...
        for (int bidx=0; bidx<nblocks; bidx++) {
          struct gkyl_gyrokinetic_app *gkbl = app->blocks[bidx];

          const struct gkyl_array *fin[gkbl->num_species];
          const struct gkyl_array *fin_neut[gkbl->num_neut_species];
          struct gkyl_array *fout[gkbl->num_species];
          struct gkyl_array *fout_neut[gkbl->num_neut_species];
          rk3_fill_fin_fout(gkbl, fidx_in, fidx_out, fin, fin_neut, fout, fout_neut);

          gkyl_gyrokinetic_dfdt(gkbl, tcurr, dt, fin, fout, fin_neut, fout_neut, &st[bidx]);
        }

        // Reduce time step across blocks and take a forward euler step.
        double dt_actual = DBL_MAX, dt_suggested = DBL_MAX;
        for (int bidx=0; bidx<nblocks; bidx++) {
          dt_actual = GKYL_MIN2(dt_actual, st[bidx].dt_actual);
          dt_suggested = GKYL_MIN2(dt_suggested, st[bidx].dt_suggested);
        }

        for (int bidx=0; bidx<nblocks; bidx++) {
          st[bidx].dt_actual = dt_actual;
          st[bidx].dt_suggested = dt_suggested;

          struct gkyl_gyrokinetic_app *gkbl = app->blocks[bidx];
          const struct gkyl_array *fin[gkbl->num_species];
          const struct gkyl_array *fin_neut[gkbl->num_neut_species];
          struct gkyl_array *fout[gkbl->num_species];
          struct gkyl_array *fout_neut[gkbl->num_neut_species];
          rk3_fill_fin_fout(gkbl, fidx_in, fidx_out, fin, fin_neut, fout, fout_neut);

          for (int i=0; i<gkbl->num_species; ++i) {
            gkyl_array_accumulate(gkyl_array_scale(fout[i], dt_actual), 1.0, fin[i]);
          }
          for (int i=0; i<gkbl->num_neut_species; ++i) {
            if (!gkbl->neut_species[i].info.is_static) {
              gkyl_array_accumulate(gkyl_array_scale(fout_neut[i], dt_actual), 1.0, fin_neut[i]);
            }
          }
        }

        // Compute the fields and apply BCs.
        calc_field_and_apply_bc(app, tcurr, fidx_out);

        dt = dt_actual;
        state = RK_STAGE_2;
        break;
      }
      case RK_STAGE_2:
      {
        // Indices in gk_species' distfs of input and output distributions.
        const int fidx_in = 1, fidx_out = 2;

        for (int bidx=0; bidx<nblocks; bidx++) {
          struct gkyl_gyrokinetic_app *gkbl = app->blocks[bidx];

          const struct gkyl_array *fin[gkbl->num_species];
          const struct gkyl_array *fin_neut[gkbl->num_neut_species];
          struct gkyl_array *fout[gkbl->num_species];
          struct gkyl_array *fout_neut[gkbl->num_neut_species];
          rk3_fill_fin_fout(gkbl, fidx_in, fidx_out, fin, fin_neut, fout, fout_neut);

          gkyl_gyrokinetic_dfdt(gkbl, tcurr+dt, dt, fin, fout, fin_neut, fout_neut, &st[bidx]);
        }

        // Reduce time step across blocks and take a forward euler step.
        double dt_actual = DBL_MAX, dt_suggested = DBL_MAX;
        for (int bidx=0; bidx<nblocks; bidx++) {
          dt_actual = GKYL_MIN2(dt_actual, st[bidx].dt_actual);
          dt_suggested = GKYL_MIN2(dt_suggested, st[bidx].dt_suggested);
        }

        for (int bidx=0; bidx<nblocks; bidx++) {
          st[bidx].dt_actual = dt_actual;
          st[bidx].dt_suggested = dt_suggested;

          struct gkyl_gyrokinetic_app *gkbl = app->blocks[bidx];
          const struct gkyl_array *fin[gkbl->num_species];
          const struct gkyl_array *fin_neut[gkbl->num_neut_species];
          struct gkyl_array *fout[gkbl->num_species];
          struct gkyl_array *fout_neut[gkbl->num_neut_species];
          rk3_fill_fin_fout(gkbl, fidx_in, fidx_out, fin, fin_neut, fout, fout_neut);

          for (int i=0; i<gkbl->num_species; ++i) {
            gkyl_array_accumulate(gkyl_array_scale(fout[i], dt_actual), 1.0, fin[i]);
          }
          for (int i=0; i<gkbl->num_neut_species; ++i) {
            if (!gkbl->neut_species[i].info.is_static) {
              gkyl_array_accumulate(gkyl_array_scale(fout_neut[i], dt_actual), 1.0, fin_neut[i]);
            }
          }
        }

        if (dt_actual < dt) {

          // Recalculate the field.
          calc_field(app, tcurr, 0);

          // collect stats
          double dt_rel_diff = (dt-dt_actual)/dt_actual;
          app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_2_fail += 1;

          dt = dt_actual;
          state = RK_STAGE_1; // restart from stage 1

        } 
        else {
          for (int bidx=0; bidx<nblocks; bidx++) {
            struct gkyl_gyrokinetic_app *gkbl = app->blocks[bidx];

            for (int i=0; i<gkbl->num_species; ++i) {
              array_combine(gkbl->species[i].distfs[1],
                3.0/4.0, gkbl->species[i].distfs[0], 1.0/4.0, gkbl->species[i].distfs[2], &gkbl->species[i].local_ext);
            }
            for (int i=0; i<gkbl->num_neut_species; ++i) {
              if (!gkbl->neut_species[i].info.is_static) {
                array_combine(gkbl->neut_species[i].distfs[1],
                  3.0/4.0, gkbl->neut_species[i].distfs[0], 1.0/4.0, gkbl->neut_species[i].distfs[2], &gkbl->neut_species[i].local_ext);
              }
            }
          }

          // Compute the fields and apply BCs.
          calc_field_and_apply_bc(app, tcurr, 1);

          state = RK_STAGE_3;
        }
        break;
      }
      case RK_STAGE_3:
      {
        // Indices in gk_species' distfs of input and output distributions.
        const int fidx_in = 1, fidx_out = 2;

        for (int bidx=0; bidx<nblocks; bidx++) {
          struct gkyl_gyrokinetic_app *gkbl = app->blocks[bidx];

          const struct gkyl_array *fin[gkbl->num_species];
          const struct gkyl_array *fin_neut[gkbl->num_neut_species];
          struct gkyl_array *fout[gkbl->num_species];
          struct gkyl_array *fout_neut[gkbl->num_neut_species];
          rk3_fill_fin_fout(gkbl, fidx_in, fidx_out, fin, fin_neut, fout, fout_neut);

          gkyl_gyrokinetic_dfdt(gkbl, tcurr+dt/2, dt, fin, fout, fin_neut, fout_neut, &st[bidx]);
        }

        // Reduce time step across blocks and take a forward euler step.
        double dt_actual = DBL_MAX, dt_suggested = DBL_MAX;
        for (int bidx=0; bidx<nblocks; bidx++) {
          dt_actual = GKYL_MIN2(dt_actual, st[bidx].dt_actual);
          dt_suggested = GKYL_MIN2(dt_suggested, st[bidx].dt_suggested);
        }

        for (int bidx=0; bidx<nblocks; bidx++) {
          st[bidx].dt_actual = dt_actual;
          st[bidx].dt_suggested = dt_suggested;

          struct gkyl_gyrokinetic_app *gkbl = app->blocks[bidx];
          const struct gkyl_array *fin[gkbl->num_species];
          const struct gkyl_array *fin_neut[gkbl->num_neut_species];
          struct gkyl_array *fout[gkbl->num_species];
          struct gkyl_array *fout_neut[gkbl->num_neut_species];
          rk3_fill_fin_fout(gkbl, fidx_in, fidx_out, fin, fin_neut, fout, fout_neut);

          for (int i=0; i<gkbl->num_species; ++i) {
            gkyl_array_accumulate(gkyl_array_scale(fout[i], dt_actual), 1.0, fin[i]);
          }
          for (int i=0; i<gkbl->num_neut_species; ++i) {
            if (!gkbl->neut_species[i].info.is_static) {
              gkyl_array_accumulate(gkyl_array_scale(fout_neut[i], dt_actual), 1.0, fin_neut[i]);
            }
          }
        }

        if (dt_actual < dt) {
          // Recalculate the field.
          calc_field(app, tcurr, 0);

          // collect stats
          double dt_rel_diff = (dt-dt_actual)/dt_actual;
          app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_3_fail += 1;

          dt = dt_actual;
          state = RK_STAGE_1; // restart from stage 1

          app->stat.nstage_2_fail += 1;
        }
        else {
          for (int bidx=0; bidx<nblocks; bidx++) {
            struct gkyl_gyrokinetic_app *gkbl = app->blocks[bidx];

            for (int i=0; i<gkbl->num_species; ++i) {
              array_combine(gkbl->species[i].distfs[1],
                1.0/3.0, gkbl->species[i].distfs[0], 2.0/3.0, gkbl->species[i].distfs[2], &gkbl->species[i].local_ext);
              gkyl_array_copy_range(gkbl->species[i].distfs[0], gkbl->species[i].distfs[1], &gkbl->species[i].local_ext);
            }
            for (int i=0; i<gkbl->num_neut_species; ++i) {
              if (!gkbl->neut_species[i].info.is_static) {
                array_combine(gkbl->neut_species[i].distfs[1],
                  1.0/3.0, gkbl->neut_species[i].distfs[0], 2.0/3.0, gkbl->neut_species[i].distfs[2], &gkbl->neut_species[i].local_ext);
                gkyl_array_copy_range(gkbl->neut_species[i].distfs[0], gkbl->neut_species[i].distfs[1], &gkbl->neut_species[i].local_ext);
              }
            }
          }

          // Compute the fields and apply BCs
          calc_field_and_apply_bc(app, tcurr, 0);

          state = RK_COMPLETE;
        }
        break;
      }
      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }

  return st[0];
}

struct gkyl_update_status
gkyl_gyrokinetic_mb_update(gkyl_gyrokinetic_mb_app* app, double dt)
{
  app->stat.nup += 1;
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_update_status status[GKYL_MAX_BLOCKS];
  status[0] = rk3(app, dt);
  app->tcurr += status[0].dt_actual;

  app->stat.total_tm += gkyl_time_diff_now_sec(wst);
  // Check for any CUDA errors during time step
  if (app->use_gpu)
    checkCuda(cudaGetLastError());
//  struct gkyl_update_status status[GKYL_MAX_BLOCKS];
//  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
//    status[bidx] = gkyl_gyrokinetic_update(app->blocks[bidx], dt);
//  }
  return status[0];
}

void
gkyl_gyrokinetic_mb_app_calc_mom(gkyl_gyrokinetic_mb_app* app)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_calc_mom(app->blocks[bidx]);
  }
}

void
gkyl_gyrokinetic_mb_app_calc_integrated_mom(gkyl_gyrokinetic_mb_app* app, double tm)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_calc_integrated_mom(app->blocks[bidx], tm);
  }
}

void
gkyl_gyrokinetic_mb_app_write_mom(gkyl_gyrokinetic_mb_app* app, double tm, int frame)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_write_mom(app->blocks[bidx], tm, frame);
  }
}

void
gkyl_gyrokinetic_mb_app_write_integrated_mom(gkyl_gyrokinetic_mb_app* app)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_write_integrated_mom(app->blocks[bidx]);
  }
}

void
gkyl_gyrokinetic_mb_app_write(gkyl_gyrokinetic_mb_app* app, double tm, int frame)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_write(app->blocks[bidx], tm, frame);
  }
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_mb_app_read_from_frame(gkyl_gyrokinetic_mb_app *app, int frame)
{
  // NYI
  struct gkyl_app_restart_status stat;
  return stat;
}

void
gkyl_gyrokinetic_mb_app_stat_write(gkyl_gyrokinetic_mb_app* app)
{
  // NYI
}

struct gkyl_gyrokinetic_stat
gkyl_gyrokinetic_mb_app_stat(gkyl_gyrokinetic_mb_app* app)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gk_species_tm(app->blocks[bidx]);
    gk_species_coll_tm(app->blocks[bidx]);
  }
  return app->blocks[0]->stat;
}

// private function to handle variable argument list for printing
static void
v_gk_app_cout(const gkyl_gyrokinetic_mb_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank, r = 0;
  gkyl_comm_get_rank(app->blocks[0]->comm, &rank);
  if ((rank == 0) && fp)
    vfprintf(fp, fmt, argp);
}

void
gkyl_gyrokinetic_mb_app_cout(const gkyl_gyrokinetic_mb_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_gk_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_gyrokinetic_mb_app_release(gkyl_gyrokinetic_mb_app* app)
{
  for (int bidx=0; bidx<app->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_release(app->blocks[bidx]);
  }

  gkyl_block_topo_release(app->btopo);
  gkyl_free(app->blocks);

  gkyl_free(app);
}
