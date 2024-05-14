#include <gkyl_amr_block_priv.h>

void
skin_ghost_ranges_init_block(struct skin_ghost_ranges_block* sgr, const struct gkyl_range* parent, const int* ghost)
{
  int ndim = parent->ndim;

  for (int d = 0; d < ndim; d++) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d], d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d], d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void
euler_transmissive_bc(double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  for (int i = 0; i < 5; i++) {
    ghost[i] = skin[i];
  }
}

void
euler_block_bc_updaters_init(struct euler_block_data* bdata, const struct gkyl_block_connections* conn)
{
  int nghost[9];
  for (int i = 0; i < 9; i++) {
    nghost[i] = 2;
  }

  for (int d = 0; d < 2; d++) {
    bdata->lower_bc[d] = bdata->upper_bc[d] = 0;

    if (conn->connections[d][0].edge == GKYL_PHYSICAL) {
      bdata->lower_bc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler, bdata->geom, d, GKYL_LOWER_EDGE, nghost,
        euler_transmissive_bc, 0);
    }

    if (conn->connections[d][1].edge == GKYL_PHYSICAL) {
      bdata->upper_bc[d] = gkyl_wv_apply_bc_new(&bdata->grid, bdata->euler, bdata->geom, d, GKYL_UPPER_EDGE, nghost,
        euler_transmissive_bc, 0);
    }
  }

  skin_ghost_ranges_init_block(&bdata->skin_ghost, &bdata->ext_range, nghost);
  long buff_sz = 0;

  for (int d = 0; d < 2; d++) {
    long vol = bdata->skin_ghost.lower_skin[d].volume;

    if (buff_sz <= vol) {
      buff_sz = vol;
    }
  }

  bdata->bc_buffer = gkyl_array_new(GKYL_DOUBLE, 5, buff_sz);
}

void
euler_block_bc_updaters_release(struct euler_block_data* bdata)
{
  for (int d = 0; d < 2; d++) {
    if (bdata->lower_bc[d]) {
      gkyl_wv_apply_bc_release(bdata->lower_bc[d]);
    }

    if (bdata->upper_bc[d]) {
      gkyl_wv_apply_bc_release(bdata->upper_bc[d]);
    }
  }

  gkyl_array_release(bdata->bc_buffer);
}

void
euler_block_bc_updaters_apply(const struct euler_block_data* bdata, double tm, struct gkyl_array* fld)
{
  for (int d = 0; d < 2; d++) {
    if (bdata->lower_bc[d]) {
      gkyl_wv_apply_bc_advance(bdata->lower_bc[d], tm, &bdata->range, fld);
    }

    if (bdata->upper_bc[d]) {
      gkyl_wv_apply_bc_advance(bdata->upper_bc[d], tm, &bdata->range, fld);
    }
  }
}

void
euler_sync_blocks(const struct gkyl_block_topo* btopo, const struct euler_block_data bdata[], struct gkyl_array* fld[])
{
  int num_blocks = btopo->num_blocks;
  int ndim = btopo->ndim;
  
  for (int i = 0; i < num_blocks; i++) {
    for (int d = 0; d < ndim; d++) {
      const struct gkyl_target_edge *te = btopo->conn[i].connections[d];

#ifdef AMR_DEBUG
      if (te[0].edge != GKYL_PHYSICAL) {
        struct gkyl_array *bc_buffer = bdata[i].bc_buffer;

        gkyl_array_copy_to_buffer(bc_buffer->data, fld[i], &(bdata[i].skin_ghost.lower_skin[d]));

        int tbid = te[0].bid;
        int tdir = te[0].dir;

        if (te[0].edge == GKYL_LOWER_POSITIVE) {
          gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
        }
        else if (te[0].edge == GKYL_UPPER_POSITIVE) {
          gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
        }
      }

      if (te[1].edge != GKYL_PHYSICAL) {
        struct gkyl_array *bc_buffer = bdata[i].bc_buffer;

        gkyl_array_copy_to_buffer(bc_buffer->data, fld[i], &(bdata[i].skin_ghost.upper_skin[d]));

        int tbid = te[1].bid;
        int tdir = te[1].dir;

        if (te[1].edge == GKYL_LOWER_POSITIVE) {
          gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
        }
        else if (te[1].edge == GKYL_UPPER_POSITIVE) {
          gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
        }
      }
#else
      if (te[0].edge != GKYL_PHYSICAL) {
        struct gkyl_array *bc_buffer = bdata[i].bc_buffer;

        gkyl_array_copy_to_buffer(bc_buffer->data, fld[i], &(bdata[i].skin_ghost.lower_skin[d]));

        int tbid = te[0].bid;
        int tdir = te[0].dir;

        if (te[0].edge == GKYL_LOWER_POSITIVE) {
          if (bdata[i].skin_ghost.lower_skin[d].volume == bdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
            gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
          }
          else if (bdata[i].skin_ghost.lower_skin[d].volume > bdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
            struct gkyl_range_iter iter;
            gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));

            int ref_factor = (int)(bdata[i].skin_ghost.lower_skin[d].volume / bdata[tbid].skin_ghost.lower_ghost[tdir].volume);

            long count = 0;
            while (gkyl_range_iter_next(&iter)) {
              long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);

              if ((bdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * (ref_factor * count++), fld[tbid]->esznc);
              }
              else {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld[tbid]->esznc);
                count += 1;
              }
            }
          }
          else if (bdata[i].skin_ghost.lower_skin[d].volume < bdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
            struct gkyl_range_iter iter;
            gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));

            double ref_factor_inv = ((double)bdata[i].skin_ghost.lower_skin[d].volume / (double)bdata[tbid].skin_ghost.lower_ghost[tdir].volume);

            long count = 0;
            while (gkyl_range_iter_next(&iter)) {
              long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);

              if ((bdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * ((int)(ref_factor_inv * count++)), fld[tbid]->esznc);
              }
              else {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld[tbid]->esznc);
                count += 1;
              }
            }
          }
        }
        else if (te[0].edge == GKYL_UPPER_POSITIVE) {
          if (bdata[i].skin_ghost.lower_skin[d].volume == bdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
            gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
          }
          else if (bdata[i].skin_ghost.lower_skin[d].volume > bdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
            struct gkyl_range_iter iter;
            gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));

            int ref_factor = (int)(bdata[i].skin_ghost.lower_skin[d].volume / bdata[tbid].skin_ghost.upper_ghost[tdir].volume);

            long count = 0;
            while (gkyl_range_iter_next(&iter)) {
              long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);

              if ((bdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * (ref_factor * count++), fld[tbid]->esznc);
              }
              else {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld[tbid]->esznc);
                count += 1;
              }
            }
          }
          else if (bdata[i].skin_ghost.lower_skin[d].volume < bdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
            struct gkyl_range_iter iter;
            gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));

            double ref_factor_inv = ((double)bdata[i].skin_ghost.lower_skin[d].volume / (double)bdata[tbid].skin_ghost.upper_ghost[tdir].volume);

            long count = 0;
            while (gkyl_range_iter_next(&iter)) {
              long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);

              if ((bdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * ((int)(ref_factor_inv * count++)), fld[tbid]->esznc);
              }
              else {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld[tbid]->esznc);
                count += 1;
              }
            }
          }
        }
      }

      if (te[1].edge != GKYL_PHYSICAL) {
        struct gkyl_array *bc_buffer = bdata[i].bc_buffer;

        gkyl_array_copy_to_buffer(bc_buffer->data, fld[i], &(bdata[i].skin_ghost.upper_skin[d]));

        int tbid = te[1].bid;
        int tdir = te[1].dir;

        if (te[1].edge == GKYL_LOWER_POSITIVE) {
          if (bdata[i].skin_ghost.upper_skin[d].volume == bdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
            gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));
          }
          else if (bdata[i].skin_ghost.upper_skin[d].volume > bdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
            struct gkyl_range_iter iter;
            gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));

            int ref_factor = (int)(bdata[i].skin_ghost.upper_skin[d].volume / bdata[tbid].skin_ghost.lower_ghost[tdir].volume);

            long count = 0;
            while (gkyl_range_iter_next(&iter)) {
              long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);

              if ((bdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * (ref_factor * count++), fld[tbid]->esznc);
              }
              else {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld[tbid]->esznc);
                count += 1;
              }
            }
          }
          else if (bdata[i].skin_ghost.upper_skin[d].volume < bdata[tbid].skin_ghost.lower_ghost[tdir].volume) {
            struct gkyl_range_iter iter;
            gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.lower_ghost[tdir]));

            double ref_factor_inv = ((double)bdata[i].skin_ghost.upper_skin[d].volume / (double)bdata[tbid].skin_ghost.lower_ghost[tdir].volume);

            long count = 0;
            while (gkyl_range_iter_next(&iter)) {
              long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.lower_ghost[tdir]), iter.idx);
              
              if ((bdata[tbid].skin_ghost.lower_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.lower_ghost[tdir].lower[0]) == 1) {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * ((int)(ref_factor_inv * count++)), fld[tbid]->esznc);
              }
              else {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld[tbid]->esznc);
                count += 1;
              }
            }
          }
        }
        else if (te[1].edge == GKYL_UPPER_POSITIVE) {
          if (bdata[i].skin_ghost.upper_skin[d].volume == bdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
            gkyl_array_copy_from_buffer(fld[tbid], bc_buffer->data, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));
          }
          else if (bdata[i].skin_ghost.upper_skin[d].volume > bdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
            struct gkyl_range_iter iter;
            gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));

            int ref_factor = (int)(bdata[i].skin_ghost.upper_skin[d].volume / bdata[tbid].skin_ghost.upper_ghost[tdir].volume);

            long count = 0;
            while (gkyl_range_iter_next(&iter)) {
              long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);
              
              if ((bdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * (ref_factor * count++), fld[tbid]->esznc);
              }
              else {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * ((ref_factor * (count - (count % 2))) + (count % 2)), fld[tbid]->esznc);
                count += 1;
              }
            }
          }
          else if (bdata[i].skin_ghost.upper_skin[d].volume < bdata[tbid].skin_ghost.upper_ghost[tdir].volume) {
            struct gkyl_range_iter iter;
            gkyl_range_iter_init(&iter, &(bdata[tbid].skin_ghost.upper_ghost[tdir]));

            double ref_factor_inv = ((double)bdata[i].skin_ghost.upper_skin[d].volume / (double)bdata[tbid].skin_ghost.upper_ghost[tdir].volume);

            long count = 0;
            while (gkyl_range_iter_next(&iter)) {
              long start = gkyl_range_idx(&(bdata[tbid].skin_ghost.upper_ghost[tdir]), iter.idx);

              if ((bdata[tbid].skin_ghost.upper_ghost[tdir].upper[0] - bdata[tbid].skin_ghost.upper_ghost[tdir].lower[0]) == 1) {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * ((int)(ref_factor_inv * count++)), fld[tbid]->esznc);
              }
              else {
                memcpy(gkyl_array_fetch(fld[tbid], start),
                  ((char*) bc_buffer->data) + fld[tbid]->esznc * (2 * (int)(0.5 * ref_factor_inv * count) + (count % 2)), fld[tbid]->esznc);
                count += 1;
              }
            }
          }
        }
      }
#endif
    }
  }
}

void
euler_block_data_write(const char* file_nm, const struct euler_block_data* bdata)
{
  gkyl_grid_sub_array_write(&bdata->grid, &bdata->range, bdata->f[0], file_nm);
}

double
euler_block_data_max_dt(const struct euler_block_data* bdata)
{
  double dt = DBL_MAX;

  for (int d = 0; d < 2; d++) {
    dt = fmin(dt, gkyl_wave_prop_max_dt(bdata->slvr[d], &bdata->range, bdata->f[0]));
  }

  return dt;
}

void
euler_update_block_job_func(void* ctx)
{
  struct euler_update_block_ctx *ub_ctx = ctx;
  const struct euler_block_data *bdata = ub_ctx->bdata;

  int d = ub_ctx->dir;
  double t_curr = ub_ctx->t_curr;
  double dt = ub_ctx->dt;

  ub_ctx->stat = gkyl_wave_prop_advance(bdata->slvr[d], t_curr, dt, &bdata->range, bdata->f[d], bdata->f[d + 1]);

  euler_block_bc_updaters_apply(bdata, t_curr, bdata->f[d + 1]);
}

struct gkyl_update_status
euler_update_all_blocks(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo,
  const struct euler_block_data bdata[], double t_curr, double dt)
{
  int num_blocks  = btopo->num_blocks;
  int ndim = btopo->ndim;

  double dt_suggested = DBL_MAX;

  for (int d = 0; d < ndim; d++) {
    struct euler_update_block_ctx euler_block_ctx[num_blocks];
    
    for (int i = 0; i < num_blocks; i++) {
      euler_block_ctx[i] = (struct euler_update_block_ctx) {
        .bdata = &bdata[i],
        .t_curr = t_curr,
        .dir = d,
        .dt = dt,
        .bidx = i,
      };
    }

#ifdef AMR_USETHREADS
    for (int i = 0; i < num_blocks; i++) {
      gkyl_job_pool_add_work(job_pool, euler_update_block_job_func, &euler_block_ctx[i]);
    }
    gkyl_job_pool_wait(job_pool);
#else
    for (int i = 0; i < num_blocks; i++) {
      euler_update_block_job_func(&euler_block_ctx[i]);
    }
#endif

    struct gkyl_array *fld[num_blocks];

    for (int i = 0; i < num_blocks; i++) {
      if (euler_block_ctx[i].stat.success == false) {
        return (struct gkyl_update_status) {
          .success = false,
          .dt_suggested = euler_block_ctx[i].stat.dt_suggested,
        };
      }

      dt_suggested = fmin(dt_suggested, euler_block_ctx[i].stat.dt_suggested);
      fld[i] = bdata[i].f[d + 1];
    }

    euler_sync_blocks(btopo, bdata, fld);
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested,
  };
}

void
euler_init_job_func_block(void* ctx)
{
  struct euler_block_data *bdata = ctx;

  gkyl_fv_proj_advance(bdata->fv_proj, 0.0, &bdata->ext_range, bdata->f[0]);
}

void
copy_job_func(void* ctx)
{
  struct copy_job_ctx *j_ctx = ctx;

  gkyl_array_copy(j_ctx->out, j_ctx->inp);
}

struct gkyl_update_status
euler_update_block(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo,
  const struct euler_block_data bdata[], double t_curr, double dt0, struct sim_stats* stats)
{
  int num_blocks = btopo->num_blocks;
  double dt_suggested = DBL_MAX;

  enum {
    UPDATE_DONE = 0,
    PRE_UPDATE,
    POST_UPDATE,
    FLUID_UPDATE,
    UPDATE_REDO,
  } state = PRE_UPDATE;

  struct copy_job_ctx euler_copy_ctx[num_blocks];
  double dt = dt0;

  while (state != UPDATE_DONE) {
    if (state == PRE_UPDATE) {
      state = FLUID_UPDATE;

      for (int i = 0; i < num_blocks; i++) {
        euler_copy_ctx[i] = (struct copy_job_ctx) {
          .bidx = i,
          .inp = bdata[i].f[0],
          .out = bdata[i].fdup,
        };
      }

#ifdef AMR_USETHREADS
      for (int i = 0; i < num_blocks; i++) {
        gkyl_job_pool_add_work(job_pool, copy_job_func, &euler_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
#else
      for (int i = 0; i < num_blocks; i++) {
        copy_job_func(&euler_copy_ctx[i]);
      }
#endif
    }
    else if (state == FLUID_UPDATE) {
      state = POST_UPDATE;

      struct gkyl_update_status s = euler_update_all_blocks(job_pool, btopo, bdata, t_curr, dt);

      if (!s.success) {
        stats->nfail += 1;
        dt = s.dt_suggested;
        state = UPDATE_REDO;
      }
      else {
        dt_suggested = fmin(dt_suggested, s.dt_suggested);
      }
    }
    else if (state == POST_UPDATE) {
      state = UPDATE_DONE;

      for (int i = 0; i < num_blocks; i++) {
        euler_copy_ctx[i] = (struct copy_job_ctx) {
          .bidx = i,
          .inp = bdata[i].f[2],
          .out = bdata[i].f[0],
        };
      }

#ifdef AMR_USETHREADS
      for (int i = 0; i < num_blocks; i++) {
        gkyl_job_pool_add_work(job_pool, copy_job_func, &euler_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
#else
      for (int i = 0; i < num_blocks; i++) {
        copy_job_func(&euler_copy_ctx[i]);
      }
#endif
    }
    else if (state == UPDATE_REDO) {
      state = PRE_UPDATE;

      for (int i = 0; i < num_blocks; i++) {
        euler_copy_ctx[i] = (struct copy_job_ctx) {
          .bidx = i,
          .inp = bdata[i].fdup,
          .out = bdata[i].f[0],
        };
      }

#ifdef AMR_USETHREADS
      for (int i = 0; i < num_blocks; i++) {
        gkyl_job_pool_add_work(job_pool, copy_job_func, &euler_copy_ctx[i]);
      }
      gkyl_job_pool_wait(job_pool);
#else
      for (int i = 0; i < num_blocks; i++) {
        copy_job_func(&euler_copy_ctx[i]);
      }
#endif
    }
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_actual = dt,
    .dt_suggested = dt_suggested,
  };
}

void
euler_write_sol_block(const char* fbase, int num_blocks, const struct euler_block_data bdata[])
{
  for (int i = 0; i < num_blocks; i++) {
    const char *fmt = "%s_b%d.gkyl";
    int sz = snprintf(0, 0, fmt, fbase, i);
    char file_nm[sz + 1];

    snprintf(file_nm, sizeof file_nm, fmt, fbase, i);
    euler_block_data_write(file_nm, &bdata[i]);
  }
}

double
euler_max_dt_block(int num_blocks, const struct euler_block_data bdata[])
{
  double dt = DBL_MAX;

  for (int i = 0; i < num_blocks; i++) {
    dt = fmin(dt, euler_block_data_max_dt(&bdata[i]));
  }

  return dt;
}

struct gkyl_block_topo*
create_block_topo()
{
  struct gkyl_block_topo *btopo = gkyl_block_topo_new(2, 9);

  btopo->conn[0] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 4, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 5, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
    .connections[1] = { { .bid = 7, .dir = 1, .edge = GKYL_UPPER_POSITIVE }, { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE } },
  };

  btopo->conn[1] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
    .connections[1] = { { .bid = 4, .dir = 1, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } },
  };

  btopo->conn[2] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 3, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
    .connections[1] = { { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } },
  };

  btopo->conn[3] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 2, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } },
    .connections[1] = { { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } },
  };

  btopo->conn[4] = (struct gkyl_block_connections) {
   .connections[0] = { { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, { .bid = 0, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
   .connections[1] = { { .bid = 6, .dir = 1, .edge = GKYL_UPPER_POSITIVE }, { .bid = 1, .dir = 1, .edge = GKYL_LOWER_POSITIVE } },
  };

  btopo->conn[5] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } },
    .connections[1] = { { .bid = 8, .dir = 1, .edge = GKYL_UPPER_POSITIVE }, { .bid = 3, .dir = 1, .edge = GKYL_LOWER_POSITIVE } },
  };

  btopo->conn[6] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, { .bid = 7, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
    .connections[1] = { { .bid = 6, .dir = 1, .edge = GKYL_PHYSICAL }, { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE } },
  };

  btopo->conn[7] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 6, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 8, .dir = 0, .edge = GKYL_LOWER_POSITIVE } },
    .connections[1] = { { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, { .bid = 0, .dir = 1, .edge = GKYL_LOWER_POSITIVE } },
  };

  btopo->conn[8] = (struct gkyl_block_connections) {
    .connections[0] = { { .bid = 7, .dir = 0, .edge = GKYL_UPPER_POSITIVE }, { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } },
    .connections[1] = { { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, { .bid = 5, .dir = 1, .edge = GKYL_LOWER_POSITIVE } },
  };

  return btopo;
}