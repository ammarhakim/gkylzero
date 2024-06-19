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

  app->num_blocks = inp->num_blocks;

  app->btopo = gkyl_block_topo_new(app->cdim, app->num_blocks);
  for (int bc=0; bc<app->num_blocks; bc++) {
    struct gkyl_gk *blinp = inp->blocks[bc];
    app->btopo->conn[bc] = blinp->block_connections;
  }

  app->use_mpi = false;
  app->use_gpu = false;
#ifdef GKYL_HAVE_MPI
  app->use_mpi = inp->use_mpi;
#endif
#ifdef GKYL_HAVE_CUDA
  app->use_gpu = inp->use_gpu;
#endif

  // Construct MB communicator.
#ifdef GKYL_HAVE_MPI
  if (app->use_gpu && app->use_mpi) {
#ifdef GKYL_HAVE_NCCL
    app->comm_mb = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) { .mpi_comm = MPI_COMM_WORLD, } );
#else
    fprintf(stderr, " Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (app->use_mpi) {
    app->comm_mb = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) { .mpi_comm = MPI_COMM_WORLD, } );
  }
  else {
    app->comm_mb = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { .use_gpu = app->use_gpu } );
  }
#else
  app->comm_mb = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { .use_gpu = app->use_gpu } );
#endif

  int comm_size, comm_rank;
  gkyl_comm_get_size(app->comm_mb, &comm_size);
  gkyl_comm_get_rank(app->comm_mb, &comm_rank);

  int decomp_vol = 0;
  for (int bidx=0; bidx<app->num_blocks; bidx++) {
    int cuts_vol = 1;
    for (int d=0; d<app->cdim; d++) cuts_vol *= inp->blocks[bidx]->cuts[d];
    decomp_vol += cuts_vol;
  }

  if (decomp_vol == comm_size) {

    // Each rank owns a single block, or a subdomain of one.
    app->num_blocks_local = 1;
    app->block_idxs = gkyl_malloc(app->num_blocks_local*sizeof(int));
    app->block_idxs[0] = -1;
    int proc_count = 0;
    for (int bidx=0; bidx<app->num_blocks; bidx++) {
      int cuts_vol = 1;
      for (int d=0; d<app->cdim; d++) cuts_vol *= inp->blocks[bidx]->cuts[d];
      proc_count += cuts_vol;

      if (app->block_idxs[0] < 0 && comm_rank < proc_count)
        app->block_idxs[0] = bidx;
    }

  }
  else {
    // This case is intended for simulations with multiple single-cut blocks
    // (scbs, w/ cuts_vol=1). Blocks with cuts_vol>1 will have one rank owning
    // a subdomain of those blocks, and that rank will own nothing else. Yet a
    // rank may handle one or more scbs.

    assert(decomp_vol > comm_size); // Can't have more ranks than decompositions.

    int *cuts_vol_per_block = gkyl_malloc(app->num_blocks * sizeof(int));
    int num_scb = 0; // Number of single-cut blocks (scb).
    int scb[GKYL_MAX_BLOCKS]; // Block ID of single-cut blocks (scb)
    int cuts_vol_max = -1;
    for (int bidx=0; bidx<app->num_blocks; bidx++) {
      int cuts_vol = 1;
      for (int d=0; d<app->cdim; d++) cuts_vol *= inp->blocks[bidx]->cuts[d];

      cuts_vol_per_block[bidx] = cuts_vol;
      cuts_vol_max = GKYL_MAX2(cuts_vol_max, cuts_vol);

      if (cuts_vol == 1) {
        scb[num_scb] = bidx;
        num_scb++;
      }
    }

    // Additional blocks that need to be assigned to a rank owning another block.
    int extra_blocks = decomp_vol - comm_size;
    // Number of ranks owning single-cut blocks (scbrank).
    int num_scbrank = num_scb - extra_blocks;
    // Distribute scb amongst scbranks: 
    int *scbrank_num_blocks = gkyl_malloc(num_scbrank * sizeof(int));
    int base = num_scb/num_scbrank;
    int rem = num_scb - num_scbrank * base;
    for (int i=0; i<num_scbrank; i++)
      scbrank_num_blocks[i] = base;
    for (int i=0; i<rem; i++)
      scbrank_num_blocks[i]++;

    // List of block IDs owned by each scbrank.
    int *scbrank_blocks = gkyl_malloc(num_scbrank * (base+1) * sizeof(int));
    int curr_scb_idx = 0;
    for (int i=0; i<num_scbrank; i++) {
      int scb_count = 0;
      for (int j=0; j<scbrank_num_blocks[i]; j++) {
        scbrank_blocks[i*(base+1)+j] = scb[scb_count];
        scb_count++;
      }
    }

    // List of ranks in each block.
    int *ranks_per_block = gkyl_malloc(app->num_blocks * cuts_vol_max * sizeof(int));
    int curr_rank_to_assign = 0, curr_scbrank = 0, scbrank_count = 0;
    for (int bidx=0; bidx<app->num_blocks; bidx++) {
      if (cuts_vol_per_block[bidx] > 1) {
        for (int i=0; i<cuts_vol_per_block[bidx]; i++) {
          ranks_per_block[bidx*cuts_vol_max+i] = curr_rank_to_assign;
          curr_rank_to_assign++;
        }
      }
      else {
        if (scbrank_count == 0) {
          // Assigning the first scb to this scbrank.
          ranks_per_block[bidx*cuts_vol_max] = curr_rank_to_assign;
          curr_rank_to_assign++;
        }
        else {
          // Assign this additional scb to a scbrank that already has one or more scbs.
          int scb_idx = scbrank_blocks[curr_scbrank*(base+1)+scbrank_count];
          ranks_per_block[bidx*cuts_vol_max] = ranks_per_block[scb_idx*cuts_vol_max];
        }
        scbrank_count++;

        if (scbrank_count == scbrank_num_blocks[curr_scbrank]) {
          scbrank_count = 0;
          curr_scbrank++;
        }
      }
    }

    // Now use the list of ranks per block to populate
    // the number and IDs of blocks owned by this rank.
    app->num_blocks_local = 0;
    int my_block_idxs[GKYL_MAX_BLOCKS];
    for (int bidx=0; bidx<app->num_blocks; bidx++) {
      for (int i=0; i<cuts_vol_per_block[bidx]; i++) {
        if (comm_rank == ranks_per_block[bidx*cuts_vol_max+i]) {
          my_block_idxs[app->num_blocks_local] = bidx;
          app->num_blocks_local++;
        }
      }
    }
    app->block_idxs = gkyl_malloc(app->num_blocks_local*sizeof(int));
    for (int i=0; i<app->num_blocks_local; i++)
      app->block_idxs[i] = my_block_idxs[i];

    gkyl_free(ranks_per_block);
    gkyl_free(scbrank_blocks);
    gkyl_free(scbrank_num_blocks);
    gkyl_free(cuts_vol_per_block);
  }

  // Allocate memory for all the block decompositions because the cross-block
  // field object needs to know the range of every block.
  app->decomp_intrab = gkyl_malloc(app->num_blocks * sizeof(struct gkyl_rect_decomp *));

  for (int bc=0; bc<app->num_blocks; bc++) {
    struct gkyl_gk *blinp = inp->blocks[bc];

    // Create intra block decompositions.
    struct gkyl_range global_range_conf;
    gkyl_create_global_range(app->cdim, blinp->cells, &global_range_conf);
    app->decomp_intrab[bc] = gkyl_rect_decomp_new_from_cuts(app->cdim, blinp->cuts, &global_range_conf);
  }

  // Only allocate memory for local blocks and their intrablock communicators.
  app->blocks = gkyl_malloc(app->num_blocks_local * sizeof(struct gkyl_gyrokinetic_app));
  app->comm_intrab = gkyl_malloc(app->num_blocks_local * sizeof(struct gkyl_comm *));

  for (int bc=0; bc<app->num_blocks_local; bc++) {
    int bidx = app->block_idxs[bc];
    struct gkyl_gk *blinp = inp->blocks[bidx];

    // Create intra blok communicators.
    int comm_color = bidx;
    struct gkyl_comm *parent_comm = bc == 0? app->comm_mb : app->comm_intrab[0];
    app->comm_intrab[bc] = gkyl_comm_split_comm(parent_comm, comm_color, app->decomp_intrab[bidx]);

    int comm_intrab_rank;
    gkyl_comm_get_rank(app->comm_intrab[bc], &comm_intrab_rank);

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

    blinp->has_low_inp = true,
    blinp->low_inp.local_range = app->decomp_intrab[bidx]->ranges[comm_intrab_rank],
    blinp->low_inp.comm = app->comm_intrab[bc],

    // Create a new app for each block.
    app->blocks[bc] = gkyl_gyrokinetic_app_new(blinp);
    app->field = gk_field_mb_new(inp, app, app->blocks[bc]);
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
  for (int bc=0; bc<app->num_blocks_local; bc++) {
    gk_field_mb_rhs(app, app->field, app->blocks[bc]);
  }
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
  for (int i=0; i<app->num_blocks_local; i++) {
    gkyl_gyrokinetic_app_release(app->blocks[i]);
  }
  gkyl_free(app->blocks);
  gkyl_free(app->block_idxs);
  
  // Release decomp and comm.
  for (int i=0; i<app->num_blocks; i++) {
    gkyl_rect_decomp_release(app->decomp_intrab[i]);
  }
  for (int i=0; i<app->num_blocks_local; i++) {
    gkyl_comm_release(app->comm_intrab[i]);
  }
  gkyl_free(app->decomp_intrab);
  gkyl_free(app->comm_intrab);

  gkyl_comm_release(app->comm_mb);

  gkyl_block_topo_release(app->btopo);
  gkyl_free(app);
}
