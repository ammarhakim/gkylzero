#include <gkyl_alloc.h>

#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>
#include <gkyl_dflt.h>

gkyl_gyrokinetic_multib_app*
gkyl_gyrokinetic_multib_app_new(struct gkyl_gk_multib *inp)
{
  disable_denorm_float();

  gkyl_gyrokinetic_multib_app *mba = gkyl_malloc(sizeof(gkyl_gyrokinetic_multib_app));

  strcpy(mba->name, inp->name);

  mba->cdim = inp->cdim;
  mba->vdim = inp->vdim;
  mba->poly_order = inp->poly_order;

  mba->cfl = inp->cfl_frac == 0 ? 1.0 : inp->cfl_frac;

  mba->num_blocks = inp->num_blocks;

  mba->btopo = gkyl_block_topo_new(mba->cdim, mba->num_blocks);
  for (int bc=0; bc<mba->num_blocks; bc++) {
    struct gkyl_gk *blinp = inp->blocks[bc];
    mba->btopo->conn[bc] = blinp->block_connections;
  }

  mba->use_mpi = false;
  mba->use_gpu = false;
#ifdef GKYL_HAVE_MPI
  mba->use_mpi = inp->use_mpi;
#endif
#ifdef GKYL_HAVE_CUDA
  mba->use_gpu = inp->use_gpu;
#endif

  // Construct MB communicator.
#ifdef GKYL_HAVE_MPI
  if (mba->use_gpu && mba->use_mpi) {
#ifdef GKYL_HAVE_NCCL
    mba->comm_multib = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) { .mpi_comm = MPI_COMM_WORLD, } );
#else
    fprintf(stderr, " Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (mba->use_mpi) {
    mba->comm_multib = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) { .mpi_comm = MPI_COMM_WORLD, } );
  }
  else {
    mba->comm_multib = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { .use_gpu = mba->use_gpu } );
  }
#else
  mba->comm_multib = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { .use_gpu = mba->use_gpu } );
#endif

  int comm_size, comm_rank;
  gkyl_comm_get_size(mba->comm_multib, &comm_size);
  gkyl_comm_get_rank(mba->comm_multib, &comm_rank);

  int cuts_vol_max = -1, cuts_vol_tot = 0;
  for (int bidx=0; bidx<mba->num_blocks; bidx++) {
    int cuts_vol = 1;
    for (int d=0; d<mba->cdim; d++) cuts_vol *= inp->blocks[bidx]->cuts[d];
    cuts_vol_tot += cuts_vol;
    cuts_vol_max = GKYL_MAX2(cuts_vol_max, cuts_vol);
  }

  int *cuts_vol_per_block = gkyl_malloc(mba->num_blocks * sizeof(int));
  int *ranks_per_block = gkyl_malloc(mba->num_blocks * cuts_vol_max * sizeof(int));
  for (int i = 0 ; i < mba->num_blocks * cuts_vol_max;i++) ranks_per_block[i] = -1;

  if (cuts_vol_tot == comm_size) {

    // Each rank owns a single block, or a subdomain of one.
    mba->num_blocks_local = 1;
    mba->block_idxs = gkyl_malloc(mba->num_blocks_local*sizeof(int));
    mba->block_idxs[0] = -1;
    int proc_count = 0;
    for (int bidx=0; bidx<mba->num_blocks; bidx++) {
      int cuts_vol = 1;
      for (int d=0; d<mba->cdim; d++) cuts_vol *= inp->blocks[bidx]->cuts[d];
      proc_count += cuts_vol;

      cuts_vol_per_block[bidx] = cuts_vol;

      if (mba->block_idxs[0] < 0 && comm_rank < proc_count)
        mba->block_idxs[0] = bidx;

    }

    int curr_rank_to_assign = 0;
    for (int bidx=0; bidx<mba->num_blocks; bidx++) {
      for (int i=0; i<cuts_vol_per_block[bidx]; i++) {
        ranks_per_block[bidx*cuts_vol_max+i] = curr_rank_to_assign;
        curr_rank_to_assign++;
      }
    }

  }
  else {
    // This case is intended for simulations with multiple single-cut blocks
    // (scb, w/ cuts_vol=1). Blocks with cuts_vol>1 will have one rank owning
    // a subdomain of those blocks, and that rank will own nothing else. Yet a
    // rank may handle one or more scbs.

    assert(cuts_vol_tot > comm_size); // Can't have more ranks than decompositions.

    int num_scb = 0; // Number of single-cut blocks (scb).
    int scb[GKYL_MAX_BLOCKS]; // Block ID of single-cut blocks (scb)
    for (int bidx=0; bidx<mba->num_blocks; bidx++) {
      int cuts_vol = 1;
      for (int d=0; d<mba->cdim; d++) cuts_vol *= inp->blocks[bidx]->cuts[d];

      cuts_vol_per_block[bidx] = cuts_vol;

      if (cuts_vol == 1) {
        scb[num_scb] = bidx;
        num_scb++;
      }
    }

    // Additional blocks that need to be assigned to a rank owning another block.
    int extra_blocks = cuts_vol_tot - comm_size;
    // Number of ranks owning single-cut blocks (scbrank).
    int num_scbrank = num_scb - extra_blocks;
    // Distribute scb amongst scbranks: 
    int scbrank_num_blocks[num_scbrank];
    int base = num_scb/num_scbrank;
    int rem = num_scb - num_scbrank * base;
    for (int i=0; i<num_scbrank; i++)
      scbrank_num_blocks[i] = base;
    for (int i=0; i<rem; i++)
      scbrank_num_blocks[i]++;

    // List of block IDs owned by each scbrank.
    int scbrank_blocks[num_scbrank * (base+1)];
    for (int i=0; i<num_scbrank; i++) {
      int scb_count = 0;
      for (int j=0; j<scbrank_num_blocks[i]; j++) {
        scbrank_blocks[i*(base+1)+j] = scb[scb_count];
        scb_count++;
      }
    }

    // List of ranks in each block.
    int curr_rank_to_assign = 0, curr_scbrank = 0, scbrank_count = 0;
    for (int bidx=0; bidx<mba->num_blocks; bidx++) {
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
          int scb_idx = scbrank_blocks[curr_scbrank*(base+1)];
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
    mba->num_blocks_local = 0;
    int my_block_idxs[GKYL_MAX_BLOCKS];
    for (int bidx=0; bidx<mba->num_blocks; bidx++) {
      for (int i=0; i<cuts_vol_per_block[bidx]; i++) {
        if (comm_rank == ranks_per_block[bidx*cuts_vol_max+i]) {
          my_block_idxs[mba->num_blocks_local] = bidx;
          mba->num_blocks_local++;
        }
      }
    }
    mba->block_idxs = gkyl_malloc(mba->num_blocks_local*sizeof(int));
    for (int i=0; i<mba->num_blocks_local; i++)
      mba->block_idxs[i] = my_block_idxs[i];

  }

  // Store the ranks_per_block for later use.
  mba->ranks_per_block = gkyl_malloc(cuts_vol_tot * sizeof(int));
  mba->cuts_vol_cum_per_block = gkyl_malloc(mba->num_blocks * sizeof(int));
  int cuts_vol_cum = 0;
  for (int bidx=0; bidx<mba->num_blocks; bidx++) {
    mba->cuts_vol_cum_per_block[bidx] = cuts_vol_cum;
    for (int i=0; i<cuts_vol_per_block[bidx]; i++) {
      mba->ranks_per_block[cuts_vol_cum+i] = ranks_per_block[bidx*cuts_vol_max+i];
    }
    cuts_vol_cum += cuts_vol_per_block[bidx];
  }

  // Allocate memory for all the block decompositions because the cross-block
  // field object needs to know the range of every block.
  mba->decomp_intrab = gkyl_malloc(mba->num_blocks * sizeof(struct gkyl_rect_decomp *));

  for (int bc=0; bc<mba->num_blocks; bc++) {
    struct gkyl_gk *blinp = inp->blocks[bc];

    // Create intra block decompositions.
    struct gkyl_range global_range_conf;
    gkyl_create_global_range(mba->cdim, blinp->cells, &global_range_conf);
    mba->decomp_intrab[bc] = gkyl_rect_decomp_new_from_cuts(mba->cdim, blinp->cuts, &global_range_conf);
  }
  gkyl_free(cuts_vol_per_block);
  gkyl_free(ranks_per_block);


  // Only allocate memory for local blocks and their intrablock communicators.
  mba->blocks = gkyl_malloc(mba->num_blocks_local * sizeof(struct gkyl_gyrokinetic_app));
  mba->comm_intrab = gkyl_malloc(mba->num_blocks_local * sizeof(struct gkyl_comm *));

  for (int bc=0; bc<mba->num_blocks_local; bc++) {
    int bidx = mba->block_idxs[bc];
    struct gkyl_gk *blinp = inp->blocks[bidx];

    // Create intra blok communicators.
    int comm_color = bidx;
    struct gkyl_comm *parent_comm = bc == 0? mba->comm_multib : mba->comm_intrab[0];
    mba->comm_intrab[bc] = gkyl_comm_split_comm(parent_comm, comm_color, mba->decomp_intrab[bidx]);

    int comm_intrab_rank;
    gkyl_comm_get_rank(mba->comm_intrab[bc], &comm_intrab_rank);

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
    blinp->low_inp.local_range = mba->decomp_intrab[bidx]->ranges[comm_intrab_rank],
    blinp->low_inp.comm = mba->comm_intrab[bc],

    // Create a new app for each block.
    mba->blocks[bc] = gkyl_gyrokinetic_app_new(blinp);
  }

  // Create a new field_multib object.
  mba->field_multib = gk_field_multib_new(inp, mba);

  return mba;
}

int
gkyl_gyrokinetic_multib_num_ranks_per_block(gkyl_gyrokinetic_multib_app *mba, int bidx)
{
  return mba->decomp_intrab[bidx]->ndecomp;
}

int
gkyl_gyrokinetic_multib_ranks_per_block(gkyl_gyrokinetic_multib_app *mba, int bidx, int *ranks)
{
  int cuts_vol = gkyl_gyrokinetic_multib_num_ranks_per_block(mba, bidx);
  int off = mba->cuts_vol_cum_per_block[bidx];
  for (int i=0; i<cuts_vol; i++)
    ranks[i] = mba->ranks_per_block[off+i];

  return cuts_vol; 
}

void
gkyl_gyrokinetic_multib_app_apply_ic(gkyl_gyrokinetic_multib_app* mba, double t0)
{
  for (int bidx=0; bidx<mba->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_apply_ic(mba->blocks[bidx], t0);
  }
  int fidx=0;
  gk_field_multib_rhs(mba, mba->field_multib, fidx);
}

static void
calc_field(gkyl_gyrokinetic_multib_app* mba, double tcurr, int fidx)
{
  // Compute fields.
  gk_field_multib_rhs(mba, mba->field_multib, fidx);
}

static void
calc_field_and_apply_bc(gkyl_gyrokinetic_multib_app* mba, double tcurr, int fidx)
{
  // Compute fields and apply BCs.

  // Compute the field.
  calc_field(mba, tcurr, fidx);

  // Apply boundary conditions.

  // NYI.

  // We need to 
  //   a) apply the BCs on physical boundaries.
  //   b) sync within each block
  //   c) sync between blocks.
  const int nblocks = mba->num_blocks_local;
  for (int bidx=0; bidx<nblocks; bidx++) {
    struct gkyl_gyrokinetic_app *app = mba->blocks[bidx];

    struct gkyl_array *f_charged[app->num_species];
    struct gkyl_array *f_neut[app->num_neut_species];
    for (int i=0; i<app->num_species; ++i) {
      f_charged[i] = app->species[i].distfs[fidx];
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      if (!app->neut_species[i].info.is_static) {
        f_neut[i] = app->neut_species[i].distfs[fidx];
      }
    }

    gkyl_gyrokinetic_apply_bc(app, f_charged, f_neut);
  }
}

static void
rk3_fill_fin_fout(const struct gkyl_gyrokinetic_app *app, int fidx_in, int fidx_out,
  const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[],
  struct gkyl_array *fout[], struct gkyl_array *fout_neut[])
{
  // Fill arrays of distributions functions with the appropriate distributions
  // specified by the fidx_in and fidx_out indices.

  for (int i=0; i<app->num_species; ++i) {
    fin[i] = app->species[i].distfs[fidx_in];
    fout[i] = app->species[i].distfs[fidx_out];
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    fin_neut[i] = app->neut_species[i].distfs[0];
    if (!app->neut_species[i].info.is_static) {
      fin_neut[i] = app->neut_species[i].distfs[fidx_in];
      fout_neut[i] = app->neut_species[i].distfs[fidx_out];
    }
  }
}

// Take time-step using the RK3 method. Also sets the status object
// which has the actual and suggested dts used. These can be different
// from the actual time-step.
static struct gkyl_update_status
rk3(gkyl_gyrokinetic_multib_app* mba, double dt0)
{
  const int nblocks = mba->num_blocks_local;

  struct gkyl_update_status st[nblocks];
  for (int bidx=0; bidx<nblocks; bidx++)
    st[bidx].success = true;

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = mba->tcurr, dt = dt0;
  int fidx_in = 0, fidx_out = 1;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
      {
        // Indices in gk_species' distfs of input and output distributions.
        const int fidx_in = 0, fidx_out = 1;

        // Compute df/dt = ...
        for (int bidx=0; bidx<nblocks; bidx++) {
          struct gkyl_gyrokinetic_app *app = mba->blocks[bidx];

          const struct gkyl_array *fin[app->num_species];
          const struct gkyl_array *fin_neut[app->num_neut_species];
          struct gkyl_array *fout[app->num_species];
          struct gkyl_array *fout_neut[app->num_neut_species];
          rk3_fill_fin_fout(app, fidx_in, fidx_out, fin, fin_neut, fout, fout_neut);

          gkyl_gyrokinetic_dfdt(app, tcurr, dt, fin, fout, fin_neut, fout_neut, &st[bidx]);
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

          struct gkyl_gyrokinetic_app *app = mba->blocks[bidx];
          const struct gkyl_array *fin[app->num_species];
          const struct gkyl_array *fin_neut[app->num_neut_species];
          struct gkyl_array *fout[app->num_species];
          struct gkyl_array *fout_neut[app->num_neut_species];
          rk3_fill_fin_fout(app, fidx_in, fidx_out, fin, fin_neut, fout, fout_neut);

          for (int i=0; i<app->num_species; ++i) {
            gkyl_array_accumulate(gkyl_array_scale(fout[i], dt_actual), 1.0, fin[i]);
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            if (!app->neut_species[i].info.is_static) {
              gkyl_array_accumulate(gkyl_array_scale(fout_neut[i], dt_actual), 1.0, fin_neut[i]);
            }
          }
        }

        // Compute the fields and apply BCs.
        calc_field_and_apply_bc(mba, tcurr, fidx_out);

        dt = dt_actual;
        state = RK_STAGE_2;
        break;
      }
      case RK_STAGE_2:
      {
        // Indices in gk_species' distfs of input and output distributions.
        const int fidx_in = 1, fidx_out = 2;

        for (int bidx=0; bidx<nblocks; bidx++) {
          struct gkyl_gyrokinetic_app *app = mba->blocks[bidx];

          const struct gkyl_array *fin[app->num_species];
          const struct gkyl_array *fin_neut[app->num_neut_species];
          struct gkyl_array *fout[app->num_species];
          struct gkyl_array *fout_neut[app->num_neut_species];
          rk3_fill_fin_fout(app, fidx_in, fidx_out, fin, fin_neut, fout, fout_neut);

          gkyl_gyrokinetic_dfdt(app, tcurr+dt, dt, fin, fout, fin_neut, fout_neut, &st[bidx]);
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

          struct gkyl_gyrokinetic_app *app = mba->blocks[bidx];
          const struct gkyl_array *fin[app->num_species];
          const struct gkyl_array *fin_neut[app->num_neut_species];
          struct gkyl_array *fout[app->num_species];
          struct gkyl_array *fout_neut[app->num_neut_species];
          rk3_fill_fin_fout(app, fidx_in, fidx_out, fin, fin_neut, fout, fout_neut);

          for (int i=0; i<app->num_species; ++i) {
            gkyl_array_accumulate(gkyl_array_scale(fout[i], dt_actual), 1.0, fin[i]);
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            if (!app->neut_species[i].info.is_static) {
              gkyl_array_accumulate(gkyl_array_scale(fout_neut[i], dt_actual), 1.0, fin_neut[i]);
            }
          }
        }

        if (dt_actual < dt) {

          // Recalculate the field.
          calc_field(mba, tcurr, 0);

          // collect stats
          double dt_rel_diff = (dt-dt_actual)/dt_actual;
          mba->stat.stage_2_dt_diff[0] = fmin(mba->stat.stage_2_dt_diff[0],
            dt_rel_diff);
          mba->stat.stage_2_dt_diff[1] = fmax(mba->stat.stage_2_dt_diff[1],
            dt_rel_diff);
          mba->stat.nstage_2_fail += 1;

          dt = dt_actual;
          state = RK_STAGE_1; // restart from stage 1

        } 
        else {
          for (int bidx=0; bidx<nblocks; bidx++) {
            struct gkyl_gyrokinetic_app *app = mba->blocks[bidx];

            for (int i=0; i<app->num_species; ++i) {
              array_combine(app->species[i].distfs[1],
                3.0/4.0, app->species[i].distfs[0], 1.0/4.0, app->species[i].distfs[2], &app->species[i].local_ext);
            }
            for (int i=0; i<app->num_neut_species; ++i) {
              if (!app->neut_species[i].info.is_static) {
                array_combine(app->neut_species[i].distfs[1],
                  3.0/4.0, app->neut_species[i].distfs[0], 1.0/4.0, app->neut_species[i].distfs[2], &app->neut_species[i].local_ext);
              }
            }
          }

          // Compute the fields and apply BCs.
          calc_field_and_apply_bc(mba, tcurr, 1);

          state = RK_STAGE_3;
        }
        break;
      }
      case RK_STAGE_3:
      {
        // Indices in gk_species' distfs of input and output distributions.
        const int fidx_in = 1, fidx_out = 2;

        for (int bidx=0; bidx<nblocks; bidx++) {
          struct gkyl_gyrokinetic_app *app = mba->blocks[bidx];

          const struct gkyl_array *fin[app->num_species];
          const struct gkyl_array *fin_neut[app->num_neut_species];
          struct gkyl_array *fout[app->num_species];
          struct gkyl_array *fout_neut[app->num_neut_species];
          rk3_fill_fin_fout(app, fidx_in, fidx_out, fin, fin_neut, fout, fout_neut);

          gkyl_gyrokinetic_dfdt(app, tcurr+dt/2, dt, fin, fout, fin_neut, fout_neut, &st[bidx]);
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

          struct gkyl_gyrokinetic_app *app = mba->blocks[bidx];
          const struct gkyl_array *fin[app->num_species];
          const struct gkyl_array *fin_neut[app->num_neut_species];
          struct gkyl_array *fout[app->num_species];
          struct gkyl_array *fout_neut[app->num_neut_species];
          rk3_fill_fin_fout(app, fidx_in, fidx_out, fin, fin_neut, fout, fout_neut);

          for (int i=0; i<app->num_species; ++i) {
            gkyl_array_accumulate(gkyl_array_scale(fout[i], dt_actual), 1.0, fin[i]);
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            if (!app->neut_species[i].info.is_static) {
              gkyl_array_accumulate(gkyl_array_scale(fout_neut[i], dt_actual), 1.0, fin_neut[i]);
            }
          }
        }

        if (dt_actual < dt) {
          // Recalculate the field.
          calc_field(mba, tcurr, 0);

          // collect stats
          double dt_rel_diff = (dt-dt_actual)/dt_actual;
          mba->stat.stage_3_dt_diff[0] = fmin(mba->stat.stage_3_dt_diff[0],
            dt_rel_diff);
          mba->stat.stage_3_dt_diff[1] = fmax(mba->stat.stage_3_dt_diff[1],
            dt_rel_diff);
          mba->stat.nstage_3_fail += 1;

          dt = dt_actual;
          state = RK_STAGE_1; // restart from stage 1

          mba->stat.nstage_2_fail += 1;
        }
        else {
          for (int bidx=0; bidx<nblocks; bidx++) {
            struct gkyl_gyrokinetic_app *app = mba->blocks[bidx];

            for (int i=0; i<app->num_species; ++i) {
              array_combine(app->species[i].distfs[1],
                1.0/3.0, app->species[i].distfs[0], 2.0/3.0, app->species[i].distfs[2], &app->species[i].local_ext);
              gkyl_array_copy_range(app->species[i].distfs[0], app->species[i].distfs[1], &app->species[i].local_ext);
            }
            for (int i=0; i<app->num_neut_species; ++i) {
              if (!app->neut_species[i].info.is_static) {
                array_combine(app->neut_species[i].distfs[1],
                  1.0/3.0, app->neut_species[i].distfs[0], 2.0/3.0, app->neut_species[i].distfs[2], &app->neut_species[i].local_ext);
                gkyl_array_copy_range(app->neut_species[i].distfs[0], app->neut_species[i].distfs[1], &app->neut_species[i].local_ext);
              }
            }
          }

          // Compute the fields and apply BCs
          calc_field_and_apply_bc(mba, tcurr, 0);

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
gkyl_gyrokinetic_multib_update(gkyl_gyrokinetic_multib_app* mba, double dt)
{
  mba->stat.nup += 1;
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_update_status status[GKYL_MAX_BLOCKS];
  status[0] = rk3(mba, dt);
  mba->tcurr += status[0].dt_actual;

  mba->stat.total_tm += gkyl_time_diff_now_sec(wst);
  // Check for any CUDA errors during time step
  if (mba->use_gpu)
    checkCuda(cudaGetLastError());
//  struct gkyl_update_status status[GKYL_MAX_BLOCKS];
//  for (int bidx=0; bidx<mba->num_blocks_local; bidx++) {
//    status[bidx] = gkyl_gyrokinetic_update(mba->blocks[bidx], dt);
//  }
  return status[0];
}

void
gkyl_gyrokinetic_multib_app_calc_mom(gkyl_gyrokinetic_multib_app* mba)
{
  for (int bidx=0; bidx<mba->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_calc_mom(mba->blocks[bidx]);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_integrated_mom(gkyl_gyrokinetic_multib_app* mba, double tm)
{
  for (int bidx=0; bidx<mba->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_calc_integrated_mom(mba->blocks[bidx], tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_mom(gkyl_gyrokinetic_multib_app* mba, double tm, int frame)
{
  for (int bidx=0; bidx<mba->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_write_mom(mba->blocks[bidx], tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_integrated_mom(gkyl_gyrokinetic_multib_app* mba)
{
  for (int bidx=0; bidx<mba->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_write_integrated_mom(mba->blocks[bidx]);
  }
}

void
gkyl_gyrokinetic_multib_app_write(gkyl_gyrokinetic_multib_app* mba, double tm, int frame)
{
  for (int bidx=0; bidx<mba->num_blocks_local; bidx++) {
    gkyl_gyrokinetic_app_write(mba->blocks[bidx], tm, frame);
  }
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_read_from_frame(gkyl_gyrokinetic_multib_app *app, int frame)
{
  // NYI
  struct gkyl_app_restart_status stat;
  return stat;
}

void
gkyl_gyrokinetic_multib_app_stat_write(gkyl_gyrokinetic_multib_app* mba)
{
  // NYI
}

struct gkyl_gyrokinetic_stat
gkyl_gyrokinetic_multib_app_stat(gkyl_gyrokinetic_multib_app* mba)
{
  for (int bidx=0; bidx<mba->num_blocks_local; bidx++) {
    gk_species_tm(mba->blocks[bidx]);
    gk_species_coll_tm(mba->blocks[bidx]);
  }
  return mba->blocks[0]->stat;
}

// private function to handle variable argument list for printing
static void
v_gk_app_cout(const gkyl_gyrokinetic_multib_app* mba, FILE *fp, const char *fmt, va_list argp)
{
  int rank, r = 0;
  gkyl_comm_get_rank(mba->comm_multib, &rank);
  if ((rank == 0) && fp)
    vfprintf(fp, fmt, argp);
}

void
gkyl_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* mba, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_gk_app_cout(mba, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_gyrokinetic_multib_app_release(gkyl_gyrokinetic_multib_app* mba)
{
  for (int i=0; i<mba->num_blocks_local; i++) {
    gkyl_gyrokinetic_app_release(mba->blocks[i]);
  }
  gkyl_free(mba->blocks);
  gkyl_free(mba->block_idxs);
  
  gkyl_free(mba->ranks_per_block);
  gkyl_free(mba->cuts_vol_cum_per_block);

  // Release decomp and comm.
  for (int i=0; i<mba->num_blocks; i++) {
    gkyl_rect_decomp_release(mba->decomp_intrab[i]);
  }
  for (int i=0; i<mba->num_blocks_local; i++) {
    gkyl_comm_release(mba->comm_intrab[i]);
  }
  gk_field_multib_release(mba, mba->field_multib);
  gkyl_free(mba->decomp_intrab);
  gkyl_free(mba->comm_intrab);

  gkyl_comm_release(mba->comm_multib);

  gkyl_block_topo_release(mba->btopo);
  gkyl_free(mba);
}
