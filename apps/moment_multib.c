#include <gkyl_moment_multib.h>
#include <gkyl_moment_multib_priv.h>
#include <gkyl_rrobin_decomp.h>

#include <mpack.h>

// compute total number of ranges specified by cuts
static inline int
calc_cuts(int ndim, const int *cuts)
{
  int tc = 1;
  for (int d=0; d<ndim; ++d) tc *= cuts[d];
  return tc;
}

// simple linear search to check if val occurs in lst
static bool
has_int(int n, int val, const int *lst)
{
  for (int i=0; i<n; ++i)
    if (val == lst[i])
      return true;
  return false;
}

// construct single-block App for given block ID
struct gkyl_moment_app *
singleb_app_new(const struct gkyl_moment_multib *mbinp, int bid,
  const struct gkyl_moment_multib_app *mbapp)
{
  int ndim = gkyl_block_geom_ndim(mbapp->block_geom);
  int num_blocks = gkyl_block_geom_num_blocks(mbapp->block_geom);

  const struct gkyl_block_geom_info *bgi =
    gkyl_block_geom_get_block(mbapp->block_geom, bid);

  // construct top-level single-block input struct
  struct gkyl_moment app_inp = { };

  strcpy(app_inp.name, mbinp->name);
  if (num_blocks > 1) {
    cstr app_name = cstr_from_fmt("%s_b%d", mbinp->name, bid);
    strcpy(app_inp.name, app_name.str);
    cstr_drop(&app_name);
  }

  app_inp.ndim = ndim;
  for (int i=0; i<ndim; ++i) {
    app_inp.lower[i] = bgi->lower[i];
    app_inp.upper[i] = bgi->upper[i];
    app_inp.cells[i] = bgi->cells[i];
  }

  app_inp.cfl_frac = mbinp->cfl_frac;

  app_inp.scheme_type = mbinp->scheme_type;
  app_inp.mp_recon = mbinp->mp_recon;
  app_inp.skip_mp_limiter = mbinp->skip_mp_limiter;
  app_inp.use_hybrid_flux_kep = mbinp->use_hybrid_flux_kep;

  app_inp.num_skip_dirs = mbinp->num_skip_dirs;
  for (int d=0; d<3; ++d)
    app_inp.skip_dirs[d] = mbinp->skip_dirs[d];

  int num_species = app_inp.num_species = mbinp->num_species;
  // construct each species input block
  for (int i=0; i<num_species; ++i) {
    const struct gkyl_moment_multib_species *sp = &mbinp->species[i];
    
    struct gkyl_moment_species species_inp = { };
    strcpy(species_inp.name, sp->name);

    species_inp.charge = sp->charge;
    species_inp.mass = sp->mass;
    
    species_inp.equation = sp->equation;
    species_inp.limiter = sp->limiter;
    species_inp.split_type = sp->split_type;

    species_inp.evolve = sp->evolve;
    species_inp.force_low_order_flux = sp->force_low_order_flux;

    // choose proper block-specific species input
    const struct gkyl_moment_multib_species_pb *sp_pb = &sp->blocks[0];
    if (!sp->duplicate_across_blocks) {
      for (int i=0; i<num_blocks; ++i)
        if (bid == sp->blocks[i].block_id) {
          sp_pb = &sp->blocks[i];
          break;
        }
    }

    species_inp.ctx = sp_pb->ctx;
    species_inp.init = sp_pb->init;

    species_inp.is_app_accel_static = sp_pb->is_app_accel_static;
    species_inp.app_accel_ctx = sp_pb->app_accel_ctx;
    species_inp.app_accel_func = sp_pb->app_accel_func;

    species_inp.nT_source_ctx = sp_pb->nT_source_ctx;
    species_inp.nT_source_func = sp_pb->nT_source_func;
    species_inp.nT_source_set_only_once = sp_pb->nT_source_set_only_once;

    // by default, skip BCs altogether
    for (int e=0; e<2; ++e) {
      species_inp.bcx[e] = GKYL_SPECIES_SKIP;
      species_inp.bcy[e] = GKYL_SPECIES_SKIP;
      species_inp.bcz[e] = GKYL_SPECIES_SKIP;
    }

    // set species physical BCs: we need to search through the list of
    // physical BCs and set the appropriate input to single-block
    // species inp
    for (int i=0; i<sp->num_physical_bcs; ++i) {
      if (bid == sp->bcs[i].bidx) {

        int e = sp->bcs[i].edge;
        if (sp->bcs[i].dir == 0)
          species_inp.bcx[e] = sp->bcs[i].bc_type;
        else if (sp->bcs[i].dir == 1)
          species_inp.bcy[e] = sp->bcs[i].bc_type;
        else
          species_inp.bcz[e] = sp->bcs[i].bc_type;
      }
    }

    // copy species input into app input
    memcpy(&app_inp.species[i], &species_inp, sizeof(struct gkyl_moment_species));
  }


  // construct field input
  if (mbinp->field.blocks) { // assumption is that blocks = 0 when there is no field
    struct gkyl_moment_field field_inp = { };
    const struct gkyl_moment_multib_field *fld = &mbinp->field;

    field_inp.epsilon0 = fld->epsilon0;
    field_inp.mu0 = fld->mu0;
    field_inp.elc_error_speed_fact = fld->elc_error_speed_fact;
    field_inp.mag_error_speed_fact = fld->mag_error_speed_fact;

    field_inp.limiter = fld->limiter;
    field_inp.evolve = fld->evolve;

    // choose proper block-specific field input
    const struct gkyl_moment_multib_field_pb *fld_pb = &fld->blocks[0];
    if (!fld->duplicate_across_blocks) {
      for (int i=0; i<num_blocks; ++i)
        if (bid == fld->blocks[i].block_id) {
          fld_pb = &fld->blocks[i];
          break;
        }
    }

    field_inp.ctx = fld_pb->ctx;
    field_inp.init = fld_pb->init;

    field_inp.app_current_ctx = fld_pb->app_current_ctx;
    field_inp.app_current_func = fld_pb->app_current_func;
    field_inp.t_ramp_curr = fld_pb->t_ramp_curr;

    field_inp.is_ext_em_static = fld_pb->is_ext_em_static;
    field_inp.ext_em_ctx = fld_pb->ext_em_ctx;
    field_inp.ext_em_func = fld_pb->ext_em_func;
    field_inp.t_ramp_ext_em = fld_pb->t_ramp_ext_em;

    field_inp.use_explicit_em_coupling = fld_pb->use_explicit_em_coupling;

    // by default, skip BCs altogether
    for (int e=0; e<2; ++e) {
      field_inp.bcx[e] = GKYL_FIELD_SKIP;
      field_inp.bcy[e] = GKYL_FIELD_SKIP;
      field_inp.bcz[e] = GKYL_FIELD_SKIP;
    }

    // set field physical BCs: we need to search through the list of
    // physical BCs and set the appropriate input to single-block
    // field inp
    for (int i=0; i<fld->num_physical_bcs; ++i) {
      if (bid == fld->bcs[i].bidx) {

        int e = fld->bcs[i].edge;
        if (fld->bcs[i].dir == 0)
          field_inp.bcx[e] = fld->bcs[i].bc_type;
        else if (fld->bcs[i].dir == 1)
          field_inp.bcy[e] = fld->bcs[i].bc_type;
        else
          field_inp.bcz[e] = fld->bcs[i].bc_type;
      }
    }    

    // copy species input into app input
    memcpy(&app_inp.field, &field_inp, sizeof(struct gkyl_moment_field));
  }

  // Set decomp information: we are recreating the decomp here but
  // this seems to be the cleanest way to do this without storing
  // decomp for each block on each rank.
  struct gkyl_range block_global_range;
  gkyl_create_global_range(ndim, bgi->cells, &block_global_range);

  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(
    ndim, bgi->cuts, &block_global_range);

  struct gkyl_comm *comm = mbapp->block_comms[bid];
  int local_rank;
  gkyl_comm_get_rank(comm, &local_rank);

  app_inp.has_low_inp = true;
  app_inp.low_inp = (struct gkyl_app_comm_low_inp) {
    .comm = comm,
    .local_range = decomp->ranges[local_rank]
  };

  gkyl_rect_decomp_release(decomp);  
  
  return gkyl_moment_app_new(&app_inp);
}

struct gkyl_moment_multib_app *
gkyl_moment_multib_app_new(const struct gkyl_moment_multib *mbinp)
{
  struct gkyl_moment_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));
  strcpy(mbapp->name, mbinp->name);
  mbapp->comm = gkyl_comm_acquire(mbinp->comm);

  mbapp->block_geom = gkyl_block_geom_acquire(mbinp->block_geom);
  int ndim = gkyl_block_geom_ndim(mbapp->block_geom);
  int num_blocks = gkyl_block_geom_num_blocks(mbapp->block_geom);

  int my_rank;
  gkyl_comm_get_rank(mbapp->comm, &my_rank);
  
  int num_ranks;
  gkyl_comm_get_size(mbapp->comm, &num_ranks);

  // construct round-robin decomposition
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *bgi = gkyl_block_geom_get_block(mbapp->block_geom, i);
    branks[i] = calc_cuts(ndim, bgi->cuts);
  }
  const struct gkyl_rrobin_decomp *rrd = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int num_local_blocks = 0;
  mbapp->local_blocks = gkyl_malloc(sizeof(int[num_blocks]));

  int lidx = 0;
  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // this is larger than needed
  
  // construct list of block communicators: there are as many
  // communictor as blocks. Not all communicators are valid on each
  // rank. The total number of valid communictors is num_local_blocks.
  mbapp->block_comms = gkyl_malloc(num_blocks*sizeof(struct gkyl_comm *));
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(rrd, i, rank_list);

    bool is_my_rank_in_decomp = has_int(branks[i], my_rank, rank_list);

    if (is_my_rank_in_decomp) {
      mbapp->local_blocks[lidx++] = i;
      num_local_blocks += 1;      
    }

    const struct gkyl_block_geom_info *bgi = gkyl_block_geom_get_block(mbapp->block_geom, i);
    struct gkyl_range block_global_range;
    gkyl_create_global_range(ndim, bgi->cells, &block_global_range);

    struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(
      ndim, bgi->cuts, &block_global_range);

    bool status;
    mbapp->block_comms[i] = gkyl_comm_create_comm_from_ranks(mbinp->comm,
      branks[i], rank_list, decomp, &status);

    gkyl_rect_decomp_release(decomp);
  }
  gkyl_free(rank_list);
  mbapp->num_local_blocks = num_local_blocks;  

  printf("Rank %d handles %d Apps\n", my_rank, num_local_blocks);
  for (int i=0; i<num_local_blocks; ++i)
    printf("  Rank %d handles block %d\n", my_rank, mbapp->local_blocks[i]);

  // create individual single-block Apps
  mbapp->apps = 0;
  if (num_local_blocks > 0)
    mbapp->apps = gkyl_malloc(num_local_blocks*sizeof(struct gkyl_moment_app*));

  for (int i=0; i<num_local_blocks; ++i)
    mbapp->apps[i] = singleb_app_new(mbinp, mbapp->local_blocks[i], mbapp);

  mbapp->stat = (struct gkyl_moment_stat) {
  };

  gkyl_free(branks);
  gkyl_rrobin_decomp_release(rrd);
  
  return mbapp;
}

double
gkyl_moment_multib_app_max_dt(gkyl_moment_multib_app *app)
{
  // TODO
  return 0;
}


void
gkyl_moment_multib_app_apply_ic(gkyl_moment_multib_app* app, double t0)
{
  // TODO
}

void
gkyl_moment_multib_app_apply_ic_field(gkyl_moment_multib_app* app, double t0)
{
  // TODO
}

void
gkyl_moment_multib_app_apply_ic_species(gkyl_moment_multib_app* app, int sidx, double t0)
{
  // TODO
}

struct gkyl_app_restart_status
gkyl_moment_multib_app_from_frame_field(gkyl_moment_multib_app *app,
  int frame)
{
  // TODO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_moment_multib_app_from_frame_species(gkyl_moment_multib_app *app,
  int sidx, int frame)
{
  // TODO
  return (struct gkyl_app_restart_status) { };  
}

void
gkyl_moment_multib_app_cout(const gkyl_moment_multib_app* app, FILE *fp, const char *fmt, ...)
{
  // TODO
}

void
gkyl_moment_multib_app_write(const gkyl_moment_multib_app* app, double tm, int frame)
{
  // TODO
}

void
gkyl_moment_multib_app_write_field(const gkyl_moment_multib_app *app, double tm, int frame)
{
  // TODO
}

void
gkyl_moment_multib_app_write_species(const gkyl_moment_multib_app* app, int sidx, double tm, int frame)
{
  // TODO
}

void
gkyl_moment_multib_app_write_field_energy(gkyl_moment_multib_app *app)
{
    // TODO
}

void
gkyl_moment_multib_app_write_integrated_mom(gkyl_moment_multib_app *app)
{
  // TODO
}

void
gkyl_moment_multib_app_stat_write(const gkyl_moment_multib_app *app)
{
  // TODO
}

struct gkyl_update_status
gkyl_moment_multib_update(gkyl_moment_multib_app *app, double dt)
{
  // TODO
  return (struct gkyl_update_status) { };
}

void
gkyl_moment_multib_app_calc_field_energy(gkyl_moment_multib_app *app, double tm)
{
  // TODO  
}

void
gkyl_moment_multib_app_get_field_energy(gkyl_moment_multib_app *app, double *vals)
{
  // TODO
}

void
gkyl_moment_multib_app_calc_integrated_mom(gkyl_moment_multib_app *app, double tm)
{
  // TODO
}

struct gkyl_moment_stat
gkyl_moment_multib_app_stat(gkyl_moment_multib_app *app)
{
  return app->stat;
}

void
gkyl_moment_multib_app_release(gkyl_moment_multib_app* mbapp)
{
  if (mbapp->apps) {
    for (int i=0; i<mbapp->num_local_blocks; ++i)
      gkyl_moment_app_release(mbapp->apps[i]);
    gkyl_free(mbapp->apps);
  }  

  int num_blocks = gkyl_block_geom_num_blocks(mbapp->block_geom);

  for (int i=0; i<num_blocks; ++i)
    gkyl_comm_release(mbapp->block_comms[i]);
  gkyl_free(mbapp->block_comms);
  gkyl_comm_release(mbapp->comm);  
  
  gkyl_block_geom_release(mbapp->block_geom);
  gkyl_free(mbapp->local_blocks);    
  
  gkyl_free(mbapp);
}
