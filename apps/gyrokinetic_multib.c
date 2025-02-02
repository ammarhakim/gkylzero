#include <gkyl_array_rio_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>

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

// compute total and maximum number of cuts
static void
calc_tot_and_max_cuts(const struct gkyl_block_geom *block_geom, int tot_max[2])
{
  int ndim = gkyl_block_geom_ndim(block_geom);
  int num_blocks = gkyl_block_geom_num_blocks(block_geom);

  int max_cuts = 0, tot_cuts = 0;
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *bgi = gkyl_block_geom_get_block(block_geom, i);
    int ncuts = calc_cuts(ndim, bgi->cuts);
    max_cuts = ncuts > max_cuts ? ncuts : max_cuts;
    tot_cuts += ncuts;
  }
  tot_max[0] = tot_cuts;
  tot_max[1] = max_cuts;
}

// construct the mpack meta-data for multi-block data files
static struct gkyl_msgpack_data *
gyrokinetic_multib_meta(struct gyrokinetic_multib_output_meta meta)
{
  struct gkyl_msgpack_data *mt = gkyl_malloc(sizeof *mt);

  mt->meta_sz = 0;
  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mt->meta, &mt->meta_sz);

  // add some data to mpack
  mpack_build_map(&writer);
  
  mpack_write_cstr(&writer, "time");
  mpack_write_double(&writer, meta.stime);

  mpack_write_cstr(&writer, "frame");
  mpack_write_i64(&writer, meta.frame);

  mpack_write_cstr(&writer, "topo_file");
  mpack_write_cstr(&writer, meta.topo_file_name);

  mpack_write_cstr(&writer, "app_name");
  mpack_write_cstr(&writer, meta.app_name);

  mpack_write_cstr(&writer, "Git_commit_hash");
  mpack_write_cstr(&writer, GIT_COMMIT_ID);

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    free(mt->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mt);
    mt = 0;
  }

  return mt;  
}

// write out multi-block data files
static int
gyrokinetic_multib_data_write(const char *fname, struct gyrokinetic_multib_output_meta meta)
{
  enum gkyl_array_rio_status status = GKYL_ARRAY_RIO_FOPEN_FAILED;
  FILE *fp = 0;
  int err;
  with_file (fp, fname, "w") {
    struct gkyl_msgpack_data *amet = gyrokinetic_multib_meta(meta);
    if (amet) {
      status = gkyl_header_meta_write_fp( &(struct gkyl_array_header_info) {
          .file_type = gkyl_file_type_int[GKYL_MULTI_BLOCK_DATA_FILE],
          .meta_size = amet->meta_sz,
          .meta = amet->meta
        },
        fp
      );
      MPACK_FREE(amet->meta);
      gkyl_free(amet);
    }
    else {
      status = GKYL_ARRAY_RIO_META_FAILED;
    }
  }
  return status;
}

// construct single-block App geometry for given block ID
static struct gkyl_gyrokinetic_app *
singleb_app_new_geom(const struct gkyl_gyrokinetic_multib *mbinp, int bid,
  const struct gkyl_gyrokinetic_multib_app *mbapp)
{
  // For kinetic simulations, block dimension defined configuration-space dimensionality.
  int cdim = gkyl_block_geom_ndim(mbapp->block_geom);
  int num_blocks = gkyl_block_geom_num_blocks(mbapp->block_geom);

  const struct gkyl_block_geom_info *bgi =
    gkyl_block_geom_get_block(mbapp->block_geom, bid);

  // construct top-level single-block input struct
  struct gkyl_gk app_inp = { };

  strcpy(app_inp.name, mbinp->name);
  if (num_blocks > 1) {
    cstr app_name = cstr_from_fmt("%s_b%d", mbinp->name, bid);
    strcpy(app_inp.name, app_name.str);
    cstr_drop(&app_name);
  }

  // Set the configuration-space extents, cells, and geometry.
  app_inp.cdim = cdim;
  for (int i=0; i<cdim; ++i) {
    app_inp.lower[i] = bgi->lower[i];
    app_inp.upper[i] = bgi->upper[i];
    app_inp.cells[i] = bgi->cells[i];
  }
  app_inp.geometry = bgi->geometry;
  int vdim = app_inp.vdim = mbinp->vdim;
  int num_species = app_inp.num_species = mbinp->num_species;
  int num_neut_species = app_inp.num_neut_species = mbinp->num_neut_species; 

  app_inp.poly_order = mbinp->poly_order;
  app_inp.basis_type = mbinp->basis_type;
  app_inp.cfl_frac = mbinp->cfl_frac; 

  struct gkyl_comm *comm = mbapp->block_comms[bid];

  struct gkyl_app_parallelism_inp parallel_inp = {};
  parallel_inp.use_gpu = mbinp->use_gpu;
  for (int d=0; d<cdim; ++d) parallel_inp.cuts[d] = bgi->cuts[d];
  parallel_inp.comm = comm;
  // Copy parallelism input into app input.
  memcpy(&app_inp.parallelism, &parallel_inp, sizeof(struct gkyl_app_parallelism_inp));

  return gkyl_gyrokinetic_app_new_geom(&app_inp);
}

// construct single-block App solver for given block ID
static void
singleb_app_new_solver(const struct gkyl_gyrokinetic_multib *mbinp, int bid,
  const struct gkyl_gyrokinetic_multib_app *mbapp, struct gkyl_gyrokinetic_app *app)
{
  // For kinetic simulations, block dimension defined configuration-space dimensionality.
  int cdim = gkyl_block_geom_ndim(mbapp->block_geom);
  int num_blocks = gkyl_block_geom_num_blocks(mbapp->block_geom);

  const struct gkyl_block_geom_info *bgi =
    gkyl_block_geom_get_block(mbapp->block_geom, bid);

  // construct top-level single-block input struct
  struct gkyl_gk app_inp = { };

  // Set the configuration-space extents, cells.
  app_inp.cdim = cdim;
  for (int i=0; i<cdim; ++i) {
    app_inp.lower[i] = bgi->lower[i];
    app_inp.upper[i] = bgi->upper[i];
    app_inp.cells[i] = bgi->cells[i];
  }

  int vdim = app_inp.vdim = mbinp->vdim;
  int num_species = app_inp.num_species = mbinp->num_species;
  int num_neut_species = app_inp.num_neut_species = mbinp->num_neut_species; 

  app_inp.poly_order = mbinp->poly_order;
  app_inp.basis_type = mbinp->basis_type;
  app_inp.cfl_frac = mbinp->cfl_frac; 

  app_inp.enforce_positivity = mbinp->enforce_positivity;

  for (int i=0; i<num_species; ++i) {
    const struct gkyl_gyrokinetic_multib_species *sp = &mbinp->species[i];
    
    struct gkyl_gyrokinetic_species species_inp = { };
    strcpy(species_inp.name, sp->name);

    species_inp.charge = sp->charge;
    species_inp.mass = sp->mass;
    species_inp.gkmodel_id = sp->gkmodel_id;
    species_inp.no_by = sp->no_by;

    // Velocity-space information
    for (int v=0; v<vdim; ++v) {
      species_inp.lower[v] = sp->lower[v];
      species_inp.upper[v] = sp->upper[v];
      species_inp.cells[v] = sp->cells[v];      
    }
    species_inp.mapc2p = sp->mapc2p;

    // Species physics modules
    species_inp.collisions = sp->collisions;
    species_inp.diffusion = sp->diffusion;
    species_inp.radiation = sp->radiation;
    species_inp.react = sp->react;
    species_inp.react_neut = sp->react_neut; 

    // Species diagnostics
    species_inp.num_diag_moments = sp->num_diag_moments;
    for (int n=0; n<species_inp.num_diag_moments; ++n) {
      strcpy(species_inp.diag_moments[n], sp->diag_moments[n]);
    }

    // choose proper block-specific species input
    const struct gkyl_gyrokinetic_multib_species_pb *sp_pb = &sp->blocks[0];
    if (!sp->duplicate_across_blocks) {
      for (int i=0; i<num_blocks; ++i) {
        if (bid == sp->blocks[i].block_id) {
          sp_pb = &sp->blocks[i];
          break;
        }
      }
    }
    species_inp.projection = sp_pb->projection;
    species_inp.source = sp_pb->source;
    species_inp.polarization_density = sp_pb->polarization_density;  

    // by default, skip BCs altogether
    species_inp.bcx.lower.type = GKYL_SPECIES_SKIP;
    species_inp.bcx.upper.type = GKYL_SPECIES_SKIP;
    species_inp.bcy.lower.type = GKYL_SPECIES_SKIP;
    species_inp.bcy.upper.type = GKYL_SPECIES_SKIP;
    species_inp.bcz.lower.type = GKYL_SPECIES_SKIP;
    species_inp.bcz.upper.type = GKYL_SPECIES_SKIP;

    // set species physical BCs: we need to search through the list of
    // physical BCs and set the appropriate input to single-block
    // species inp
    for (int i=0; i<sp->num_physical_bcs; ++i) {
      if (bid == sp->bcs[i].bidx) {

        int e = sp->bcs[i].edge;
        if (sp->bcs[i].dir == 0) {
          if (e == 0) {
            species_inp.bcx.lower.type = sp->bcs[i].bc_type;
          }
          else {
            species_inp.bcx.upper.type = sp->bcs[i].bc_type;
          }
        }
        else if (sp->bcs[i].dir == 1) {
          if (e == 0) {
            species_inp.bcy.lower.type = sp->bcs[i].bc_type;
          }
          else {
            species_inp.bcy.upper.type = sp->bcs[i].bc_type;
          }
        }
        else {
          if (e == 0) {
            species_inp.bcz.lower.type = sp->bcs[i].bc_type;
          }
          else {
            species_inp.bcz.upper.type = sp->bcs[i].bc_type;
          }
        }
      }
    }

    // copy species input into app input
    memcpy(&app_inp.species[i], &species_inp, sizeof(struct gkyl_gyrokinetic_species));
  }

  for (int i=0; i<num_neut_species; ++i) {
    const struct gkyl_gyrokinetic_multib_neut_species *nsp = &mbinp->neut_species[i];
    
    struct gkyl_gyrokinetic_neut_species neut_species_inp = { };
    strcpy(neut_species_inp.name, nsp->name);

    neut_species_inp.mass = nsp->mass; 
    neut_species_inp.is_static = nsp->is_static; 

    // Velocity space information (neutrals are 3V)
    for (int v=0; v<3; ++v) {
      neut_species_inp.lower[v] = nsp->lower[v];
      neut_species_inp.upper[v] = nsp->upper[v];
      neut_species_inp.cells[v] = nsp->cells[v];
    }
    neut_species_inp.mapc2p = nsp->mapc2p;

    // Neutral species physics modules
    neut_species_inp.react_neut = nsp->react_neut;

    // Neutral species diagnostics
    neut_species_inp.num_diag_moments = nsp->num_diag_moments;
    for (int n=0; n<neut_species_inp.num_diag_moments; ++n) {
      strcpy(neut_species_inp.diag_moments[n], nsp->diag_moments[n]);
    }

    // choose proper block-specific species input
    const struct gkyl_gyrokinetic_multib_neut_species_pb *nsp_pb = &nsp->blocks[0];
    if (!nsp->duplicate_across_blocks) {
      for (int i=0; i<num_blocks; ++i) {
        if (bid == nsp->blocks[i].block_id) {
          nsp_pb = &nsp->blocks[i];
          break;
        }
      }
    }
    neut_species_inp.projection = nsp_pb->projection;
    neut_species_inp.source = nsp_pb->source;

    // by default, skip BCs altogether
    neut_species_inp.bcx.lower.type = GKYL_SPECIES_SKIP;
    neut_species_inp.bcx.upper.type = GKYL_SPECIES_SKIP;
    neut_species_inp.bcy.lower.type = GKYL_SPECIES_SKIP;
    neut_species_inp.bcy.upper.type = GKYL_SPECIES_SKIP;
    neut_species_inp.bcz.lower.type = GKYL_SPECIES_SKIP;
    neut_species_inp.bcz.upper.type = GKYL_SPECIES_SKIP;

    // set species physical BCs: we need to search through the list of
    // physical BCs and set the appropriate input to single-block
    // species inp
    for (int i=0; i<nsp->num_physical_bcs; ++i) {
      if (bid == nsp->bcs[i].bidx) {

        int e = nsp->bcs[i].edge;
        if (nsp->bcs[i].dir == 0) {
          if (e == 0) {
            neut_species_inp.bcx.lower.type = nsp->bcs[i].bc_type;
          }
          else {
            neut_species_inp.bcx.upper.type = nsp->bcs[i].bc_type;
          }
        }
        else if (nsp->bcs[i].dir == 1) {
          if (e == 0) {
            neut_species_inp.bcy.lower.type = nsp->bcs[i].bc_type;
          }
          else {
            neut_species_inp.bcy.upper.type = nsp->bcs[i].bc_type;
          }
        }
        else {
          if (e == 0) {
            neut_species_inp.bcz.lower.type = nsp->bcs[i].bc_type;
          }
          else {
            neut_species_inp.bcz.upper.type = nsp->bcs[i].bc_type;
          }
        }
      }
    }

    // copy neutral species input into app input
    memcpy(&app_inp.neut_species[i], &neut_species_inp, sizeof(struct gkyl_gyrokinetic_neut_species));
  } 

  // Initialize field
  const struct gkyl_gyrokinetic_multib_field *fld = &mbinp->field;
  struct gkyl_gyrokinetic_field field_inp = { };
  field_inp.gkfield_id = fld->gkfield_id;
  field_inp.kperpSq = fld->kperpSq; 
  field_inp.xLCFS = fld->xLCFS; 

  // Adiabatic electron inputs
  field_inp.electron_mass = fld->electron_mass;
  field_inp.electron_charge = fld->electron_charge;
  field_inp.electron_density = fld->electron_density; 
  field_inp.electron_temp = fld->electron_temp; 

  // BCs
  // MF 2024/10/20: hardcode the BC value here because input file infra doesn't
  // support passing a value yet.
  for (int d=0; d<cdim-1; d++) {
    for (int k=0; k<fld->num_physical_bcs; k++) { 
      if (fld->bcs[k].edge == GKYL_LOWER_EDGE)
        field_inp.poisson_bcs.lo_type[d] = fld->bcs[k].bc_type;
      if (fld->bcs[k].edge == GKYL_UPPER_EDGE)
        field_inp.poisson_bcs.up_type[d] = fld->bcs[k].bc_type;
    }
  }

  const struct gkyl_gyrokinetic_multib_field_pb *fld_pb = &fld->blocks[0];
  // choose proper block-specific field input
  if (!fld->duplicate_across_blocks) {
    for (int i=0; i<num_blocks; ++i) {
      if (bid == fld->blocks[i].block_id) {
        const struct gkyl_gyrokinetic_multib_field_pb *fld_pb = &fld->blocks[i];
        field_inp.fem_parbc = fld_pb->fem_parbc;
        break;
      }
    }
  }
  else {
    field_inp.fem_parbc = fld_pb->fem_parbc;
  }

  if (!fld->duplicate_across_blocks) {
    for (int i=0; i<num_blocks; ++i) {
      if (bid == fld->blocks[i].block_id) {
        fld_pb = &fld->blocks[i];
        break;
      }
    }
  }

  field_inp.polarization_bmag = fld_pb->polarization_bmag ? fld_pb->polarization_bmag : mbapp->bmag_ref;

  field_inp.phi_wall_lo_ctx = fld_pb->phi_wall_lo_ctx; 
  field_inp.phi_wall_lo = fld_pb->phi_wall_lo; 
  field_inp.phi_wall_lo_evolve = fld_pb->phi_wall_lo_evolve; 

  field_inp.phi_wall_up_ctx = fld_pb->phi_wall_up_ctx; 
  field_inp.phi_wall_up = fld_pb->phi_wall_up; 
  field_inp.phi_wall_up_evolve = fld_pb->phi_wall_up_evolve;   

  // Copy field input into app input
  memcpy(&app_inp.field, &field_inp, sizeof(struct gkyl_gyrokinetic_field));  
  
  gkyl_gyrokinetic_app_new_solver(&app_inp, app);
}

gkyl_gyrokinetic_multib_app* gkyl_gyrokinetic_multib_app_new(const struct gkyl_gyrokinetic_multib *mbinp)
{
  int my_rank, num_ranks;
  gkyl_comm_get_rank(mbinp->comm, &my_rank);
  gkyl_comm_get_size(mbinp->comm, &num_ranks);

  int tot_max[2];
  calc_tot_and_max_cuts(mbinp->block_geom, tot_max);
  if ((num_ranks > tot_max[0]) || (num_ranks < tot_max[1])) {
    fprintf(stderr, "\nSpecified %d total cuts but provided %d processes, \
and the maximum number of cuts in a block is %d\n\n", tot_max[0], num_ranks, tot_max[1]);
    return 0;
  }

  struct gkyl_gyrokinetic_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));

  strcpy(mbapp->name, mbinp->name);
  mbapp->comm = gkyl_comm_acquire(mbinp->comm);  
  mbapp->use_gpu = mbinp->use_gpu;
  
  mbapp->block_geom = gkyl_block_geom_acquire(mbinp->block_geom);
  mbapp->block_topo = gkyl_block_geom_topo(mbinp->block_geom);
  
  int cdim = gkyl_block_geom_ndim(mbapp->block_geom);
  int num_blocks = gkyl_block_geom_num_blocks(mbapp->block_geom);

  // construct round-robin decomposition
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *bgi = gkyl_block_geom_get_block(mbapp->block_geom, i);
    branks[i] = calc_cuts(cdim, bgi->cuts);
  }
  mbapp->round_robin = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int num_local_blocks = 0;
  mbapp->local_blocks = gkyl_malloc(sizeof(int[num_blocks]));

  int lidx = 0;
  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // this is larger than needed

  mbapp->decomp = gkyl_malloc(num_blocks*sizeof(struct gkyl_rect_decomp*));
  
  // construct list of block communicators: there are as many
  // communicators as blocks. Not all communicators are valid on each
  // rank. The total number of valid communicators is
  // num_local_blocks.
  mbapp->block_comms = gkyl_malloc(num_blocks*sizeof(struct gkyl_comm *));
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(mbapp->round_robin, i, rank_list);

    bool is_my_rank_in_decomp = has_int(branks[i], my_rank, rank_list);

    if (is_my_rank_in_decomp) {
      mbapp->local_blocks[lidx++] = i;
      num_local_blocks += 1;      
    }

    const struct gkyl_block_geom_info *bgi = gkyl_block_geom_get_block(mbapp->block_geom, i);
    struct gkyl_range block_global_range;
    gkyl_create_global_range(cdim, bgi->cells, &block_global_range);

    mbapp->decomp[i] = gkyl_rect_decomp_new_from_cuts(
      cdim, bgi->cuts, &block_global_range);

    bool status;
    mbapp->block_comms[i] = gkyl_comm_create_comm_from_ranks(mbinp->comm,
      branks[i], rank_list, mbapp->decomp[i], &status);
  }

  mbapp->num_local_blocks = num_local_blocks;  

  printf("Rank %d handles %d Apps\n", my_rank, num_local_blocks);
  for (int i=0; i<num_local_blocks; ++i)
    printf("  Rank %d handles block %d\n", my_rank, mbapp->local_blocks[i]);

  mbapp->num_species = 0;
  mbapp->num_neut_species = 0;
  mbapp->update_field = 0;
  mbapp->singleb_apps = 0;

  if (num_local_blocks > 0) {
    mbapp->num_species = mbinp->num_species;
    mbapp->num_neut_species = mbinp->num_neut_species;
    mbapp->update_field = !mbinp->skip_field; // note inversion of truth value (default: update field)


    mbapp->singleb_apps = gkyl_malloc(num_local_blocks*sizeof(struct gkyl_gyrokinetic_app*));
  }

  for (int i=0; i<mbinp->num_species; ++i)
    strcpy(mbapp->species_name[i], mbinp->species[i].name);

  for (int i=0; i<mbinp->num_neut_species; ++i)
    strcpy(mbapp->neut_species_name[i], mbinp->neut_species[i].name);  

  // Create single-block grids and geometries.
  for (int i=0; i<num_local_blocks; ++i)
    mbapp->singleb_apps[i] = singleb_app_new_geom(mbinp, mbapp->local_blocks[i], mbapp);

  // Set bmag_ref
  double bmag_min_local = DBL_MAX;
  double bmag_min_global;
  for (int i=0; i<num_local_blocks; ++i) {
    double bmag_min = gkyl_gk_geometry_reduce_bmag(mbapp->singleb_apps[i]->gk_geom, GKYL_MIN);
    bmag_min_local = GKYL_MIN2(bmag_min_local, bmag_min);
  }
  gkyl_comm_allreduce_host(mbapp->comm, GKYL_DOUBLE, GKYL_MIN, 1, &bmag_min_local, &bmag_min_global);

  double bmag_max_local = -DBL_MAX;
  double bmag_max_global;
  for (int i=0; i<num_local_blocks; ++i) {
    double bmag_max = gkyl_gk_geometry_reduce_bmag(mbapp->singleb_apps[i]->gk_geom, GKYL_MAX);
    bmag_max_local = GKYL_MAX2(bmag_max_local, bmag_max);
  }
  gkyl_comm_allreduce_host(mbapp->comm, GKYL_DOUBLE, GKYL_MAX, 1, &bmag_max_local, &bmag_max_global);

  mbapp->bmag_ref = (bmag_max_global + bmag_min_global)/2.0;

  // Create the rest of the single-block solvers.
  for (int i=0; i<num_local_blocks; ++i)
    singleb_app_new_solver(mbinp, mbapp->local_blocks[i], mbapp, mbapp->singleb_apps[i]);

  // Create the MB field app.
  mbapp->field = gk_multib_field_new(mbinp, mbapp);

  // Create connections needed for conf-space syncs.
  int ghost[] = { 1, 1, 1 };
  mbapp->mbcc_sync_conf = gkyl_malloc(sizeof(struct gkyl_mbcc_sr));
  mbapp->mbcc_sync_conf->send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  mbapp->mbcc_sync_conf->recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));

  for (int bI=0; bI<num_local_blocks; ++bI) {
    int bid = mbapp->local_blocks[bI];

    gkyl_rrobin_decomp_getranks(mbapp->round_robin, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    mbapp->mbcc_sync_conf->recv[bI] = gkyl_multib_comm_conn_new_recv(bid, brank, 
      ghost, &mbapp->block_topo->conn[bid], mbapp->decomp);
    mbapp->mbcc_sync_conf->send[bI] = gkyl_multib_comm_conn_new_send(bid, brank,
      ghost, &mbapp->block_topo->conn[bid], mbapp->decomp);

    struct gkyl_multib_comm_conn *mbcc_s = mbapp->mbcc_sync_conf->send[bI], *mbcc_r = mbapp->mbcc_sync_conf->recv[bI];
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];

    for (int ns=0; ns<mbcc_s->num_comm_conn; ++ns) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      struct gkyl_comm_conn *ccs = &mbcc_s->comm_conn[ns];
      int rankIdx = ccs->rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, ccs->block_id, rank_list);
      ccs->rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&ccs->range, &sbapp->local_ext, ccs->range.lower, ccs->range.upper);
    }
    for (int nr=0; nr<mbcc_r->num_comm_conn; ++nr) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      struct gkyl_comm_conn *ccr = &mbcc_r->comm_conn[nr];
      int rankIdx = ccr->rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, ccr->block_id, rank_list);
      ccr->rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&ccr->range, &sbapp->local_ext, ccr->range.lower, ccr->range.upper);
    }

    // Sort connections according to rank and block ID (needed by NCCL).
    gkyl_multib_comm_conn_sort(mbcc_r);
    gkyl_multib_comm_conn_sort(mbcc_s);
  }

  // Create connections needed for syncing charged species phase-space quantities.
  mbapp->mbcc_sync_charged = gkyl_malloc(mbapp->num_species * sizeof(struct gkyl_mbcc_sr));
  for (int i=0; i<mbinp->num_species; ++i) {
    mbapp->mbcc_sync_charged[i].send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
    mbapp->mbcc_sync_charged[i].recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  }
  for (int bI=0; bI<num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    struct gkyl_multib_comm_conn *mbcc_cs = mbapp->mbcc_sync_conf->send[bI], *mbcc_cr = mbapp->mbcc_sync_conf->recv[bI];
    for (int i=0; i<mbinp->num_species; ++i) {
      mbapp->mbcc_sync_charged[i].send[bI] = gkyl_multib_comm_conn_new(mbcc_cs->num_comm_conn, mbcc_cs->comm_conn);
      mbapp->mbcc_sync_charged[i].recv[bI] = gkyl_multib_comm_conn_new(mbcc_cr->num_comm_conn, mbcc_cr->comm_conn);
      struct gkyl_multib_comm_conn *mbcc_ps = mbapp->mbcc_sync_charged[i].send[bI],
                                   *mbcc_pr = mbapp->mbcc_sync_charged[i].recv[bI];
      // Extend ranges to include velocity space, and make them a sub range.
      struct gk_species *gks = &sbapp->species[i];
      for (int ns=0; ns<mbcc_cs->num_comm_conn; ++ns) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cs->comm_conn[ns];
        struct gkyl_comm_conn *ccs_phase = &mbcc_ps->comm_conn[ns];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gks->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gks->local_ext, phase_r.lower, phase_r.upper);
      }
      for (int nr=0; nr<mbcc_cr->num_comm_conn; ++nr) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cr->comm_conn[nr];
        struct gkyl_comm_conn *ccs_phase = &mbcc_pr->comm_conn[nr];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gks->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gks->local_ext, phase_r.lower, phase_r.upper);
      }
    }
  }

  // Create connections needed for syncing neutral species phase-space quantities.
  mbapp->mbcc_sync_neut = gkyl_malloc(mbapp->num_neut_species * sizeof(struct gkyl_mbcc_sr));
  for (int i=0; i<mbinp->num_neut_species; ++i) {
    mbapp->mbcc_sync_neut[i].send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
    mbapp->mbcc_sync_neut[i].recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  }
  for (int bI=0; bI<num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    struct gkyl_multib_comm_conn *mbcc_cs = mbapp->mbcc_sync_conf->send[bI], *mbcc_cr = mbapp->mbcc_sync_conf->recv[bI];
    for (int i=0; i<mbinp->num_neut_species; ++i) {
      mbapp->mbcc_sync_neut[i].send[bI] = gkyl_multib_comm_conn_new(mbcc_cs->num_comm_conn, mbcc_cs->comm_conn);
      mbapp->mbcc_sync_neut[i].recv[bI] = gkyl_multib_comm_conn_new(mbcc_cr->num_comm_conn, mbcc_cr->comm_conn);
      struct gkyl_multib_comm_conn *mbcc_ps = mbapp->mbcc_sync_neut[i].send[bI],
                                   *mbcc_pr = mbapp->mbcc_sync_neut[i].recv[bI];
      // Extend ranges to include velocity space, and make them a sub range.
      struct gk_neut_species *gkns = &sbapp->neut_species[i];
      for (int ns=0; ns<mbcc_cs->num_comm_conn; ++ns) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cs->comm_conn[ns];
        struct gkyl_comm_conn *ccs_phase = &mbcc_ps->comm_conn[ns];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gkns->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gkns->local_ext, phase_r.lower, phase_r.upper);
      }
      for (int nr=0; nr<mbcc_cr->num_comm_conn; ++nr) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cr->comm_conn[nr];
        struct gkyl_comm_conn *ccs_phase = &mbcc_pr->comm_conn[nr];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gkns->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gkns->local_ext, phase_r.lower, phase_r.upper);
      }
    }
  }

  // Sync the conf-space Jacobian needed for syncing quantities that include a
  // jacobgeo factor in them.
  struct gkyl_array *jacs[mbapp->num_local_blocks];
  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
    jacs[b] = sbapp->gk_geom->jacobgeo;
  }
  // Sync across blocks.
  gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
    mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs, jacs);

  // Allocate updaters to rescale jac*f by the jacobians in the skin/ghost cells.
  while (true) {
    struct gkyl_gyrokinetic_app *sbapp0 = mbapp->singleb_apps[0];
    for (int d=0; d<cdim; d++) {
      for (int e=0; e<2; e++) {
        mbapp->jf_rescale_charged[d*2+e] = gkyl_rescale_ghost_jacf_new(d,e,&sbapp0->confBasis,
          &sbapp0->basis, mbapp->use_gpu);
        if ((mbapp->num_neut_species > 0) && (!sbapp0->neut_species[0].info.is_static)) {
          mbapp->jf_rescale_neut[d*2+e] = gkyl_rescale_ghost_jacf_new(d,e,&sbapp0->confBasis,
            &sbapp0->neut_basis, mbapp->use_gpu);
        }
      }
    }
    break;
  }

  mbapp->stat = (struct gkyl_gyrokinetic_stat) {};

  gkyl_free(rank_list);
  gkyl_free(branks);

  return mbapp;
}

void
gyrokinetic_multib_calc_field(struct gkyl_gyrokinetic_multib_app* app, double tcurr, const struct gkyl_array *fin[])
{
  // Compute fields.
  if (app->update_field) {
    // Solve the field equation.
    gk_multib_field_rhs(app, app->field, fin);
  }
}

static void
gyrokinetic_multib_apply_bc(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array *distf_neut[])
{
  // Apply boundary conditions in each block (including intrablock sync).
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int li_charged = b * app->num_species;
    int li_neut = b * app->num_neut_species;
    for (int i=0; i<app->num_species; ++i) {
      gk_species_apply_bc(sbapp, &sbapp->species[i], distf[li_charged+i]);
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      gk_neut_species_apply_bc(sbapp, &sbapp->neut_species[i], distf_neut[li_neut+i]);
    }
  }

  // Sync blocks.
  int cdim = gkyl_block_geom_ndim(app->block_geom);
  for (int i=0; i<app->num_species; ++i) {
    // Sync charged species.
    struct gkyl_array *fs[app->num_local_blocks];
    for (int b=0; b<app->num_local_blocks; ++b) {
      int li_charged = b * app->num_species;
      fs[b] = distf[li_charged+i];
    }
    gkyl_multib_comm_conn_array_transfer(app->comm, app->num_local_blocks, app->local_blocks,
      app->mbcc_sync_charged[i].send, app->mbcc_sync_charged[i].recv, fs, fs);

    // Divide and multiply by the appropriate jacobians.
    for (int b=0; b<app->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
      int bid = app->local_blocks[b];
      int li_charged = b * app->num_species;

      struct gk_species *gks = &sbapp->species[i];

      for (int dir=0; dir<cdim; ++dir) {
        for (int e=0; e<2; ++e) {
          if (app->block_topo->conn[bid].connections[dir][e].edge != GKYL_PHYSICAL) {
            gkyl_rescale_ghost_jacf_advance(app->jf_rescale_charged[dir*2+e],
              e==0? &sbapp->global_lower_skin[dir] : &sbapp->global_upper_skin[dir],
              e==0? &sbapp->global_lower_ghost[dir] : &sbapp->global_upper_ghost[dir],
              e==0? &gks->global_lower_ghost[dir] : &gks->global_upper_ghost[dir],
              sbapp->gk_geom->jacobgeo, distf[li_charged+i]);
          }
        }
      }
    }
  }

  struct gkyl_gyrokinetic_app *sbapp0 = app->singleb_apps[0];
  for (int i=0; i<app->num_neut_species; ++i) {
    // Sync neutral species.
    if (!sbapp0->neut_species[i].info.is_static) {
      struct gkyl_array *fs[app->num_local_blocks];
      for (int b=0; b<app->num_local_blocks; ++b) {
        int li_neut = b * app->num_neut_species;
        fs[b] = distf_neut[li_neut+i];
      }
      gkyl_multib_comm_conn_array_transfer(app->comm, app->num_local_blocks, app->local_blocks,
        app->mbcc_sync_neut[i].send, app->mbcc_sync_neut[i].recv, fs, fs);
      // Divide and multiply by the appropriate jacobians.
      for (int b=0; b<app->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
        int bid = app->local_blocks[b];
        int li_neut = b * app->num_neut_species;

        for (int dir=0; dir<cdim; ++dir) {
          for (int e=0; e<2; ++e) {
            if (app->block_topo->conn[bid].connections[dir][e].edge != GKYL_PHYSICAL) {
              struct gk_neut_species *gkns = &sbapp->neut_species[i];
              gkyl_rescale_ghost_jacf_advance(app->jf_rescale_neut[dir*2+e],
                e==0? &sbapp->global_lower_skin[dir] : &sbapp->global_upper_skin[dir],
                e==0? &sbapp->global_lower_ghost[dir] : &sbapp->global_upper_ghost[dir],
                e==0? &gkns->global_lower_ghost[dir] : &gkns->global_upper_ghost[dir],
                sbapp->gk_geom->jacobgeo, distf_neut[li_neut+i]);
            }
          }
        }
      }
    }
  }
}

void
gyrokinetic_multib_calc_field_and_apply_bc(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array *distf_neut[])
{
  // Compute fields and apply BCs.

  // Compute the field.
  // MF 2024/09/27/: Need the cast here for consistency. Fixing
  // this may require removing 'const' from a lot of places.
  gyrokinetic_multib_calc_field(app, tcurr, (const struct gkyl_array **) distf);

  // Apply boundary conditions.
  gyrokinetic_multib_apply_bc(app, tcurr, distf, distf_neut);
}

void
gkyl_gyrokinetic_multib_app_apply_ic(gkyl_gyrokinetic_multib_app* app, double t0)
{
  app->tcurr = t0;
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_apply_ic_species(app, i, t0);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_apply_ic_neut_species(app, i, t0);
  }  

  // Compute the fields and apply BCs.
  struct gkyl_array *distf[app->num_species * app->num_local_blocks];
  struct gkyl_array *distf_neut[app->num_neut_species * app->num_local_blocks];
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int li_charged = b * app->num_species;
    int li_neut = b * app->num_neut_species;
    for (int i=0; i<app->num_species; ++i) {
      distf[li_charged+i] = sbapp->species[i].f;
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      distf_neut[li_neut+i] = sbapp->neut_species[i].f;
    }
  }
  gyrokinetic_multib_calc_field_and_apply_bc(app, t0, distf, distf_neut);
}

void
gkyl_gyrokinetic_multib_app_apply_ic_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  app->tcurr = t0;
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_apply_ic_species(app->singleb_apps[b], sidx, t0);
  }
  gkyl_comm_barrier(app->comm);
}

void
gkyl_gyrokinetic_multib_app_apply_ic_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  app->tcurr = t0;
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_apply_ic_neut_species(app->singleb_apps[b], sidx, t0);
  }
  gkyl_comm_barrier(app->comm);
}


struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_file_field(gkyl_gyrokinetic_multib_app *app, const char *fname)
{
 // TO DO
  return (struct gkyl_app_restart_status) { };

}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_read_from_frame(gkyl_gyrokinetic_multib_app *app, int frame)
{
  struct gkyl_app_restart_status rstat;
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    for (int i=0; i<app->num_neut_species; i++) {
      int neut_frame = frame;
      if (sbapp->neut_species[i].info.is_static) {
        neut_frame = 0;
      }
      rstat = gkyl_gyrokinetic_app_from_frame_neut_species(sbapp, i, neut_frame);
    }
    for (int i=0; i<app->num_species; i++) {
      rstat = gkyl_gyrokinetic_app_from_frame_species(sbapp, i, frame);
    }
  }
  
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    // Compute the fields and apply BCs.
    struct gkyl_array *distf[app->num_species * app->num_local_blocks];
    struct gkyl_array *distf_neut[app->num_neut_species * app->num_local_blocks];
    for (int b=0; b<app->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
      int li_charged = b * app->num_species;
      int li_neut = b * app->num_neut_species;
      for (int i=0; i<app->num_species; ++i) {
        distf[li_charged+i] = sbapp->species[i].f;
      }
      for (int i=0; i<app->num_neut_species; ++i) {
        distf_neut[li_neut+i] = sbapp->neut_species[i].f;
      }
    }
//    if (app->update_field && app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
//      for (int i=0; i<app->num_species; ++i) {
//        struct gk_species *s = &app->species[i];
//
//        // Compute advection speeds so we can compute the initial boundary flux.
//        gkyl_dg_calc_gyrokinetic_vars_alpha_surf(s->calc_gk_vars, 
//          &app->local, &s->local, &s->local_ext, app->field->phi_smooth,
//          s->alpha_surf, s->sgn_alpha_surf, s->const_sgn_alpha);
//
//        // Compute and store (in the ghost cell of of out) the boundary fluxes.
//        // NOTE: this overwrites ghost cells that may be used for sourcing.
//        gk_species_bflux_rhs(app, s, &s->bflux, distf[i], distf[i]);
//      }
//    }
    gyrokinetic_multib_calc_field_and_apply_bc(app, rstat.stime, distf, distf_neut);
  }

  assert(app->num_local_blocks == 1);
  struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[0];
  sbapp->field->is_first_energy_write_call = false; // Append to existing diagnostic.

  return rstat;
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_field(gkyl_gyrokinetic_multib_app *app, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

// private function to handle variable argument list for printing
static void
v_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp)
    vfprintf(fp, fmt, argp);
}

void
gkyl_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_gyrokinetic_multib_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_gyrokinetic_multib_app_write_topo(const gkyl_gyrokinetic_multib_app* app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (0 == rank) {
    cstr file_name = cstr_from_fmt("%s_btopo.gkyl", app->name);
    gkyl_block_topo_write(app->block_topo, file_name.str);
    cstr_drop(&file_name);
  }
}

void
gkyl_gyrokinetic_multib_app_write_field(gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_field(app->singleb_apps[b], tm, frame);
  }

// MF 2024/10/20: This stuff is corrupting the file.
//  if (app->update_field) {
//    int rank;
//    gkyl_comm_get_rank(app->comm, &rank);
//    if (0 == rank) {
//      cstr file_name = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "field", frame);
//      cstr topo_file_name = cstr_from_fmt("%s_btopo.gkyl", app->name);
//      
//      gyrokinetic_multib_data_write(file_name.str, (struct gyrokinetic_multib_output_meta) {
//          .frame = frame,
//          .stime = tm,
//          .topo_file_name = topo_file_name.str,
//          .app_name = app->name
//        }
//      );
//      
//      cstr_drop(&topo_file_name);
//      cstr_drop(&file_name);
//    }
//  }
//
//  gkyl_comm_barrier(app->comm);
}

void
gkyl_gyrokinetic_multib_app_calc_field_energy(gkyl_gyrokinetic_multib_app* app, double tm)
{
  if (app->update_field) {
    for (int b=0; b<app->num_local_blocks; ++b) {
      gkyl_gyrokinetic_app_calc_field_energy(app->singleb_apps[b], tm);
    }
  }
}

void
gkyl_gyrokinetic_multib_app_write_field_energy(gkyl_gyrokinetic_multib_app* app)
{
  if (app->update_field) {
    for (int b=0; b<app->num_local_blocks; ++b) {
      gkyl_gyrokinetic_app_write_field_energy(app->singleb_apps[b]);
    }
  }
}

//
// ............. Species outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species(app->singleb_apps[b], sidx, tm, frame);
  }

// MF 2024/10/20: This stuff is corrupting the file.
//  int rank;
//  gkyl_comm_get_rank(app->comm, &rank);
//  if (0 == rank) {
//    cstr file_name = cstr_from_fmt("%s-%s_%d.gkyl", app->name, app->species_name[sidx], frame);
//    cstr topo_file_name = cstr_from_fmt("%s_btopo.gkyl", app->name);
//      
//    gyrokinetic_multib_data_write(file_name.str, (struct gyrokinetic_multib_output_meta) {
//        .frame = frame,
//        .stime = tm,
//        .topo_file_name = topo_file_name.str,
//        .app_name = app->name
//      }
//    );
//    
//    cstr_drop(&topo_file_name);
//    cstr_drop(&file_name);
//  }
//
//  gkyl_comm_barrier(app->comm);
}

void
gkyl_gyrokinetic_multib_app_write_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species(app->singleb_apps[b], sidx, tm, frame);
  }

// MF 2024/10/20: This stuff is corrupting the file.
//  int rank;
//  gkyl_comm_get_rank(app->comm, &rank);
//  if (0 == rank) {
//    cstr file_name = cstr_from_fmt("%s-%s_%d.gkyl", app->name, app->neut_species_name[sidx], frame);
//    cstr topo_file_name = cstr_from_fmt("%s_btopo.gkyl", app->name);
//      
//    gyrokinetic_multib_data_write(file_name.str, (struct gyrokinetic_multib_output_meta) {
//        .frame = frame,
//        .stime = tm,
//        .topo_file_name = topo_file_name.str,
//        .app_name = app->name
//      }
//    );
//    
//    cstr_drop(&topo_file_name);
//    cstr_drop(&file_name);
//  }
//
//  gkyl_comm_barrier(app->comm);
}

void
gkyl_gyrokinetic_multib_app_write_species_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_neut_species_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_neut_species_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_integrated_mom(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_integrated_mom(app->singleb_apps[b], sidx);
  }
}

//
// ............. Source outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_source(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_source(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_source(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_source_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_source_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_source_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_source_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_source_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_neut_species_source_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_neut_species_source_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_source_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source_integrated_mom(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_source_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_source_integrated_mom(app->singleb_apps[b], sidx);
  }
}

//
// ............. Collision outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_lbo_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_lbo_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_max_corr_status(gkyl_gyrokinetic_multib_app* app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_max_corr_status(app->singleb_apps[b], sidx);
  }
}

//
// ............. Radiation outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_rad_drag(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_rad_drag(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_rad_emissivity(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_rad_emissivity(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_rad_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_rad_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_rad_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_rad_integrated_mom(app->singleb_apps[b], sidx);
  }
}

//
// ............. Functions that group several outputs for a single species ............... //
//
void
gkyl_gyrokinetic_multib_app_write_species_phase(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_species(app, sidx, tm, frame);

  gkyl_gyrokinetic_multib_app_write_species_source(app, sidx, tm, frame);

  gkyl_gyrokinetic_multib_app_write_species_rad_drag(app, sidx, tm, frame);
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_phase(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_neut_species(app, sidx, tm, frame);

  gkyl_gyrokinetic_multib_app_write_neut_species_source(app, sidx, tm, frame);
}

void
gkyl_gyrokinetic_multib_app_write_species_conf(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_conf(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_conf(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_conf(app->singleb_apps[b], sidx, tm, frame);
  }
}

//
// ............. Functions that group several species outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_mom(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_source_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_lbo_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_rad_emissivity(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_neut_species_source_mom(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_integrated_mom(gkyl_gyrokinetic_multib_app* app, double tm)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_calc_species_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_source_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_rad_integrated_mom(app, i, tm);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_calc_neut_species_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_neut_species_source_integrated_mom(app, i, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_integrated_mom(gkyl_gyrokinetic_multib_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_species_source_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_species_max_corr_status(app, i);
    gkyl_gyrokinetic_multib_app_write_species_rad_integrated_mom(app, i);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_neut_species_source_integrated_mom(app, i);
  }
}

void
gkyl_gyrokinetic_multib_app_write_conf(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_field(app, tm, frame);

  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_conf(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_conf(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_phase(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_phase(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_phase(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_phase(app, tm, frame);

  gkyl_gyrokinetic_multib_app_write_conf(app, tm, frame);
}

//
// ............. End of write functions ............... //
// 

struct gkyl_update_status
gkyl_gyrokinetic_multib_update(gkyl_gyrokinetic_multib_app* app, double dt)
{
  app->stat.nup += 1;
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_update_status status = gyrokinetic_multib_update_ssp_rk3(app, dt);
  app->tcurr += status.dt_actual;

  app->stat.total_tm += gkyl_time_diff_now_sec(wst);

  // Check for any CUDA errors during time step
  if (app->use_gpu)
    checkCuda(cudaGetLastError());

  return status;
}

struct gkyl_gyrokinetic_stat
gkyl_gyrokinetic_multib_app_stat(gkyl_gyrokinetic_multib_app* app)
{
  app->stat.species_rhs_tm = 0.0;
  app->stat.field_rhs_tm = 0.0;
  app->stat.species_coll_mom_tm = 0.0;
  app->stat.species_coll_tm = 0.0;
  app->stat.species_bc_tm = 0.0;
  app->stat.nspecies_omega_cfl = 0;
  app->stat.species_omega_cfl_tm = 0.0;
  app->stat.ndiag = 0;
  app->stat.diag_tm = 0.0;
  app->stat.nio = 0;
  app->stat.io_tm = 0.0;
  for (int i=0; i<app->num_species; ++i) {
    app->stat.species_lbo_coll_diff_tm[i] = 0.0;
    app->stat.species_lbo_coll_drag_tm[i] = 0.0;
  }

  for (int b=0; b<app->num_local_blocks; ++b) {
    // Add time spent on various operations for each local block.
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    struct gkyl_gyrokinetic_stat sb_stat = gkyl_gyrokinetic_app_stat(sbapp);

    app->stat.species_rhs_tm += sb_stat.species_rhs_tm;
    app->stat.field_rhs_tm += sb_stat.field_rhs_tm;
    app->stat.species_coll_mom_tm += sb_stat.species_coll_mom_tm;
    app->stat.species_coll_tm += sb_stat.species_coll_tm;
    app->stat.species_bc_tm += sb_stat.species_bc_tm;
    app->stat.nspecies_omega_cfl += sb_stat.nspecies_omega_cfl;
    app->stat.species_omega_cfl_tm += sb_stat.species_omega_cfl_tm;
    app->stat.ndiag += sb_stat.ndiag;
    app->stat.diag_tm += sb_stat.diag_tm;
    app->stat.nio += sb_stat.nio;
    app->stat.io_tm += sb_stat.io_tm;
    for (int i=0; i<app->num_species; ++i) {
      app->stat.species_lbo_coll_diff_tm[i] += sbapp->stat.species_lbo_coll_diff_tm[i];
      app->stat.species_lbo_coll_drag_tm[i] += sbapp->stat.species_lbo_coll_drag_tm[i];
    }
  }
  return app->stat;
}

void gkyl_gyrokinetic_multib_app_species_ktm_rhs(gkyl_gyrokinetic_multib_app* app, int update_vol_term)
{
  // TO DO
}

void
gkyl_gyrokinetic_multib_app_stat_write(gkyl_gyrokinetic_multib_app* app)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    gkyl_gyrokinetic_app_stat_write(sbapp);
  }
}

void gkyl_gyrokinetic_multib_app_release(gkyl_gyrokinetic_multib_app* mbapp)
{
  while (true) {
    struct gkyl_gyrokinetic_app *sbapp0 = mbapp->singleb_apps[0];
    for (int d=0; d<sbapp0->cdim; d++) {
      for (int e=0; e<2; e++) {
        gkyl_rescale_ghost_jacf_release(mbapp->jf_rescale_charged[d*2+e]);
        if ((mbapp->num_neut_species > 0) && (!sbapp0->neut_species[0].info.is_static)) {
          gkyl_rescale_ghost_jacf_release(mbapp->jf_rescale_neut[d*2+e]);
        }
      }
    }
    break;
  }

  for (int i=0; i<mbapp->num_neut_species; ++i) {
    for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_neut[i].send[bI]);
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_neut[i].recv[bI]);
    }
    gkyl_free(mbapp->mbcc_sync_neut[i].send);
    gkyl_free(mbapp->mbcc_sync_neut[i].recv);
  }
  gkyl_free(mbapp->mbcc_sync_neut);

  for (int i=0; i<mbapp->num_species; ++i) {
    for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_charged[i].send[bI]);
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_charged[i].recv[bI]);
    }
    gkyl_free(mbapp->mbcc_sync_charged[i].send);
    gkyl_free(mbapp->mbcc_sync_charged[i].recv);
  }
  gkyl_free(mbapp->mbcc_sync_charged);

  for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->send[bI]);
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->recv[bI]);
  }
  gkyl_free(mbapp->mbcc_sync_conf->send);
  gkyl_free(mbapp->mbcc_sync_conf->recv);
  gkyl_free(mbapp->mbcc_sync_conf);

  if (mbapp->singleb_apps) {
    for (int i=0; i<mbapp->num_local_blocks; ++i)
      gkyl_gyrokinetic_app_release(mbapp->singleb_apps[i]);
    gkyl_free(mbapp->singleb_apps);
  }  

  gk_multib_field_release(mbapp->field);

  int num_blocks = gkyl_block_geom_num_blocks(mbapp->block_geom);

  for (int i=0; i<num_blocks; ++i)
    gkyl_comm_release(mbapp->block_comms[i]);
  gkyl_free(mbapp->block_comms);

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(mbapp->decomp[i]);
  gkyl_free(mbapp->decomp);

  gkyl_free(mbapp->local_blocks);    

  gkyl_rrobin_decomp_release(mbapp->round_robin);
  
  gkyl_block_geom_release(mbapp->block_geom);
  gkyl_block_topo_release(mbapp->block_topo);
  
  gkyl_comm_release(mbapp->comm);

  gkyl_free(mbapp);
}
