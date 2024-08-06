#include <gkyl_array_rio_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>
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
static struct gkyl_array_meta *
gyrokinetic_multib_meta(struct gyrokinetic_multib_output_meta meta)
{
  struct gkyl_array_meta *mt = gkyl_malloc(sizeof *mt);

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
    struct gkyl_array_meta *amet = gyrokinetic_multib_meta(meta);
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

// construct single-block App for given block ID
static struct gkyl_gyrokinetic_app *
singleb_app_new(const struct gkyl_gyrokinetic_multib *mbinp, int bid,
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
  app_inp.use_gpu = mbinp->use_gpu;
  app_inp.cfl_frac = mbinp->cfl_frac;  

  for (int i=0; i<num_species; ++i) {
    const struct gkyl_gyrokinetic_multib_species *sp = &mbinp->species[i];
    
    struct gkyl_gyrokinetic_species species_inp = { };
    strcpy(species_inp.name, sp->name);

    species_inp.charge = sp->charge;
    species_inp.mass = sp->mass;
    species_inp.gkmodel_id = sp->gkmodel_id;
    species_inp.no_by = sp->no_by;
    species_inp.enforce_positivity = sp->enforce_positivity;

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

  // choose proper block-specific field input
  const struct gkyl_gyrokinetic_multib_field_pb *fld_pb = &fld->blocks[0];
  if (!fld->duplicate_across_blocks) {
    for (int i=0; i<num_blocks; ++i) {
      if (bid == fld->blocks[i].block_id) {
        fld_pb = &fld->blocks[i];
        break;
      }
    }
  }

  field_inp.polarization_bmag = fld_pb->polarization_bmag;

  field_inp.phi_wall_lo_ctx = fld_pb->phi_wall_lo_ctx; 
  field_inp.phi_wall_lo = fld_pb->phi_wall_lo; 
  field_inp.phi_wall_lo_evolve = fld_pb->phi_wall_lo_evolve; 

  field_inp.phi_wall_up_ctx = fld_pb->phi_wall_up_ctx; 
  field_inp.phi_wall_up = fld_pb->phi_wall_up; 
  field_inp.phi_wall_up_evolve = fld_pb->phi_wall_up_evolve;   

  // copy field input into app input
  memcpy(&app_inp.field, &field_inp, sizeof(struct gkyl_gyrokinetic_field));  

  struct gkyl_comm *comm = mbapp->block_comms[bid];
  int local_rank;
  gkyl_comm_get_rank(comm, &local_rank);

  app_inp.has_low_inp = true;
  app_inp.low_inp = (struct gkyl_app_comm_low_inp) {
    .comm = comm,
    .local_range = mbapp->decomp[bid]->ranges[local_rank]
  };
  
  return gkyl_gyrokinetic_app_new(&app_inp);    
}

gkyl_gyrokinetic_multib_app* gkyl_gyrokinetic_multib_app_new(const struct gkyl_gyrokinetic_multib *mbinp)
{
  int my_rank;
  gkyl_comm_get_rank(mbinp->comm, &my_rank);
  int num_ranks;
  gkyl_comm_get_size(mbinp->comm, &num_ranks);

  int tot_max[2];
  calc_tot_and_max_cuts(mbinp->block_geom, tot_max);
  if ((num_ranks > tot_max[0]) || (num_ranks < tot_max[1]))
    return 0;

  struct gkyl_gyrokinetic_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));
  strcpy(mbapp->name, mbinp->name);
  mbapp->comm = gkyl_comm_acquire(mbinp->comm);  
  
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
  gkyl_free(rank_list);
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

  for (int i=0; i<num_local_blocks; ++i)
    mbapp->singleb_apps[i] = singleb_app_new(mbinp, mbapp->local_blocks[i], mbapp);

  mbapp->stat = (struct gkyl_gyrokinetic_stat) {
  };

  gkyl_free(branks);
  
  return mbapp;
}


void gkyl_gyrokinetic_multib_app_apply_ic(gkyl_gyrokinetic_multib_app* app, double t0)
{
  app->tcurr = t0;
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_apply_ic_species(app, i, t0);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_apply_ic_neut_species(app, i, t0);
  }  
}

void gkyl_gyrokinetic_multib_app_apply_ic_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  app->tcurr = t0;
  for (int i=0; i<app->num_local_blocks; ++i) {
    gkyl_gyrokinetic_app_apply_ic_species(app->singleb_apps[i], sidx, t0);
  }
  gkyl_comm_barrier(app->comm);
}

void gkyl_gyrokinetic_multib_app_apply_ic_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  app->tcurr = t0;
  for (int i=0; i<app->num_local_blocks; ++i) {
    gkyl_gyrokinetic_app_apply_ic_neut_species(app->singleb_apps[i], sidx, t0);
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
  // TO DO
  return (struct gkyl_app_restart_status) { };
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
v_gyrokinetic_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, va_list argp)
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
  v_gyrokinetic_app_cout(app, fp, fmt, argp);
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

void gkyl_gyrokinetic_multib_app_calc_mom(gkyl_gyrokinetic_multib_app *app)
{
  for (int i=0; i<app->num_local_blocks; ++i) {
    gkyl_gyrokinetic_app_calc_mom(app->singleb_apps[i]);
  }
}

void gkyl_gyrokinetic_multib_app_calc_integrated_mom(gkyl_gyrokinetic_multib_app* app, double tm)
{
  for (int i=0; i<app->num_local_blocks; ++i) {
    gkyl_gyrokinetic_app_calc_integrated_mom(app->singleb_apps[i], tm);
  }
  // TO DO: REDUCE ACROSS BLOCKS
}

void gkyl_gyrokinetic_multib_app_calc_integrated_neut_mom(gkyl_gyrokinetic_multib_app* app, double tm)
{
  for (int i=0; i<app->num_local_blocks; ++i) {
    gkyl_gyrokinetic_app_calc_integrated_neut_mom(app->singleb_apps[i], tm);
  }  
  // TO DO: REDUCE ACROSS BLOCKS
}

void gkyl_gyrokinetic_multib_app_calc_field_energy(gkyl_gyrokinetic_multib_app* app, double tm)
{
  for (int i=0; i<app->num_local_blocks; ++i) {
    gkyl_gyrokinetic_app_calc_field_energy(app->singleb_apps[i], tm);
  }    
  // TO DO: REDUCE ACROSS BLOCKS
}

void
gkyl_gyrokinetic_multib_app_write(const gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species(app, i, tm, frame);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species(app, i, tm, frame);
  }  
}

void
gkyl_gyrokinetic_multib_app_write_field(const gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{
  for (int i=0; i<app->num_local_blocks; ++i) {
    gkyl_gyrokinetic_app_write_field(app->singleb_apps[i], tm, frame);
  }

  if (app->update_field) {
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (0 == rank) {
      cstr file_name = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "field", frame);
      cstr topo_file_name = cstr_from_fmt("%s_btopo.gkyl", app->name);
      
      gyrokinetic_multib_data_write(file_name.str, (struct gyrokinetic_multib_output_meta) {
          .frame = frame,
          .stime = tm,
          .topo_file_name = topo_file_name.str,
          .app_name = app->name
        }
      );
      
      cstr_drop(&topo_file_name);
      cstr_drop(&file_name);
    }
  }

  gkyl_comm_barrier(app->comm);
}

void
gkyl_gyrokinetic_multib_app_write_species(const gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int i=0; i<app->num_local_blocks; ++i) {
    gkyl_gyrokinetic_app_write_species(app->singleb_apps[i], sidx, tm, frame);
  }

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (0 == rank) {
    cstr file_name = cstr_from_fmt("%s-%s_%d.gkyl", app->name, app->species_name[sidx], frame);
    cstr topo_file_name = cstr_from_fmt("%s_btopo.gkyl", app->name);
      
    gyrokinetic_multib_data_write(file_name.str, (struct gyrokinetic_multib_output_meta) {
        .frame = frame,
        .stime = tm,
        .topo_file_name = topo_file_name.str,
        .app_name = app->name
      }
    );
    
    cstr_drop(&topo_file_name);
    cstr_drop(&file_name);
  }

  gkyl_comm_barrier(app->comm);
}

void
gkyl_gyrokinetic_multib_app_write_neut_species(const gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int i=0; i<app->num_local_blocks; ++i) {
    gkyl_gyrokinetic_app_write_neut_species(app->singleb_apps[i], sidx, tm, frame);
  }

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (0 == rank) {
    cstr file_name = cstr_from_fmt("%s-%s_%d.gkyl", app->name, app->neut_species_name[sidx], frame);
    cstr topo_file_name = cstr_from_fmt("%s_btopo.gkyl", app->name);
      
    gyrokinetic_multib_data_write(file_name.str, (struct gyrokinetic_multib_output_meta) {
        .frame = frame,
        .stime = tm,
        .topo_file_name = topo_file_name.str,
        .app_name = app->name
      }
    );
    
    cstr_drop(&topo_file_name);
    cstr_drop(&file_name);
  }

  gkyl_comm_barrier(app->comm);
}

void gkyl_gyrokinetic_multib_app_write_source_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_source_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_coll_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_rad_drag(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_rad_emissivity(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_rad_integrated_moms(gkyl_gyrokinetic_multib_app *app, int sidx, double tm)
{
  // TO DO
}


void gkyl_gyrokinetic_multib_app_write_iz_react(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_recomb_react(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_iz_react_neut(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_recomb_react_neut(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_mom(gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_source_mom(gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_integrated_mom(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_integrated_source_mom(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_field_energy(gkyl_gyrokinetic_multib_app* app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_max_corr_status(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_write_geometry(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_read_geometry(gkyl_gyrokinetic_multib_app *app)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_stat_write(gkyl_gyrokinetic_multib_app* app)
{
  // TO DO
}

struct gkyl_update_status gkyl_gyrokinetic_multib_update(gkyl_gyrokinetic_multib_app* app, double dt)
{
  // TO DO
  return (struct gkyl_update_status) { };
}

struct gkyl_gyrokinetic_stat gkyl_gyrokinetic_multib_app_stat(gkyl_gyrokinetic_multib_app* app)
{
  return app->stat;
}

void gkyl_gyrokinetic_multib_app_species_ktm_rhs(gkyl_gyrokinetic_multib_app* app, int update_vol_term)
{
  // TO DO
}

void gkyl_gyrokinetic_multib_app_release(gkyl_gyrokinetic_multib_app* mbapp)
{
  if (mbapp->singleb_apps) {
    for (int i=0; i<mbapp->num_local_blocks; ++i)
      gkyl_gyrokinetic_app_release(mbapp->singleb_apps[i]);
    gkyl_free(mbapp->singleb_apps);
  }  

  int num_blocks = gkyl_block_geom_num_blocks(mbapp->block_geom);

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(mbapp->decomp[i]);
  gkyl_free(mbapp->decomp);

  for (int i=0; i<num_blocks; ++i)
    gkyl_comm_release(mbapp->block_comms[i]);
  gkyl_free(mbapp->block_comms);
  gkyl_comm_release(mbapp->comm);

  gkyl_rrobin_decomp_release(mbapp->round_robin);
  
  gkyl_block_geom_release(mbapp->block_geom);
  gkyl_block_topo_release(mbapp->block_topo);
  gkyl_free(mbapp->local_blocks);    
  
  gkyl_free(mbapp);
}
