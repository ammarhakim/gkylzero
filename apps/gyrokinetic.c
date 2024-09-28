#include <stdarg.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_basis.h>
#include <gkyl_comm_io.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>
#include <gkyl_null_comm.h>

#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_app_priv.h>

#include <mpack.h>

void
gyrokinetic_cuts_check(struct gkyl_gyrokinetic_app* app, struct gkyl_comm *comm, const int *cuts, FILE *iostream)
{
  // A temporary function that checks the consistency of the communicator and
  // cuts provided for a simulation (e.g at the moment we only decompose along
  // z, and we need to check that cuts meets that requirement).
  int cdim = app->cdim;

  // Create decomposition.
  int cuts_used[cdim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < cdim; d++)
    cuts_used[d] = cuts[d];
#else
  for (int d = 0; d < cdim; d++) cuts_used[d] = 1;
#endif

  int comm_rank, comm_size;
  gkyl_comm_get_rank(comm, &comm_rank);
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < cdim; d++) ncuts *= cuts_used[d];

  if (ncuts != comm_size) {
    if (comm_rank == 0)
      fprintf(iostream, "\n*** Number of ranks, %d, does not match total cuts, %d!\n\n", comm_size, ncuts);
    assert(false);
  }

  for (int d = 0; d < cdim - 1; d++) {
    if (cuts_used[d] > 1) {
      if (comm_rank == 0)
        fprintf(iostream,
          "\n*** Parallelization only allowed in z. Number of ranks, %d, in direction %d cannot be > 1!\n\n", cuts_used[d], d);
      assert(false);
    }
  }
}

// returned gkyl_array_meta must be freed using gyrokinetic_array_meta_release
static struct gkyl_array_meta*
gyrokinetic_array_meta_new(struct gyrokinetic_output_meta meta)
{
  struct gkyl_array_meta *mt = gkyl_malloc(sizeof(*mt));

  mt->meta_sz = 0;
  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mt->meta, &mt->meta_sz);

  // add some data to mpack
  mpack_build_map(&writer);
  
  mpack_write_cstr(&writer, "time");
  mpack_write_double(&writer, meta.stime);

  mpack_write_cstr(&writer, "frame");
  mpack_write_i64(&writer, meta.frame);

  mpack_write_cstr(&writer, "polyOrder");
  mpack_write_i64(&writer, meta.poly_order);

  mpack_write_cstr(&writer, "basisType");
  mpack_write_cstr(&writer, meta.basis_type);

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

static void
gyrokinetic_array_meta_release(struct gkyl_array_meta *mt)
{
  if (!mt) return;
  MPACK_FREE(mt->meta);
  gkyl_free(mt);
}

static struct gyrokinetic_output_meta
gyrokinetic_meta_from_mpack(struct gkyl_array_meta *mt)
{
  struct gyrokinetic_output_meta meta = { .frame = 0, .stime = 0.0 };

  if (mt->meta_sz > 0) {
    mpack_tree_t tree;
    mpack_tree_init_data(&tree, mt->meta, mt->meta_sz);
    mpack_tree_parse(&tree);
    mpack_node_t root = mpack_tree_root(&tree);

    mpack_node_t tm_node = mpack_node_map_cstr(root, "time");
    meta.stime = mpack_node_double(tm_node);

    mpack_node_t fr_node = mpack_node_map_cstr(root, "frame");
    meta.frame = mpack_node_i64(fr_node);

    mpack_node_t po_node = mpack_node_map_cstr(root, "polyOrder");
    meta.poly_order = mpack_node_i64(po_node);

    mpack_node_t bt_node = mpack_node_map_cstr(root, "basisType");
    char *basis_type = mpack_node_cstr_alloc(bt_node, 64);
    strcpy(meta.basis_type_nm, basis_type);
    meta.basis_type = meta.basis_type_nm;
    MPACK_FREE(basis_type);

    mpack_tree_destroy(&tree);
  }
  return meta;
}

gkyl_gyrokinetic_app*
gkyl_gyrokinetic_app_new(struct gkyl_gk *gk)
{
  disable_denorm_float();

  assert(gk->num_species <= GKYL_MAX_SPECIES);

  gkyl_gyrokinetic_app *app = gkyl_malloc(sizeof(gkyl_gyrokinetic_app));

  int cdim = app->cdim = gk->cdim;
  int vdim = app->vdim = gk->vdim;
  int pdim = cdim+vdim;
  int poly_order = app->poly_order = gk->poly_order;
  int ns = app->num_species = gk->num_species;
  int neuts = app->num_neut_species = gk->num_neut_species;

  double cfl_frac = gk->cfl_frac == 0 ? 1.0 : gk->cfl_frac;
  app->cfl = cfl_frac;

#ifdef GKYL_HAVE_CUDA
  app->use_gpu = gk->parallelism.use_gpu;
#else
  app->use_gpu = false; // can't use GPUs if we don't have them!
#endif

  app->enforce_positivity = gk->enforce_positivity;

  app->num_periodic_dir = gk->num_periodic_dir;
  for (int d=0; d<cdim; ++d)
    app->periodic_dirs[d] = gk->periodic_dirs[d];

  strcpy(app->name, gk->name);
  app->tcurr = 0.0; // reset on init

  if (app->use_gpu) {
    // allocate device basis if we are using GPUs
    app->basis_on_dev.basis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    app->basis_on_dev.neut_basis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    app->basis_on_dev.confBasis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  }
  else {
    app->basis_on_dev.basis = &app->basis;
    app->basis_on_dev.neut_basis = &app->neut_basis;
    app->basis_on_dev.confBasis = &app->confBasis;
  }

  // basis functions
  switch (gk->basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);
      if (poly_order > 1) {
        gkyl_cart_modal_serendip(&app->basis, pdim, poly_order);
        gkyl_cart_modal_serendip(&app->neut_basis, pdim+1, poly_order); // neutral species are 3v
      }
      else if (poly_order == 1) {
        gkyl_cart_modal_gkhybrid(&app->basis, cdim, vdim); // p=2 in vparallel
        gkyl_cart_modal_hybrid(&app->neut_basis, cdim, vdim+1); // p=2 in v for neutral species
      }

      if (app->use_gpu) {
        gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.confBasis, cdim, poly_order);
        if (poly_order > 1) {
          gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.basis, pdim, poly_order);
          gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.neut_basis, pdim+1, poly_order); // neutral species are 3v
        }
        else if (poly_order == 1) {
          gkyl_cart_modal_gkhybrid_cu_dev(app->basis_on_dev.basis, cdim, vdim); // p=2 in vparallel
          gkyl_cart_modal_hybrid_cu_dev(app->basis_on_dev.neut_basis, cdim, vdim+1); // p=2 in v for neutral species
        }
      }
      break;
    default:
      assert(false);
      break;
  }

  gkyl_rect_grid_init(&app->grid, cdim, gk->lower, gk->upper, gk->cells);

  int ghost[] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&app->grid, ghost, &app->global_ext, &app->global);

  if (gk->parallelism.comm == 0) {
    int cuts[3] = { 1, 1, 1 };
    app->decomp = gkyl_rect_decomp_new_from_cuts(cdim, cuts, &app->global);
    
    app->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = app->decomp,
        .use_gpu = app->use_gpu
      }
    );
    
    // Global and local ranges are same, and so just copy them.
    memcpy(&app->local, &app->global, sizeof(struct gkyl_range));
    memcpy(&app->local_ext, &app->global_ext, sizeof(struct gkyl_range));
  }
  else {
    gyrokinetic_cuts_check(app, gk->parallelism.comm, gk->parallelism.cuts, stdout);

    // Create decomp.
    app->decomp = gkyl_rect_decomp_new_from_cuts(app->cdim, gk->parallelism.cuts, &app->global);

    // Create a new communicator with the decomposition in it.
    app->comm = gkyl_comm_split_comm(gk->parallelism.comm, 0, app->decomp);

    // Create local and local_ext.
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    gkyl_create_ranges(&app->decomp->ranges[rank], ghost, &app->local_ext, &app->local);
  }

  // local skin and ghost ranges for configuration space fields
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&app->lower_skin[dir], &app->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &app->local_ext, ghost); 
    gkyl_skin_ghost_ranges(&app->upper_skin[dir], &app->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &app->local_ext, ghost);
  }

  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);

  // Configuration space geometry initialization

  // Initialize the input struct from user side input struct
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = gk->geometry.geometry_id,
    .c2p_ctx = gk->geometry.c2p_ctx,
    .mapc2p = gk->geometry.mapc2p,
    .bmag_ctx = gk->geometry.bmag_ctx,
    .bmag_func = gk->geometry.bmag_func,
    .efit_info = gk->geometry.efit_info,
    .tok_grid_info = gk->geometry.tok_grid_info,
    .mirror_grid_info = gk->geometry.mirror_grid_info,
    .grid = app->grid,
    .local = app->local,
    .local_ext = app->local_ext,
    .global = app->global,
    .global_ext = app->global_ext,
    .basis = app->confBasis,
  };
  for(int i = 0; i<3; i++)
    geometry_inp.world[i] = gk->geometry.world[i];

  if(app->cdim < 3){
    geometry_inp.geo_grid = gkyl_gk_geometry_augment_grid(app->grid, geometry_inp);
    switch (gk->basis_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
        gkyl_cart_modal_serendip(&geometry_inp.geo_basis, 3, poly_order);
        break;
      default:
        assert(false);
        break;
    }

    int ghost[] = { 1, 1, 1 };
    gkyl_create_grid_ranges(&geometry_inp.geo_grid, ghost, &geometry_inp.geo_global_ext, &geometry_inp.geo_global);
    if (comm_sz > 1) {
      // create local and local_ext from user-supplied local range
      gkyl_gk_geometry_augment_local(&app->local, ghost, &geometry_inp.geo_local_ext, &geometry_inp.geo_local);
    }
    else {
      // global and local ranges are same, and so just copy
      memcpy(&geometry_inp.geo_local, &geometry_inp.geo_global, sizeof(struct gkyl_range));
      memcpy(&geometry_inp.geo_local_ext, &geometry_inp.geo_global_ext, sizeof(struct gkyl_range));
    }

  }
  else{
    geometry_inp.geo_grid = app->grid;
    geometry_inp.geo_local = app->local;
    geometry_inp.geo_local_ext = app->local_ext;
    geometry_inp.geo_global = app->global;
    geometry_inp.geo_global_ext = app->global_ext;
    geometry_inp.geo_basis = app->confBasis;
  }

  struct gk_geometry* gk_geom_3d;
  switch (geometry_inp.geometry_id) {
    case GKYL_GEOMETRY_FROMFILE:
      gk_geom_3d = gkyl_gk_geometry_new(app->gk_geom, &geometry_inp, false);
      break;
    case GKYL_TOKAMAK:
      gk_geom_3d = gkyl_gk_geometry_tok_new(&geometry_inp);
      break;
    case GKYL_MIRROR:
      gk_geom_3d = gkyl_gk_geometry_mirror_new(&geometry_inp);
      break;
    case GKYL_MAPC2P:
      gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_inp);
      break;
  }

  // deflate geometry if necessary
  if (geometry_inp.geometry_id != GKYL_GEOMETRY_FROMFILE) {
    if(app->cdim < 3)
      app->gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_inp);
    else
      app->gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);
  }
  else {
    app->gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);
    gkyl_gyrokinetic_app_read_geometry(app);
  }

  gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry

  gkyl_gk_geometry_bmag_mid(app->gk_geom); // set bmag mid
  int bcast_rank = comm_sz/2;
  gkyl_comm_array_bcast_host(app->comm, app->gk_geom->bmag_mid, app->gk_geom->bmag_mid, bcast_rank);
  
  // If we are on the gpu, copy from host
  if (app->use_gpu) {
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(app->gk_geom, &geometry_inp, app->use_gpu);
    gkyl_gk_geometry_release(app->gk_geom);
    app->gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  gkyl_gyrokinetic_app_write_geometry(app);

  // allocate space to store species and neutral species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct gk_species[ns])) : 0;
  app->neut_species = neuts>0 ? gkyl_malloc(sizeof(struct gk_neut_species[neuts])) : 0;

  // set info for each species: this needs to be done here as we need
  // to access species name from gk_species_init
  for (int i=0; i<ns; ++i)
    app->species[i].info = gk->species[i];

  // set info for each neutral species: this needs to be done here as we need
  // to access species name from gk_species_init
  for (int i=0; i<neuts; ++i)
    app->neut_species[i].info = gk->neut_species[i];

  app->update_field = !gk->skip_field; // note inversion of truth value (default: update field)
  app->field = gk_field_new(gk, app); // initialize field, even if we are skipping field updates

  // initialize each species
  for (int i=0; i<ns; ++i) 
    gk_species_init(gk, app, &app->species[i]);

  // initialize each neutral species
  for (int i=0; i<neuts; ++i) 
    gk_neut_species_init(gk, app, &app->neut_species[i]);

  // initialize each species cross-collisions terms: this has to be done here
  // as need pointers to colliding species' collision objects
  // allocated in gk_species_init and gk_neut_species_init
  for (int i=0; i<ns; ++i) {
    struct gk_species *gk_s = &app->species[i];

    // initialize cross-species collisions (e.g, LBO or BGK)
    if (gk_s->lbo.collision_id == GKYL_LBO_COLLISIONS) {
      if (gk_s->lbo.num_cross_collisions) {
        gk_species_lbo_cross_init(app, &app->species[i], &gk_s->lbo);
      }
    }
    if (gk_s->bgk.collision_id == GKYL_BGK_COLLISIONS) {
      if (gk_s->bgk.num_cross_collisions) {
        gk_species_bgk_cross_init(app, &app->species[i], &gk_s->bgk);
      }
    }
    // initialize cross-species reactions with plasma species (e.g., ionization, recombination, or charge exchange)
    if (gk_s->react.num_react) {
      gk_species_react_cross_init(app, &app->species[i], &gk_s->react);
    }
    // initialize cross-species reactions with neutral species (e.g., ionization, recombination, or charge exchange)
    if (gk_s->react_neut.num_react) {
      gk_species_react_cross_init(app, &app->species[i], &gk_s->react_neut);
    }
    // initial radiation (e.g., line radiation from cross-collisions of electrons with ions)
    if (gk_s->info.radiation.radiation_id == GKYL_GK_RADIATION) {
      gk_species_radiation_init(app, &app->species[i], &gk_s->rad);
    }
  }

  // initialize neutral species cross-species reactions with plasma species
  for (int i=0; i<neuts; ++i) {
    if (app->neut_species[i].react_neut.num_react) {
      gk_neut_species_react_cross_init(app, &app->neut_species[i], &app->neut_species[i].react_neut);
    }
  }

  // initialize each plasma species and neutral species source terms
  // This has to be done here as sources may initialize a boundary 
  // flux updater for their source species
  for (int i=0; i<ns; ++i) {
    gk_species_source_init(app, &app->species[i], &app->species[i].src);
  }
  for (int i=0; i<neuts; ++i) {
    gk_neut_species_source_init(app, &app->neut_species[i], &app->neut_species[i].src);
  }

  if (app->enforce_positivity) {
    // Number of density of the positivity shift added over all the ions.
    app->ps_delta_m0_ions = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  }

  // initialize stat object
  app->stat = (struct gkyl_gyrokinetic_stat) {
    .use_gpu = app->use_gpu,
    .stage_2_dt_diff = { DBL_MAX, 0.0 },
    .stage_3_dt_diff = { DBL_MAX, 0.0 },
  };

  return app;
}

// Compute fields.
static void
calc_field(gkyl_gyrokinetic_app* app, double tcurr, const struct gkyl_array *fin[])
{
  if (app->update_field) {
    // Compute electrostatic potential from gyrokinetic Poisson's equation.
    gk_field_accumulate_rho_c(app, app->field, fin);

    // Compute ambipolar potential sheath values if using adiabatic electrons
    // done here as the RHS update for all species should be complete before
    // boundary fluxes are computed (ion fluxes needed for sheath values) 
    // and these boundary fluxes are stored temporarily in ghost cells of RHS
    if (app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN)
      gk_field_calc_ambi_pot_sheath_vals(app, app->field);

    // Compute biased wall potential if present and time-dependent.
    // Note: biased wall potential use eval_on_nodes. 
    // so does copy to GPU every call if app->use_gpu = true.
    if (app->field->phi_wall_lo_evolve || app->field->phi_wall_up_evolve)
      gk_field_calc_phi_wall(app, app->field, tcurr);

    // Solve the field equation.
    gk_field_rhs(app, app->field);
  }
}

// Compute fields and apply BCs.
static void
calc_field_and_apply_bc(gkyl_gyrokinetic_app* app, double tcurr, struct gkyl_array *distf[], struct gkyl_array *distf_neut[])
{

  // Compute the field.
  calc_field(app, tcurr, (const struct gkyl_array **) distf);

  // Apply boundary conditions.
  for (int i=0; i<app->num_species; ++i) {
    gk_species_apply_bc(app, &app->species[i], distf[i]);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    if (!app->neut_species[i].info.is_static) {
      gk_neut_species_apply_bc(app, &app->neut_species[i], distf_neut[i]);
    }
  }

}

struct gk_species *
gk_find_species(const gkyl_gyrokinetic_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return &app->species[i];
  return 0;
}

int
gk_find_species_idx(const gkyl_gyrokinetic_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return i;
  return -1;
}

struct gk_neut_species *
gk_find_neut_species(const gkyl_gyrokinetic_app *app, const char *nm)
{
  for (int i=0; i<app->num_neut_species; ++i)
    if (strcmp(nm, app->neut_species[i].info.name) == 0)
      return &app->neut_species[i];
  return 0;
}

int
gk_find_neut_species_idx(const gkyl_gyrokinetic_app *app, const char *nm)
{
  for (int i=0; i<app->num_neut_species; ++i)
    if (strcmp(nm, app->neut_species[i].info.name) == 0)
      return i;
  return -1;
}

void
gkyl_gyrokinetic_app_apply_ic(gkyl_gyrokinetic_app* app, double t0)
{
  app->tcurr = t0;
  for (int i=0; i<app->num_species; ++i)
    gkyl_gyrokinetic_app_apply_ic_species(app, i, t0);
  for (int i=0; i<app->num_neut_species; ++i)
    gkyl_gyrokinetic_app_apply_ic_neut_species(app, i, t0);

  for (int i=0; i<app->num_species; ++i)
    gkyl_gyrokinetic_app_apply_ic_cross_species(app, i, t0);

  // Compute the fields and apply BCs.
  struct gkyl_array *distf[app->num_species];
  struct gkyl_array *distf_neut[app->num_neut_species];
  for (int i=0; i<app->num_species; ++i) {
    distf[i] = app->species[i].f;
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    distf_neut[i] = app->neut_species[i].f;
  }
  if (app->update_field && app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
    for (int i=0; i<app->num_species; ++i) {
      struct gk_species *s = &app->species[i];

      // Compute advection speeds so we can compute the initial boundary flux.
      gkyl_dg_calc_gyrokinetic_vars_alpha_surf(s->calc_gk_vars, 
        &app->local, &s->local, &s->local_ext, app->field->phi_smooth,
        s->alpha_surf, s->sgn_alpha_surf, s->const_sgn_alpha);

      // Compute and store (in the ghost cell of of out) the boundary fluxes.
      // NOTE: this overwrites ghost cells that may be used for sourcing.
      gk_species_bflux_rhs(app, s, &s->bflux, distf[i], distf[i]);
    }
  }
  calc_field_and_apply_bc(app, 0., distf, distf_neut);
}

void
gkyl_gyrokinetic_app_apply_ic_species(gkyl_gyrokinetic_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  struct gk_species *gk_s = &app->species[sidx];

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  gk_species_apply_ic(app, gk_s, t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);
}

void
gkyl_gyrokinetic_app_apply_ic_neut_species(gkyl_gyrokinetic_app* app, int sidx, double t0)
{
  assert(sidx < app->num_neut_species);

  struct gk_neut_species *gk_ns = &app->neut_species[sidx];

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  gk_neut_species_apply_ic(app, gk_ns, t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);
}

void
gkyl_gyrokinetic_app_apply_ic_cross_species(gkyl_gyrokinetic_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  struct gk_species *gk_s = &app->species[sidx];

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  gk_species_apply_ic_cross(app, gk_s, t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);
}

//
// ............. Geometry outputs ............... //
// 
static void
gyrokinetic_app_geometry_copy_and_write(gkyl_gyrokinetic_app* app, struct gkyl_array *arr,
  struct gkyl_array *arr_host, char *varNm, struct gkyl_array_meta *mt)
{
  gkyl_array_copy(arr_host, arr);

  const char *fmt = "%s-%s.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, varNm);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, varNm);

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, arr_host, fileNm);
}

void
gkyl_gyrokinetic_app_write_geometry(gkyl_gyrokinetic_app* app)
{
  struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = 0,
      .stime = 0,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  // Gather geo into a global array
  struct gkyl_array* arr_ho1 = mkarr(false,   app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array* arr_ho3 = mkarr(false, 3*app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array* arr_ho6 = mkarr(false, 6*app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array* arr_ho9 = mkarr(false, 9*app->confBasis.num_basis, app->local_ext.volume);

  struct timespec wtm = gkyl_wall_clock();
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->mc2p        , arr_ho3, "mapc2p", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->bmag        , arr_ho1, "bmag", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->g_ij        , arr_ho6, "g_ij", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->dxdz        , arr_ho9, "dxdz", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->dzdx        , arr_ho9, "dzdx", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->normals     , arr_ho9, "normals", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->jacobgeo    , arr_ho1, "jacobgeo", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->jacobgeo_inv, arr_ho1, "jacobgeo_inv", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->gij         , arr_ho6, "gij", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->b_i         , arr_ho3, "b_i", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->bcart       , arr_ho3, "bcart", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->cmag        , arr_ho1, "cmag", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->jacobtot    , arr_ho1, "jacobtot", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->jacobtot_inv, arr_ho1, "jacobtot_inv", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->bmag_inv    , arr_ho1, "bmag_inv", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->bmag_inv_sq , arr_ho1, "bmag_inv_sq", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->gxxj        , arr_ho1, "gxxj", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->gxyj        , arr_ho1, "gxyj", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->gyyj        , arr_ho1, "gyyj", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->gxzj        , arr_ho1, "gxzj", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->eps2        , arr_ho1, "eps2", mt);
  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.nio += 21;

  // Write out global geometry on rank 0
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Create Nodal Range and Grid and Write Nodal Coordinates
    struct gkyl_range nrange;
    gkyl_gk_geometry_init_nodal_range(&nrange, &app->global, app->poly_order);
    struct gkyl_array* mc2p_nodal = mkarr(false, 3, nrange.volume);
    struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&app->confBasis, &app->grid, false);
    gkyl_array_copy(arr_ho3, app->gk_geom->mc2p);
    gkyl_nodal_ops_m2n(n2m, &app->confBasis, &app->grid, &nrange, &app->global, 3, mc2p_nodal, arr_ho3);
    struct gkyl_rect_grid ngrid;
    gkyl_gk_geometry_init_nodal_grid(&ngrid, &app->grid, &nrange);

    const char *fmt = "%s-%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, "nodes");
    char fileNm[sz+1]; // ensures no buffer overflow
    sprintf(fileNm, fmt, app->name, "nodes");

    struct timespec wtn = gkyl_wall_clock();
    gkyl_grid_sub_array_write(&ngrid, &nrange, 0,  mc2p_nodal, fileNm);
    app->stat.io_tm += gkyl_time_diff_now_sec(wtn);
    app->stat.nio += 1;

    gkyl_nodal_ops_release(n2m);
    gkyl_array_release(mc2p_nodal);
  }

  gkyl_array_release(arr_ho1);
  gkyl_array_release(arr_ho3);
  gkyl_array_release(arr_ho6);
  gkyl_array_release(arr_ho9);

  gyrokinetic_array_meta_release(mt);
}

//
// ............. Field outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_field(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  if (app->update_field) {
    struct timespec wst = gkyl_wall_clock();
    // Copy data from device to host before writing it out.
    if (app->use_gpu) {
      gkyl_array_copy(app->field->phi_host, app->field->phi_smooth);
    }

    struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->confBasis.id
      }
    );

    const char *fmt = "%s-field_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);

    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, app->field->phi_host, fileNm);
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;

    gyrokinetic_array_meta_release(mt);
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_calc_field_energy(gkyl_gyrokinetic_app* app, double tm)
{
  if (app->update_field) {
    struct timespec wst = gkyl_wall_clock();
    gk_field_calc_energy(app, tm, app->field);
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_write_field_energy(gkyl_gyrokinetic_app* app)
{
  if (app->update_field) {
    struct timespec wst = gkyl_wall_clock();
    // write out diagnostic moments
    const char *fmt = "%s-field_energy.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name);

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);

    if (rank == 0) {
      struct timespec wtm = gkyl_wall_clock();
      if (app->field->is_first_energy_write_call) {
        // write to a new file (this ensure previous output is removed)
        gkyl_dynvec_write(app->field->integ_energy, fileNm);
        app->field->is_first_energy_write_call = false;
      }
      else {
        // append to existing file
        gkyl_dynvec_awrite(app->field->integ_energy, fileNm);
      }
      app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.nio += 1;
    }
    gkyl_dynvec_clear(app->field->integ_energy);
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

//
// ............. Species outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gk_species *gk_s = &app->species[sidx];

  struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }
  );

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gk_s->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_s->info.name, frame);

  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(gk_s->f_host, gk_s->f);
  }

  struct timespec wtm = gkyl_wall_clock();
  gkyl_comm_array_write(gk_s->comm, &gk_s->grid, &gk_s->local, mt, gk_s->f_host, fileNm);
  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.nio += 1;

  gyrokinetic_array_meta_release(mt);  
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_gyrokinetic_app_write_neut_species(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_neut_species *gk_ns = &app->neut_species[sidx];

  if (!gk_ns->info.is_static || frame == 0) {
    struct timespec wst = gkyl_wall_clock();
    struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->basis.id
      }
    );

    const char *fmt = "%s-%s_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gk_ns->info.name, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_ns->info.name, frame);

    // copy data from device to host before writing it out
    if (app->use_gpu) {
      gkyl_array_copy(gk_ns->f_host, gk_ns->f);
    }

    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(gk_ns->comm, &gk_ns->grid, &gk_ns->local, mt, gk_ns->f_host, fileNm);
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;

    gyrokinetic_array_meta_release(mt);  
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_write_species_mom(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gk_species *gks = &app->species[sidx];

  struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  for (int m=0; m<gks->info.num_diag_moments; ++m) {
    gk_species_moment_calc(&gks->moms[m], gks->local, app->local, gks->f);
    app->stat.nmom += 1;

    const char *fmt = "%s-%s_%s_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name,
      gks->info.diag_moments[m], frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name,
      gks->info.diag_moments[m], frame);

    // Rescale moment by inverse of Jacobian if not already re-scaled 
    if (!gks->moms[m].is_bimaxwellian_moms && !gks->moms[m].is_maxwellian_moms) {
      gkyl_dg_div_op_range(gks->moms[m].mem_geo, app->confBasis, 
        0, gks->moms[m].marr, 0, gks->moms[m].marr, 0, 
        app->gk_geom->jacobgeo, &app->local);  
    }    

    if (app->use_gpu) {
      gkyl_array_copy(gks->moms[m].marr_host, gks->moms[m].marr);
    }

    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
      gks->moms[m].marr_host, fileNm);
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;
  }

  if (app->enforce_positivity) {
    // We placed the change in f from the positivity shift in fnew.
    gk_species_moment_calc(&gks->ps_moms, gks->local, app->local, gks->fnew);
    app->stat.nmom += 1;

    const char *fmt = "%s-%s_positivity_shift_FourMoments_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

    // Rescale moment by inverse of Jacobian.
    gkyl_dg_div_op_range(gks->ps_moms.mem_geo, app->confBasis, 
      0, gks->ps_moms.marr, 0, gks->ps_moms.marr, 0, 
      app->gk_geom->jacobgeo, &app->local);  

    if (app->use_gpu) {
      gkyl_array_copy(gks->ps_moms.marr_host, gks->ps_moms.marr);
    }

    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
      gks->ps_moms.marr_host, fileNm);
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;
  }
  gyrokinetic_array_meta_release(mt);   
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_gyrokinetic_app_write_neut_species_mom(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_neut_species *gk_ns = &app->neut_species[sidx];

  if (!gk_ns->info.is_static || frame == 0) {
    struct timespec wst = gkyl_wall_clock();

    struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->confBasis.id
      }
    );

    for (int m=0; m<gk_ns->info.num_diag_moments; ++m) {
      gk_neut_species_moment_calc(&gk_ns->moms[m], gk_ns->local, app->local, gk_ns->f);
      app->stat.nmom += 1;

      const char *fmt = "%s-%s_%s_%d.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gk_ns->info.name,
        gk_ns->info.diag_moments[m], frame);
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_ns->info.name,
        gk_ns->info.diag_moments[m], frame);

      // Rescale moment by inverse of Jacobian
      gkyl_dg_div_op_range(gk_ns->moms[m].mem_geo, app->confBasis, 
        0, gk_ns->moms[m].marr, 0, gk_ns->moms[m].marr, 0, 
        app->gk_geom->jacobgeo, &app->local);      

      if (app->use_gpu) {
        gkyl_array_copy(gk_ns->moms[m].marr_host, gk_ns->moms[m].marr);
      }

      struct timespec wtm = gkyl_wall_clock();
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
        gk_ns->moms[m].marr_host, fileNm);
      app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.nio += 1;
    }
    gyrokinetic_array_meta_release(mt);   
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_calc_species_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  struct gk_species *gk_s = &app->species[sidx];

  int vdim = app->vdim;
  double avals_global[2+vdim];

  gk_species_moment_calc(&gk_s->integ_moms, gk_s->local, app->local, gk_s->f); 
  // Reduce (sum) over whole domain, append to diagnostics.
  gkyl_array_reduce_range(gk_s->red_integ_diag, gk_s->integ_moms.marr, GKYL_SUM, &app->local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, 
    gk_s->red_integ_diag, gk_s->red_integ_diag_global);
  if (app->use_gpu) {
    gkyl_cu_memcpy(avals_global, gk_s->red_integ_diag_global, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(avals_global, gk_s->red_integ_diag_global, sizeof(double[2+vdim]));
  }
  gkyl_dynvec_append(gk_s->integ_diag, tm, avals_global);

  if (app->enforce_positivity) {
    // The change in f from the positivity shift is in fnew.
    gk_species_moment_calc(&gk_s->integ_moms, gk_s->local, app->local, gk_s->fnew); 
    // Reduce (sum) over whole domain, append to diagnostics.
    gkyl_array_reduce_range(gk_s->red_integ_diag, gk_s->integ_moms.marr, GKYL_SUM, &app->local);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, 
      gk_s->red_integ_diag, gk_s->red_integ_diag_global);
    if (app->use_gpu) {
      gkyl_cu_memcpy(avals_global, gk_s->red_integ_diag_global, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      memcpy(avals_global, gk_s->red_integ_diag_global, sizeof(double[2+vdim]));
    }
    gkyl_dynvec_append(gk_s->ps_integ_diag, tm, avals_global);
  }

  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_gyrokinetic_app_calc_neut_species_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm)
{
  struct gk_neut_species *gk_ns = &app->neut_species[sidx];

  if (!gk_ns->info.is_static) {
    struct timespec wst = gkyl_wall_clock();

    int vdim = app->vdim+1; // Neutrals are always 3V
    double avals_global[2+vdim];

    gk_neut_species_moment_calc(&gk_ns->integ_moms, gk_ns->local, app->local, gk_ns->f); 
    // reduce to compute sum over whole domain, append to diagnostics
    gkyl_array_reduce_range(gk_ns->red_integ_diag, gk_ns->integ_moms.marr, GKYL_SUM, &app->local);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, 
      gk_ns->red_integ_diag, gk_ns->red_integ_diag_global);
    if (app->use_gpu) {
      gkyl_cu_memcpy(avals_global, gk_ns->red_integ_diag_global, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      memcpy(avals_global, gk_ns->red_integ_diag_global, sizeof(double[2+vdim]));
    }
    gkyl_dynvec_append(gk_ns->integ_diag, tm, avals_global);
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_write_species_integrated_mom(gkyl_gyrokinetic_app *app, int sidx)
{
  struct timespec wst = gkyl_wall_clock();
  struct gk_species *gks = &app->species[sidx];

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s_%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "integrated_moms");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "integrated_moms");

    struct timespec wtm = gkyl_wall_clock();
    if (gks->is_first_integ_write_call) {
      gkyl_dynvec_write(gks->integ_diag, fileNm);
      gks->is_first_integ_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(gks->integ_diag, fileNm);
    }
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;
  }
  gkyl_dynvec_clear(gks->integ_diag);

  if (app->enforce_positivity) {
    if (rank == 0) {
      // Write integrated diagnostic moments.
      const char *fmt = "%s-%s_positivity_shift_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "integrated_moms");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "integrated_moms");

      struct timespec wtm = gkyl_wall_clock();
      if (gks->is_first_ps_integ_write_call) {
        gkyl_dynvec_write(gks->ps_integ_diag, fileNm);
        gks->is_first_ps_integ_write_call = false;
      }
      else {
        gkyl_dynvec_awrite(gks->ps_integ_diag, fileNm);
      }
      app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.nio += 1;
    }
    gkyl_dynvec_clear(gks->ps_integ_diag);
  }

  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_gyrokinetic_app_write_neut_species_integrated_mom(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_neut_species *gk_ns = &app->neut_species[sidx];

  if (!gk_ns->info.is_static) {
    struct timespec wst = gkyl_wall_clock();

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);

    if (rank == 0) {
      // Write integrated diagnostic moments.
      const char *fmt = "%s-%s_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gk_ns->info.name, "integrated_moms");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_ns->info.name, "integrated_moms");

      struct timespec wtm = gkyl_wall_clock();
      if (gk_ns->is_first_integ_write_call) {
        gkyl_dynvec_write(gk_ns->integ_diag, fileNm);
        gk_ns->is_first_integ_write_call = false;
      }
      else {
        gkyl_dynvec_awrite(gk_ns->integ_diag, fileNm);
      }
      app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.nio += 1;
    }
    gkyl_dynvec_clear(gk_ns->integ_diag);
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

//
// ............. Source outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_source(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gk_s = &app->species[sidx];

  if (gk_s->src.source_id && (gk_s->src.evolve || frame == 0)) {
    struct timespec wst = gkyl_wall_clock();
    struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->basis.id
      }
    );

    // Write out the source distribution function
    const char *fmt = "%s-%s_source_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gk_s->info.name, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_s->info.name, frame);

    // copy data from device to host before writing it out
    if (app->use_gpu) {
      gkyl_array_copy(gk_s->src.source_host, gk_s->src.source);
    }

    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(gk_s->comm, &gk_s->grid, &gk_s->local, mt, gk_s->src.source_host, fileNm);
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;

    gyrokinetic_array_meta_release(mt);   
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}


void
gkyl_gyrokinetic_app_write_neut_species_source(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_neut_species *gk_ns = &app->neut_species[sidx];

  if (gk_ns->src.source_id && (gk_ns->src.evolve || frame == 0)) {
    struct timespec wst = gkyl_wall_clock();
    struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->basis.id
      }
    );

    // Write out the source distribution function
    const char *fmt = "%s-%s_source_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gk_ns->info.name, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_ns->info.name, frame);

    // copy data from device to host before writing it out
    if (app->use_gpu) {
      gkyl_array_copy(gk_ns->src.source_host, gk_ns->src.source);
    }

    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(gk_ns->comm, &gk_ns->grid, &gk_ns->local, mt, gk_ns->src.source_host, fileNm);
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;

    gyrokinetic_array_meta_release(mt);   
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_write_species_source_mom(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gk_s = &app->species[sidx];

  if (gk_s->src.source_id && (gk_s->src.evolve || frame == 0)) {
    struct timespec wst = gkyl_wall_clock();
    struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->confBasis.id
      }
    );

    for (int m=0; m<gk_s->info.num_diag_moments; ++m) {
      gk_species_moment_calc(&gk_s->src.moms[m], gk_s->local, app->local, gk_s->src.source);
      app->stat.nmom += 1;

      const char *fmt = "%s-%s_source_%s_%d.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gk_s->info.name,
        gk_s->info.diag_moments[m], frame);
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_s->info.name,
        gk_s->info.diag_moments[m], frame);

      // Rescale moment by inverse of Jacobian
      gkyl_dg_div_op_range(gk_s->moms[m].mem_geo, app->confBasis, 
        0, gk_s->src.moms[m].marr, 0, gk_s->src.moms[m].marr, 0, 
        app->gk_geom->jacobgeo, &app->local);      

      if (app->use_gpu) {
        gkyl_array_copy(gk_s->src.moms[m].marr_host, gk_s->src.moms[m].marr);
      }

      struct timespec wtm = gkyl_wall_clock();
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
        gk_s->src.moms[m].marr_host, fileNm);
      app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.nio += 1;
    }

    gyrokinetic_array_meta_release(mt);   
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_write_neut_species_source_mom(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_neut_species *gk_ns = &app->neut_species[sidx];

  if (gk_ns->src.source_id && (gk_ns->src.evolve || frame == 0)) {
    struct timespec wst = gkyl_wall_clock();

    struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->confBasis.id
      }
    );

    for (int m=0; m<gk_ns->info.num_diag_moments; ++m) {
      gk_neut_species_moment_calc(&gk_ns->src.moms[m], gk_ns->local, app->local, gk_ns->src.source);
      app->stat.nmom += 1;

      const char *fmt = "%s-%s_source_%s_%d.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gk_ns->info.name,
        gk_ns->info.diag_moments[m], frame);
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_ns->info.name,
        gk_ns->info.diag_moments[m], frame);

      // Rescale moment by inverse of Jacobian
      gkyl_dg_div_op_range(gk_ns->moms[m].mem_geo, app->confBasis, 
        0, gk_ns->src.moms[m].marr, 0, gk_ns->src.moms[m].marr, 0, 
        app->gk_geom->jacobgeo, &app->local);      

      if (app->use_gpu) {
        gkyl_array_copy(gk_ns->src.moms[m].marr_host, gk_ns->src.moms[m].marr);
      }

      struct timespec wtm = gkyl_wall_clock();
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
        gk_ns->src.moms[m].marr_host, fileNm);
      app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.nio += 1;
    }

    gyrokinetic_array_meta_release(mt);   
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_calc_species_source_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm)
{
  struct gk_species *gk_s = &app->species[sidx];

  if (gk_s->src.source_id && gk_s->src.evolve) {
    struct timespec wst = gkyl_wall_clock();

    int vdim = app->vdim;
    double avals_global[2+vdim];

    gk_species_moment_calc(&gk_s->src.integ_moms, gk_s->local, app->local, gk_s->src.source); 
    // reduce to compute sum over whole domain, append to diagnostics
    gkyl_array_reduce_range(gk_s->src.red_integ_diag, gk_s->src.integ_moms.marr, GKYL_SUM, &app->local);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, 
      gk_s->src.red_integ_diag, gk_s->src.red_integ_diag_global);
    if (app->use_gpu) {
      gkyl_cu_memcpy(avals_global, gk_s->src.red_integ_diag_global, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      memcpy(avals_global, gk_s->src.red_integ_diag_global, sizeof(double[2+vdim]));
    }
    gkyl_dynvec_append(gk_s->src.integ_diag, tm, avals_global);
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_calc_neut_species_source_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm)
{
  struct gk_neut_species *gk_ns = &app->neut_species[sidx];

  if (gk_ns->src.source_id && gk_ns->src.evolve) {
    struct timespec wst = gkyl_wall_clock();

    int vdim = app->vdim+1; // Neutrals are always 3V
    double avals_global[2+vdim];

    gk_neut_species_moment_calc(&gk_ns->src.integ_moms, gk_ns->local, app->local, gk_ns->src.source); 
    // reduce to compute sum over whole domain, append to diagnostics
    gkyl_array_reduce_range(gk_ns->src.red_integ_diag, gk_ns->src.integ_moms.marr, GKYL_SUM, &app->local);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, 
      gk_ns->src.red_integ_diag, gk_ns->src.red_integ_diag_global);
    if (app->use_gpu) {
      gkyl_cu_memcpy(avals_global, gk_ns->src.red_integ_diag_global, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      memcpy(avals_global, gk_ns->src.red_integ_diag_global, sizeof(double[2+vdim]));
    }
    gkyl_dynvec_append(gk_ns->src.integ_diag, tm, avals_global);
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_write_species_source_integrated_mom(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_species *gks = &app->species[sidx];

  if (gks->src.source_id && gks->src.evolve) {
    struct timespec wst = gkyl_wall_clock();

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out integrated diagnostic moments
      const char *fmt = "%s-%s_source_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "integrated_moms");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "integrated_moms");

      struct timespec wtm = gkyl_wall_clock();
      if (gks->src.is_first_integ_write_call) {
        gkyl_dynvec_write(gks->src.integ_diag, fileNm);
        gks->src.is_first_integ_write_call = false;
      }
      else {
        gkyl_dynvec_awrite(gks->src.integ_diag, fileNm);
      }
      app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.nio += 1;
    }
    gkyl_dynvec_clear(gks->src.integ_diag);
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_write_neut_species_source_integrated_mom(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_neut_species *gk_ns = &app->neut_species[sidx];

  if (gk_ns->src.source_id && gk_ns->src.evolve) {
    struct timespec wst = gkyl_wall_clock();

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out integrated diagnostic moments
      const char *fmt = "%s-%s_source_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gk_ns->info.name, "integrated_moms");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_ns->info.name, "integrated_moms");

      struct timespec wtm = gkyl_wall_clock();
      if (gk_ns->src.is_first_integ_write_call) {
        gkyl_dynvec_write(gk_ns->src.integ_diag, fileNm);
        gk_ns->src.is_first_integ_write_call = false;
      }
      else {
        gkyl_dynvec_awrite(gk_ns->src.integ_diag, fileNm);
      }
      app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.nio += 1;
    }
    gkyl_dynvec_clear(gk_ns->src.integ_diag);
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

//
// ............. Collision outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_lbo_mom(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gk_s = &app->species[sidx];

  if (gk_s->lbo.collision_id == GKYL_LBO_COLLISIONS && gk_s->lbo.write_diagnostics) {
    struct timespec wst = gkyl_wall_clock();

    struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->confBasis.id
      }
    );

    // Construct the file handles for collision frequency and primitive moments
    const char *fmt_prim = "%s-%s_prim_moms_%d.gkyl";
    int sz_prim = gkyl_calc_strlen(fmt_prim, app->name, gk_s->info.name, frame);
    char fileNm_prim[sz_prim+1]; // ensures no buffer overflow
    snprintf(fileNm_prim, sizeof fileNm_prim, fmt_prim, app->name, gk_s->info.name, frame);

    // Compute primitive moments
    const struct gkyl_array *fin[app->num_species];
    gk_species_lbo_moms(app, gk_s, &gk_s->lbo, gk_s->f);

    // copy data from device to host before writing it out
    if (app->use_gpu) {  
      gkyl_array_copy(gk_s->lbo.prim_moms_host, gk_s->lbo.prim_moms);
    }

    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gk_s->lbo.prim_moms_host, fileNm_prim);
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 1;

    // Uncomment the following to write out nu_sum and nu_prim_moms
    /*const char *fmt = "%s-%s_nu_sum_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gk_s->info.name, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_s->info.name, frame);
    
    const char *fmt_nu_prim = "%s-%s_nu_prim_moms_%d.gkyl";
    int sz_nu_prim = gkyl_calc_strlen(fmt_nu_prim, app->name, gk_s->info.name, frame);
    char fileNm_nu_prim[sz_nu_prim+1]; // ensures no buffer overflow
    snprintf(fileNm_nu_prim, sizeof fileNm_nu_prim, fmt_nu_prim, app->name, gk_s->info.name, frame);
    
    if (gk_s->lbo.num_cross_collisions)
      gk_species_lbo_cross_moms(app, gk_s, &gk_s->lbo, gk_s->f);
    
    // copy data from device to host before writing it out
    if (app->use_gpu) {
      gkyl_array_copy(gk_s->lbo.nu_sum_host, gk_s->lbo.nu_sum);
      gkyl_array_copy(gk_s->lbo.nu_prim_moms_host, gk_s->lbo.nu_prim_moms);
    }
    
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gk_s->lbo.nu_sum_host, fileNm);
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gk_s->lbo.nu_prim_moms_host, fileNm_nu_prim);*/

    gyrokinetic_array_meta_release(mt); 
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_write_species_bgk_max_corr_status(gkyl_gyrokinetic_app* app, int sidx)
{
  struct gk_species *gks = &app->species[sidx];

  if (gks->bgk.collision_id == GKYL_BGK_COLLISIONS && gks->bgk.write_diagnostics) {
    struct timespec wst = gkyl_wall_clock();

    // write out diagnostic moments
    const char *fmt = "%s-%s-%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "corr-max-stat");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "corr-max-stat");

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);

    if (rank == 0) {
      struct timespec wtm = gkyl_wall_clock();
      if (gks->bgk.is_first_corr_status_write_call) {
        // write to a new file (this ensure previous output is removed)
        gkyl_dynvec_write(gks->bgk.corr_stat, fileNm);
        gks->bgk.is_first_corr_status_write_call = false;
      }
      else {
        // append to existing file
        gkyl_dynvec_awrite(gks->bgk.corr_stat, fileNm);
      }
      app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.nio += 1;
    }
    gkyl_dynvec_clear(gks->bgk.corr_stat);
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

//
// ............. Radiation outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_rad_drag(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gk_s = &app->species[sidx];

  if (gk_s->rad.radiation_id == GKYL_GK_RADIATION) {
    struct timespec wst = gkyl_wall_clock();
    struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->basis.id      
      }
    );

    // Construct the file handles for vparallel and mu drag
    const char *fmt_nvnu_surf = "%s-%s_radiation_nvnu_surf_%d.gkyl";
    int sz_nvnu_surf = gkyl_calc_strlen(fmt_nvnu_surf, app->name, gk_s->info.name, frame);
    char fileNm_nvnu_surf[sz_nvnu_surf+1]; // ensures no buffer overflow
    snprintf(fileNm_nvnu_surf, sizeof fileNm_nvnu_surf, fmt_nvnu_surf, app->name, gk_s->info.name, frame);

    const char *fmt_nvnu = "%s-%s_radiation_nvnu_%d.gkyl";
    int sz_nvnu = gkyl_calc_strlen(fmt_nvnu, app->name, gk_s->info.name, frame);
    char fileNm_nvnu[sz_nvnu+1]; // ensures no buffer overflow
    snprintf(fileNm_nvnu, sizeof fileNm_nvnu, fmt_nvnu, app->name, gk_s->info.name, frame);

    const char *fmt_nvsqnu_surf = "%s-%s_radiation_nvsqnu_surf_%d.gkyl";
    int sz_nvsqnu_surf = gkyl_calc_strlen(fmt_nvsqnu_surf, app->name, gk_s->info.name, frame);
    char fileNm_nvsqnu_surf[sz_nvsqnu_surf+1]; // ensures no buffer overflow
    snprintf(fileNm_nvsqnu_surf, sizeof fileNm_nvsqnu_surf, fmt_nvsqnu_surf, app->name, gk_s->info.name, frame);

    const char *fmt_nvsqnu = "%s-%s_radiation_nvsqnu_%d.gkyl";
    int sz_nvsqnu = gkyl_calc_strlen(fmt_nvsqnu, app->name, gk_s->info.name, frame);
    char fileNm_nvsqnu[sz_nvsqnu+1]; // ensures no buffer overflow
    snprintf(fileNm_nvsqnu, sizeof fileNm_nvsqnu, fmt_nvsqnu, app->name, gk_s->info.name, frame);

    // Compute radiation drag coefficients
    const struct gkyl_array *fin_neut[app->num_neut_species];
    const struct gkyl_array *fin[app->num_species];
    for (int i=0; i<app->num_species; ++i) 
      fin[i] = app->species[i].f;
    for (int i=0; i<app->num_neut_species; ++i)
      fin_neut[i] = app->neut_species[i].f;

    gk_species_radiation_moms(app, gk_s, &gk_s->rad, fin, fin_neut);

    // copy data from device to host before writing it out
    if (app->use_gpu) {
      gkyl_array_copy(gk_s->rad.nvnu_surf_host, gk_s->rad.nvnu_surf);
      gkyl_array_copy(gk_s->rad.nvnu_host, gk_s->rad.nvnu);
      gkyl_array_copy(gk_s->rad.nvsqnu_surf_host, gk_s->rad.nvsqnu_surf);
      gkyl_array_copy(gk_s->rad.nvsqnu_host, gk_s->rad.nvsqnu);
    }

    struct timespec wtm = gkyl_wall_clock();
    gkyl_comm_array_write(gk_s->comm, &gk_s->grid, &gk_s->local, mt, gk_s->rad.nvnu_surf_host, fileNm_nvnu_surf);
    gkyl_comm_array_write(gk_s->comm, &gk_s->grid, &gk_s->local, mt, gk_s->rad.nvnu_host, fileNm_nvnu);
    gkyl_comm_array_write(gk_s->comm, &gk_s->grid, &gk_s->local, mt, gk_s->rad.nvsqnu_surf_host, fileNm_nvsqnu_surf);
    gkyl_comm_array_write(gk_s->comm, &gk_s->grid, &gk_s->local, mt, gk_s->rad.nvsqnu_host, fileNm_nvsqnu);
    app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.nio += 4;

    gyrokinetic_array_meta_release(mt);   
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_write_species_rad_emissivity(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gk_s = &app->species[sidx];

  if (gk_s->rad.radiation_id == GKYL_GK_RADIATION) {
    struct timespec wst = gkyl_wall_clock();
    struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->confBasis.id
      }
    );
    
    const struct gkyl_array *fin_neut[app->num_neut_species];
    const struct gkyl_array *fin[app->num_species];
    for (int i=0; i<app->num_species; ++i) 
      fin[i] = app->species[i].f;
    for (int i=0; i<app->num_neut_species; ++i)
      fin_neut[i] = app->neut_species[i].f;

    gk_species_radiation_emissivity(app, gk_s, &gk_s->rad, fin, fin_neut);
    for (int i=0; i<gk_s->rad.num_cross_collisions; i++) {
      // copy data from device to host before writing it out
      if (app->use_gpu) {
        gkyl_array_copy(gk_s->rad.emissivity_host[i], gk_s->rad.emissivity[i]);
      }
      // Construct the file handles for vparallel and mu drag
      const char *fmt_emissivity = "%s-%s_radiation_emissivity_%s_%d.gkyl";  
      if (gk_s->rad.is_neut_species[i]) {
        int sz_emissivity = gkyl_calc_strlen(fmt_emissivity, app->name, gk_s->info.name,
          app->neut_species[gk_s->rad.collide_with_idx[i]].info.name, frame);
        char fileNm_emissivity[sz_emissivity+1]; // ensures no buffer overflow
        snprintf(fileNm_emissivity, sizeof fileNm_emissivity, fmt_emissivity, app->name,
          gk_s->info.name, app->neut_species[gk_s->rad.collide_with_idx[i]].info.name, frame);
        struct timespec wtm = gkyl_wall_clock();
        gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gk_s->rad.emissivity_host[i], fileNm_emissivity);
        app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      } else {
        int sz_emissivity = gkyl_calc_strlen(fmt_emissivity, app->name, gk_s->info.name,
          app->species[gk_s->rad.collide_with_idx[i]].info.name, frame);
        char fileNm_emissivity[sz_emissivity+1]; // ensures no buffer overflow
        snprintf(fileNm_emissivity, sizeof fileNm_emissivity, fmt_emissivity, app->name,
          gk_s->info.name, app->species[gk_s->rad.collide_with_idx[i]].info.name, frame);
        struct timespec wtm = gkyl_wall_clock();
        gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gk_s->rad.emissivity_host[i], fileNm_emissivity);
        app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      }  
      app->stat.nio += 1;
    }

    gyrokinetic_array_meta_release(mt);   
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_calc_species_rad_integrated_mom(gkyl_gyrokinetic_app *app, int sidx, double tm)
{
  struct gk_species *gk_s = &app->species[sidx];

  if (gk_s->rad.radiation_id == GKYL_GK_RADIATION) {
    struct timespec wst = gkyl_wall_clock();

    int vdim = app->vdim;
    double avals_global[2+vdim];

    // Compute radiation drag coefficients
    const struct gkyl_array *fin_neut[app->num_neut_species];
    const struct gkyl_array *fin[app->num_species];
    for (int i=0; i<app->num_species; ++i) 
      fin[i] = app->species[i].f;
    for (int i=0; i<app->num_neut_species; ++i)
      fin_neut[i] = app->neut_species[i].f;
    gk_species_radiation_moms(app, gk_s, &gk_s->rad, fin, fin_neut);
    gk_species_radiation_integrated_moms(app, gk_s, &gk_s->rad, fin, fin_neut);
  
    // reduce to compute sum over whole domain, append to diagnostics
    gkyl_array_reduce_range(gk_s->rad.red_integ_diag, gk_s->rad.integ_moms.marr, GKYL_SUM, &app->local);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, 
      gk_s->rad.red_integ_diag, gk_s->rad.red_integ_diag_global);
    if (app->use_gpu) {
      gkyl_cu_memcpy(avals_global, gk_s->rad.red_integ_diag_global, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      memcpy(avals_global, gk_s->rad.red_integ_diag_global, sizeof(double[2+vdim]));
    }
    gkyl_dynvec_append(gk_s->rad.integ_diag, tm, avals_global);

    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

void
gkyl_gyrokinetic_app_write_species_rad_integrated_mom(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_species *gk_s = &app->species[sidx];

  if (gk_s->rad.radiation_id == GKYL_GK_RADIATION) {
    struct timespec wst = gkyl_wall_clock();
    // Write from rank 0
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out integrated diagnostic moments
      const char *fmt = "%s-%s_radiation_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gk_s->info.name, "integrated_moms");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_s->info.name, "integrated_moms");

      struct timespec wtm = gkyl_wall_clock();
      if (gk_s->rad.is_first_integ_write_call) {
        gkyl_dynvec_write(gk_s->rad.integ_diag, fileNm);
        gk_s->rad.is_first_integ_write_call = false;
      }
      else {
        gkyl_dynvec_awrite(gk_s->rad.integ_diag, fileNm);
      }
      app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.nio += 1;
    }
    gkyl_dynvec_clear(gk_s->rad.integ_diag);
    app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
    app->stat.ndiag += 1;
  }
}

//
// ............. Ionization outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_iz_react(gkyl_gyrokinetic_app* app, int sidx, int ridx, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gk_species *gk_s = &app->species[sidx];

  struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  // Compute reaction rate
  const struct gkyl_array *fin[app->num_species];
  const struct gkyl_array *fin_neut[app->num_neut_species];
  for (int i=0; i<app->num_species; ++i) 
    fin[i] = app->species[i].f;
  for (int i=0; i<app->num_neut_species; ++i)
    fin_neut[i] = app->neut_species[i].f;
  gk_species_react_cross_moms(app, gk_s, &gk_s->react, fin[sidx], fin, fin_neut);

  const char *fmt = "%s-%s_%s_%s_iz_react_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gk_s->info.name,
    gk_s->react.react_type[ridx].ion_nm, gk_s->react.react_type[ridx].donor_nm, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_s->info.name,
    gk_s->react.react_type[ridx].ion_nm, gk_s->react.react_type[ridx].donor_nm, frame);
  
  if (app->use_gpu) {
    gkyl_array_copy(gk_s->react.coeff_react_host[ridx], gk_s->react.coeff_react[ridx]);
  }
  struct timespec wtm = gkyl_wall_clock();
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gk_s->react.coeff_react_host[ridx], fileNm);
  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.nio += 1;

  gyrokinetic_array_meta_release(mt); 
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_gyrokinetic_app_write_species_iz_react_neut(gkyl_gyrokinetic_app* app, int sidx, int ridx, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  struct gk_species *gk_s = &app->species[sidx];

  // Compute reaction rate
  const struct gkyl_array *fin[app->num_species];
  const struct gkyl_array *fin_neut[app->num_neut_species];
  for (int i=0; i<app->num_species; ++i) 
    fin[i] = app->species[i].f;
  for (int i=0; i<app->num_neut_species; ++i)
    fin_neut[i] = app->neut_species[i].f;
  gk_species_react_cross_moms(app, gk_s, &gk_s->react_neut, fin[sidx], fin, fin_neut);

  const char *fmt = "%s-%s_%s_%s_iz_react_neut_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gk_s->info.name,
    gk_s->react_neut.react_type[ridx].ion_nm, gk_s->react_neut.react_type[ridx].donor_nm, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_s->info.name,
    gk_s->react_neut.react_type[ridx].ion_nm, gk_s->react_neut.react_type[ridx].donor_nm, frame);
  
  if (app->use_gpu) {
    gkyl_array_copy(gk_s->react_neut.coeff_react_host[ridx], gk_s->react_neut.coeff_react[ridx]);
  }
  struct timespec wtm = gkyl_wall_clock();
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gk_s->react_neut.coeff_react_host[ridx], fileNm);
  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.nio += 1;

  gyrokinetic_array_meta_release(mt); 
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

//
// ............. Recombination outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_recomb_react(gkyl_gyrokinetic_app* app, int sidx, int ridx, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  struct gk_species *gk_s = &app->species[sidx];

  // Compute reaction rate
  const struct gkyl_array *fin[app->num_species];
  const struct gkyl_array *fin_neut[app->num_neut_species];
  for (int i=0; i<app->num_species; ++i) 
    fin[i] = app->species[i].f;
  for (int i=0; i<app->num_neut_species; ++i) 
    fin_neut[i] = app->neut_species[i].f;  
  gk_species_react_cross_moms(app, gk_s, &gk_s->react, fin[sidx], fin, fin_neut);

  const char *fmt = "%s-%s_%s_%s_recomb_react_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gk_s->info.name,
    gk_s->react.react_type[ridx].ion_nm, gk_s->react.react_type[ridx].recvr_nm, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_s->info.name,
    gk_s->react.react_type[ridx].ion_nm, gk_s->react.react_type[ridx].recvr_nm, frame);
  
  if (app->use_gpu) {
    gkyl_array_copy(gk_s->react.coeff_react_host[ridx], gk_s->react.coeff_react[ridx]);
  }
  struct timespec wtm = gkyl_wall_clock();
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gk_s->react.coeff_react_host[ridx], fileNm);
  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.nio += 1;

  gyrokinetic_array_meta_release(mt); 
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_gyrokinetic_app_write_species_recomb_react_neut(gkyl_gyrokinetic_app* app, int sidx, int ridx, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  struct gk_species *gk_s = &app->species[sidx];

  // Compute reaction rate
  const struct gkyl_array *fin[app->num_species];
  const struct gkyl_array *fin_neut[app->num_neut_species];
  for (int i=0; i<app->num_species; ++i) 
    fin[i] = app->species[i].f;
  for (int i=0; i<app->num_neut_species; ++i) 
    fin_neut[i] = app->neut_species[i].f;  
  gk_species_react_cross_moms(app, gk_s, &gk_s->react_neut, fin[sidx], fin, fin_neut);

  const char *fmt = "%s-%s_%s_%s_recomb_react_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gk_s->info.name,
    gk_s->react_neut.react_type[ridx].ion_nm, gk_s->react_neut.react_type[ridx].recvr_nm, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_s->info.name,
    gk_s->react_neut.react_type[ridx].ion_nm, gk_s->react_neut.react_type[ridx].recvr_nm, frame);
  
  if (app->use_gpu) {
    gkyl_array_copy(gk_s->react_neut.coeff_react_host[ridx], gk_s->react_neut.coeff_react[ridx]);
  }
  struct timespec wtm = gkyl_wall_clock();
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gk_s->react_neut.coeff_react_host[ridx], fileNm);
  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.nio += 1;

  gyrokinetic_array_meta_release(mt); 
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

//
// ............. Charge exchange outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_cx_react_neut(gkyl_gyrokinetic_app* app, int sidx, int ridx, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_array_meta *mt = gyrokinetic_array_meta_new( (struct gyrokinetic_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  struct gk_species *gk_s = &app->species[sidx];

  // Compute reaction rate
  const struct gkyl_array *fin[app->num_species];
  const struct gkyl_array *fin_neut[app->num_neut_species];
  for (int i=0; i<app->num_species; ++i) 
    fin[i] = app->species[i].f;
  for (int i=0; i<app->num_neut_species; ++i) 
    fin_neut[i] = app->neut_species[i].f;  
  gk_species_react_cross_moms(app, gk_s, &gk_s->react_neut, fin[sidx], fin, fin_neut);

  const char *fmt = "%s-%s_%s_cx_react_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gk_s->info.name,
    gk_s->react_neut.react_type[ridx].partner_nm, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gk_s->info.name,
    gk_s->react_neut.react_type[ridx].partner_nm, frame);
  
  if (app->use_gpu) {
    gkyl_array_copy(gk_s->react_neut.coeff_react_host[ridx], gk_s->react_neut.coeff_react[ridx]);
  }
  struct timespec wtm = gkyl_wall_clock();
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gk_s->react_neut.coeff_react_host[ridx], fileNm);
  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.nio += 1;

  gyrokinetic_array_meta_release(mt); 
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

//
// ............. Functions that group several outputs for a single species ............... //
// 
void
gkyl_gyrokinetic_app_write_species_phase(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_app_write_species(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_source(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_rad_drag(app, sidx, tm, frame);
}

void
gkyl_gyrokinetic_app_write_neut_species_phase(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_app_write_neut_species(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_neut_species_source(app, sidx, tm, frame);
}

void
gkyl_gyrokinetic_app_write_species_conf(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_app_write_species_mom(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_source_mom(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_lbo_mom(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_rad_emissivity(app, sidx, tm, frame);

  struct gk_species *gk_s = &app->species[sidx];
  for (int j=0; j<gk_s->react.num_react; ++j) {
    if ((gk_s->react.react_id[j] == GKYL_REACT_IZ) 
      && (gk_s->react.type_self[j] == GKYL_SELF_ELC)) {
      gkyl_gyrokinetic_app_write_species_iz_react(app, sidx, j, tm, frame);
    }
    if ((gk_s->react.react_id[j] == GKYL_REACT_RECOMB) 
      && (gk_s->react.type_self[j] == GKYL_SELF_ELC)) {
      gkyl_gyrokinetic_app_write_species_recomb_react(app, sidx, j, tm, frame);
    }
  }
  for (int j=0; j<gk_s->react_neut.num_react; ++j) {
    if ((gk_s->react_neut.react_id[j] == GKYL_REACT_IZ) 
      && (gk_s->react_neut.type_self[j] == GKYL_SELF_ELC)) {
      gkyl_gyrokinetic_app_write_species_iz_react_neut(app, sidx, j, tm, frame);
    }
    if ((gk_s->react_neut.react_id[j] == GKYL_REACT_RECOMB) 
      && (gk_s->react_neut.type_self[j] == GKYL_SELF_ELC)) {
      gkyl_gyrokinetic_app_write_species_recomb_react_neut(app, sidx, j, tm, frame);
    }
    if ((gk_s->react_neut.react_id[j] == GKYL_REACT_CX)
      && (gk_s->react_neut.type_self[j] == GKYL_SELF_ION)) {
      gkyl_gyrokinetic_app_write_species_cx_react_neut(app, sidx, j, tm, frame);
    }
  }
}

void
gkyl_gyrokinetic_app_write_neut_species_conf(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_app_write_neut_species_mom(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_neut_species_source_mom(app, sidx, tm, frame);
}

//
// ............. Functions that group several species outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_mom(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_app_write_species_mom(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_species_source_mom(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_species_lbo_mom(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_species_rad_emissivity(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_app_write_neut_species_mom(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_neut_species_source_mom(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_app_calc_integrated_mom(gkyl_gyrokinetic_app* app, double tm)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_app_calc_species_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_app_calc_species_source_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_app_calc_species_rad_integrated_mom(app, i, tm);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_app_calc_neut_species_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_app_calc_neut_species_source_integrated_mom(app, i, tm);
  }
}

void
gkyl_gyrokinetic_app_write_integrated_mom(gkyl_gyrokinetic_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_app_write_species_integrated_mom(app, i);
    gkyl_gyrokinetic_app_write_species_source_integrated_mom(app, i);
    gkyl_gyrokinetic_app_write_species_bgk_max_corr_status(app, i);
    gkyl_gyrokinetic_app_write_species_rad_integrated_mom(app, i);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_app_write_neut_species_integrated_mom(app, i);
    gkyl_gyrokinetic_app_write_neut_species_source_integrated_mom(app, i);
  }
}

void
gkyl_gyrokinetic_app_write_conf(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  gkyl_gyrokinetic_app_write_field(app, tm, frame);

  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_app_write_species_conf(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_app_write_neut_species_conf(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_app_write_phase(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_app_write_species_phase(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_app_write_neut_species_phase(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_app_write(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  gkyl_gyrokinetic_app_write_phase(app, tm, frame);

  gkyl_gyrokinetic_app_write_conf(app, tm, frame);
}
//
// ............. End of write functions ............... //
// 

// Take a forward Euler step with the suggested time-step dt. This may
// not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by
// stability. The actual time-step and dt_suggested are returned in
// the status object.
static void
forward_euler(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], struct gkyl_array *fout[], 
  const struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[], 
  struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;

  double dtmin = DBL_MAX;

  // Compute necessary moments and boundary corrections for collisions.
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].lbo.collision_id == GKYL_LBO_COLLISIONS) {
      gk_species_lbo_moms(app, &app->species[i], 
        &app->species[i].lbo, fin[i]);
    }
    if (app->species[i].bgk.collision_id == GKYL_BGK_COLLISIONS) {
      gk_species_bgk_moms(app, &app->species[i], 
        &app->species[i].bgk, fin[i]);
    }
  }

  // Compute necessary moments for cross-species collisions.
  // Needs to be done after self-collisions moments, so separate loop over species.
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gk_s = &app->species[i];

    if (gk_s->lbo.collision_id == GKYL_LBO_COLLISIONS) { 
      if (gk_s->lbo.num_cross_collisions) {
        gk_species_lbo_cross_moms(app, &app->species[i], 
          &gk_s->lbo, fin[i]);        
      }
    }
    if (gk_s->bgk.collision_id == GKYL_BGK_COLLISIONS) {
      if (gk_s->bgk.num_cross_collisions) {
        gk_species_bgk_cross_moms(app, &app->species[i], 
          &gk_s->bgk, fin[i]);        
      }
    }
    // Compute reaction rates (e.g., ionization, recombination, or charge exchange).
    if (gk_s->react.num_react) {
      gk_species_react_cross_moms(app, &app->species[i], 
        &gk_s->react, fin[i], fin, fin_neut);
    }
    if (gk_s->react_neut.num_react) {
      gk_species_react_cross_moms(app, &app->species[i], 
        &gk_s->react_neut, fin[i], fin, fin_neut);
    }
    // Compute necessary drag coefficients for radiation operator.
    if (gk_s->rad.radiation_id == GKYL_GK_RADIATION) {
      gk_species_radiation_moms(app, &app->species[i], 
        &gk_s->rad, fin, fin_neut);
    }
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    // Compute reaction cross moments (e.g., ionization, recombination, or charge exchange).
    if (app->neut_species[i].react_neut.num_react) {
      gk_neut_species_react_cross_moms(app, &app->neut_species[i], 
        &app->neut_species[i].react_neut, fin, fin_neut);
    }
  }

  // Compute RHS of Gyrokinetic equation.
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    double dt1 = gk_species_rhs(app, s, fin[i], fout[i]);
    dtmin = fmin(dtmin, dt1);

    // Compute and store (in the ghost cell of of out) the boundary fluxes.
    // NOTE: this overwrites ghost cells that may be used for sourcing.
    if (app->update_field && app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN)
      gk_species_bflux_rhs(app, s, &s->bflux, fin[i], fout[i]);
  }

  // Compute RHS of neutrals.
  for (int i=0; i<app->num_neut_species; ++i) {
    double dt1 = gk_neut_species_rhs(app, &app->neut_species[i], fin_neut[i], fout_neut[i]);
    dtmin = fmin(dtmin, dt1);
  }

  // Compute plasma source term.
  // Done here as the RHS update for all species should be complete before
  // in case we are using boundary fluxes as a component of our source function
  for (int i=0; i<app->num_species; ++i) {
    gk_species_source_rhs(app, &app->species[i], 
      &app->species[i].src, fin[i], fout[i]);
  }

  // Compute neutral source term.
  // Done here as the RHS update for all species should be complete before
  // in case we are using boundary fluxes as a component of our source function.
  for (int i=0; i<app->num_neut_species; ++i) {
    gk_neut_species_source_rhs(app, &app->neut_species[i], 
      &app->neut_species[i].src, fin_neut[i], fout_neut[i]);
  }

  double dt_max_rel_diff = 0.01;
  // Check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // Compute minimum time-step across all processors.
  double dtmin_local = dtmin, dtmin_global;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dtmin_local, &dtmin_global);
  dtmin = dtmin_global;
  
  // Don't take a time-step larger that input dt.
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // Complete update of distribution functions.
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout[i], dta), 1.0, fin[i]);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    if (!app->neut_species[i].info.is_static) {
      gkyl_array_accumulate(gkyl_array_scale(fout_neut[i], dta), 1.0, fin_neut[i]);
    }
  }

}

// Take time-step using the RK3 method. Also sets the status object
// which has the actual and suggested dts used. These can be different
// from the actual time-step.
static struct gkyl_update_status
rk3(gkyl_gyrokinetic_app* app, double dt0)
{
  const struct gkyl_array *fin[app->num_species];
  struct gkyl_array *fout[app->num_species];
  const struct gkyl_array *fin_neut[app->num_neut_species];
  struct gkyl_array *fout_neut[app->num_neut_species];
  struct gkyl_update_status st = { .success = true };

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f;
          fout[i] = app->species[i].f1;
        }
        for (int i=0; i<app->num_neut_species; ++i) {
          fin_neut[i] = app->neut_species[i].f;
          if (!app->neut_species[i].info.is_static) {
            fout_neut[i] = app->neut_species[i].f1;
          }
        }

        forward_euler(app, tcurr, dt, fin, fout, fin_neut, fout_neut, &st);
        // Compute the fields and apply BCs.
        calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

        dt = st.dt_actual;
        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
        }
        for (int i=0; i<app->num_neut_species; ++i) {
          if (!app->neut_species[i].info.is_static) {
            fin_neut[i] = app->neut_species[i].f1;
            fout_neut[i] = app->neut_species[i].fnew;
          }
          else {
            fin_neut[i] = app->neut_species[i].f;
          }
        }

        forward_euler(app, tcurr+dt, dt, fin, fout, fin_neut, fout_neut, &st);

        if (st.dt_actual < dt) {

          // Recalculate the field.
          for (int i=0; i<app->num_species; ++i)
            fin[i] = app->species[i].f;
          calc_field(app, tcurr, fin);

          // collect stats
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_2_fail += 1;

          dt = st.dt_actual;
          state = RK_STAGE_1; // restart from stage 1

        } 
        else {
          for (int i=0; i<app->num_species; ++i) {
            array_combine(app->species[i].f1,
              3.0/4.0, app->species[i].f, 1.0/4.0, app->species[i].fnew, &app->species[i].local_ext);
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            if (!app->neut_species[i].info.is_static) {
              array_combine(app->neut_species[i].f1,
                3.0/4.0, app->neut_species[i].f, 1.0/4.0, app->neut_species[i].fnew, &app->neut_species[i].local_ext);
            }
          }

          // Compute the fields and apply BCs.
          for (int i=0; i<app->num_species; ++i) {
            fout[i] = app->species[i].f1;
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            fout_neut[i] = app->neut_species[i].f1;
          }
          calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
        }
        for (int i=0; i<app->num_neut_species; ++i) {
          if (!app->neut_species[i].info.is_static) {
            fin_neut[i] = app->neut_species[i].f1;
            fout_neut[i] = app->neut_species[i].fnew;
          }
          else {
            fin_neut[i] = app->neut_species[i].f;
          }          
        }

        forward_euler(app, tcurr+dt/2, dt, fin, fout, fin_neut, fout_neut, &st);

        if (st.dt_actual < dt) {
          // Recalculate the field.
          for (int i=0; i<app->num_species; ++i)
            fin[i] = app->species[i].f;
          calc_field(app, tcurr, fin);

          // collect stats
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_3_fail += 1;

          dt = st.dt_actual;
          state = RK_STAGE_1; // restart from stage 1

          app->stat.nstage_2_fail += 1;
        }
        else {
          for (int i=0; i<app->num_species; ++i) {
            array_combine(app->species[i].f1,
              1.0/3.0, app->species[i].f, 2.0/3.0, app->species[i].fnew, &app->species[i].local_ext);
            gkyl_array_copy_range(app->species[i].f, app->species[i].f1, &app->species[i].local_ext);
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            if (!app->neut_species[i].info.is_static) {
              array_combine(app->neut_species[i].f1,
                1.0/3.0, app->neut_species[i].f, 2.0/3.0, app->neut_species[i].fnew, &app->neut_species[i].local_ext);
              gkyl_array_copy_range(app->neut_species[i].f, app->neut_species[i].f1, &app->neut_species[i].local_ext);
            }
          }

          if (app->enforce_positivity) {
            // Apply positivity shift if requested.
            int elc_idx = -1;
            gkyl_array_clear(app->ps_delta_m0_ions, 0.0);
            for (int i=0; i<app->num_species; ++i) {
              struct gk_species *gks = &app->species[i];

              // Copy f so we can calculate the moments of the change later. 
              gkyl_array_set(gks->fnew, -1.0, gks->f);

              // Shift each species.
              gkyl_positivity_shift_gyrokinetic_advance(gks->pos_shift_op, &app->local, &gks->local,
                gks->f, gks->m0.marr, gks->ps_delta_m0);

              // Accumulate the shift density of all ions:
              if (gks->info.charge > 0.0)
                gkyl_array_accumulate(app->ps_delta_m0_ions, 1.0, gks->ps_delta_m0);
              else if (gks->info.charge < 0.0) 
                elc_idx = i;
            }

            // Rescale each species to enforce quasineutrality.
            for (int i=0; i<app->num_species; ++i) {
              struct gk_species *gks = &app->species[i];
              if (gks->info.charge > 0.0) {
                struct gk_species *gkelc = &app->species[elc_idx];
                gkyl_positivity_shift_gyrokinetic_quasineutrality_scale(gks->pos_shift_op, &app->local, &gks->local,
                  gks->ps_delta_m0, app->ps_delta_m0_ions, gkelc->ps_delta_m0, gks->m0.marr, gks->f);
              }
              else {
                gkyl_positivity_shift_gyrokinetic_quasineutrality_scale(gks->pos_shift_op, &app->local, &gks->local,
                  gks->ps_delta_m0, gks->ps_delta_m0, app->ps_delta_m0_ions, gks->m0.marr, gks->f);
              }

              gkyl_array_accumulate(gks->fnew, 1.0, gks->f);
            }
          }

          // Compute the fields and apply BCs
          for (int i=0; i<app->num_species; ++i) {
            fout[i] = app->species[i].f;
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            fout_neut[i] = app->neut_species[i].f;
          }
          calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

          state = RK_COMPLETE;
        }
        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }

  return st;
}

struct gkyl_update_status
gkyl_gyrokinetic_update(gkyl_gyrokinetic_app* app, double dt)
{
  app->stat.nup += 1;
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_update_status status = rk3(app, dt);
  app->tcurr += status.dt_actual;

  app->stat.total_tm += gkyl_time_diff_now_sec(wst);
  // Check for any CUDA errors during time step
  if (app->use_gpu)
    checkCuda(cudaGetLastError());
  return status;
}

struct gkyl_gyrokinetic_stat
gkyl_gyrokinetic_app_stat(gkyl_gyrokinetic_app* app)
{
  gk_species_tm(app);
  gk_species_coll_tm(app);
  return app->stat;
}

void
gkyl_gyrokinetic_app_species_ktm_rhs(gkyl_gyrokinetic_app* app, int update_vol_term)
{
  for (int i=0; i<app->num_species; ++i) {

    struct gk_species *species = &app->species[i];

    const struct gkyl_array *fin = species->f;
    struct gkyl_array *rhs = species->f1;

    gkyl_array_clear(rhs, 0.0);
    gkyl_dg_updater_gyrokinetic_advance(species->slvr, &species->local, 
      fin, species->cflrate, rhs); 
  }
}

static void
range_stat_write(gkyl_gyrokinetic_app* app, const char *nm, const struct gkyl_range *r, FILE *fp)
{
  gkyl_gyrokinetic_app_cout(app, fp, " %s_cells : [ ", nm);
  for (int i=0; i<r->ndim; ++i)
    gkyl_gyrokinetic_app_cout(app, fp, " %d, ", gkyl_range_shape(r, i));
  gkyl_gyrokinetic_app_cout(app, fp, " ],\n");
}

// ensure stats across processors are made consistent
static void
comm_reduce_app_stat(const gkyl_gyrokinetic_app* app,
  const struct gkyl_gyrokinetic_stat *local, struct gkyl_gyrokinetic_stat *global)
{
  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);
  if (comm_sz == 1) {
    memcpy(global, local, sizeof(struct gkyl_gyrokinetic_stat));
    return;
  }

  global->use_gpu = local->use_gpu;

  enum { NUP, NFEULER, NSTAGE_2_FAIL, NSTAGE_3_FAIL, L_END };
  int64_t l_red[] = {
    [NUP] = local->nup,
    [NFEULER] = local->nfeuler,
    [NSTAGE_2_FAIL] = local->nstage_2_fail,
    [NSTAGE_3_FAIL] = local->nstage_3_fail
  };

  int64_t l_red_global[L_END];
  gkyl_comm_allreduce_host(app->comm, GKYL_INT_64, GKYL_MAX, L_END, l_red, l_red_global);

  global->nup = l_red_global[NUP];
  global->nfeuler = l_red_global[NFEULER];
  global->nstage_2_fail = l_red_global[NSTAGE_2_FAIL];
  global->nstage_3_fail = l_red_global[NSTAGE_3_FAIL];  

  enum {
    TOTAL_TM, INIT_SPECIES_TM, SPECIES_RHS_TM, 
    SPECIES_COLL_MOM_TM, SPECIES_COL_TM, FIELD_RHS_TM, 
    SPECIES_OMEGA_CFL_TM, MOM_TM, DIAG_TM, IO_TM, SPECIES_BC_TM, 
    D_END
  };

  double d_red[D_END] = {
    [TOTAL_TM] = local->total_tm,
    [INIT_SPECIES_TM] = local->init_species_tm,
    [SPECIES_RHS_TM] = local->species_rhs_tm,
    [SPECIES_COLL_MOM_TM] = local->species_coll_mom_tm,
    [SPECIES_COL_TM] = local->species_coll_tm,
    [FIELD_RHS_TM] = local->field_rhs_tm,
    [SPECIES_OMEGA_CFL_TM] = local->species_omega_cfl_tm,
    [MOM_TM] = local->mom_tm,
    [DIAG_TM] = local->diag_tm,
    [IO_TM] = local->io_tm,
    [SPECIES_BC_TM] = local->species_bc_tm,
  };

  double d_red_global[D_END];
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, D_END, d_red, d_red_global);
  
  global->total_tm = d_red_global[TOTAL_TM];
  global->init_species_tm = d_red_global[INIT_SPECIES_TM];
  global->species_rhs_tm = d_red_global[SPECIES_RHS_TM];
  global->species_coll_mom_tm = d_red_global[SPECIES_COLL_MOM_TM];
  global->species_coll_tm = d_red_global[SPECIES_COL_TM];
  global->field_rhs_tm = d_red_global[FIELD_RHS_TM];
  global->species_omega_cfl_tm = d_red_global[SPECIES_OMEGA_CFL_TM];
  global->mom_tm = d_red_global[MOM_TM];
  global->diag_tm = d_red_global[DIAG_TM];
  global->io_tm = d_red_global[IO_TM];
  global->species_bc_tm = d_red_global[SPECIES_BC_TM];

  // misc data needing reduction

  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_2_dt_diff,
    global->stage_2_dt_diff);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_3_dt_diff,
    global->stage_3_dt_diff);

  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, GKYL_MAX_SPECIES, local->species_lbo_coll_drag_tm,
    global->species_lbo_coll_drag_tm);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, GKYL_MAX_SPECIES, local->species_lbo_coll_diff_tm,
    global->species_lbo_coll_diff_tm);
}

void
gkyl_gyrokinetic_app_stat_write(gkyl_gyrokinetic_app* app)
{
  const char *fmt = "%s-%s";
  int sz = gkyl_calc_strlen(fmt, app->name, "stat.json");
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, "stat.json");

  int num_ranks;
  gkyl_comm_get_size(app->comm, &num_ranks);

  char buff[70];
  time_t t = time(NULL);
  struct tm curr_tm = *localtime(&t);

  gk_species_coll_tm(app);
  gk_species_tm(app);

  struct gkyl_gyrokinetic_stat stat = { };
  comm_reduce_app_stat(app, &app->stat, &stat);
  
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  if (rank == 0) fp = fopen(fileNm, "a");

  gkyl_gyrokinetic_app_cout(app, fp, "{\n");

  if (strftime(buff, sizeof buff, "%c", &curr_tm))
    gkyl_gyrokinetic_app_cout(app, fp, " date : %s,\n", buff);

  gkyl_gyrokinetic_app_cout(app, fp, " use_gpu : %d,\n", stat.use_gpu);
  gkyl_gyrokinetic_app_cout(app, fp, " num_ranks : %d,\n", num_ranks); 
  
  for (int s=0; s<app->num_species; ++s)
    range_stat_write(app, app->species[s].info.name, &app->species[s].global, fp);
  
  gkyl_gyrokinetic_app_cout(app, fp, " nup : %ld,\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, fp, " nfeuler : %ld,\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, fp, " nstage_2_fail : %ld,\n", stat.nstage_2_fail);
  gkyl_gyrokinetic_app_cout(app, fp, " nstage_3_fail : %ld,\n", stat.nstage_3_fail);

  gkyl_gyrokinetic_app_cout(app, fp, " stage_2_dt_diff : [ %lg, %lg ],\n",
    stat.stage_2_dt_diff[0], stat.stage_2_dt_diff[1]);
  gkyl_gyrokinetic_app_cout(app, fp, " stage_3_dt_diff : [ %lg, %lg ],\n",
    stat.stage_3_dt_diff[0], stat.stage_3_dt_diff[1]);

  gkyl_gyrokinetic_app_cout(app, fp, " total_tm : %lg,\n", stat.total_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " init_species_tm : %lg,\n", stat.init_species_tm);

  gkyl_gyrokinetic_app_cout(app, fp, " species_rhs_tm : %lg,\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " field_rhs_tm : %lg,\n", stat.field_rhs_tm);

  for (int s=0; s<app->num_species; ++s) {
    gkyl_gyrokinetic_app_cout(app, fp, " species_coll_drag_tm[%d] : %lg,\n", s,
      stat.species_lbo_coll_drag_tm[s]);
    gkyl_gyrokinetic_app_cout(app, fp, " species_coll_diff_tm[%d] : %lg,\n", s,
      stat.species_lbo_coll_diff_tm[s]);
  }

  gkyl_gyrokinetic_app_cout(app, fp, " species_coll_mom_tm : %lg,\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_coll_tm : %lg,\n", stat.species_coll_tm);

  gkyl_gyrokinetic_app_cout(app, fp, " species_bc_tm : %lg,\n", stat.species_bc_tm);

  gkyl_gyrokinetic_app_cout(app, fp, " nmom : %ld,\n", stat.nmom);
  gkyl_gyrokinetic_app_cout(app, fp, " mom_tm : %lg\n", stat.mom_tm);

  gkyl_gyrokinetic_app_cout(app, fp, " ndiag : %ld,\n", stat.ndiag);
  gkyl_gyrokinetic_app_cout(app, fp, " diag_tm : %lg\n", stat.diag_tm);
  
  gkyl_gyrokinetic_app_cout(app, fp, " nspecies_omega_cfl : %ld,\n", stat.nspecies_omega_cfl);
  gkyl_gyrokinetic_app_cout(app, fp, " species_omega_cfl_tm : %lg\n", stat.species_omega_cfl_tm);

  gkyl_gyrokinetic_app_cout(app, fp, " nio : %ld,\n", stat.nio);
  gkyl_gyrokinetic_app_cout(app, fp, " io_tm : %lg\n", stat.io_tm);
  
  gkyl_gyrokinetic_app_cout(app, fp, "}\n");

  if (rank == 0)
    fclose(fp);  

}

static struct gkyl_app_restart_status
header_from_file(gkyl_gyrokinetic_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = { .io_status = 0 };
  
  FILE *fp = 0;
  with_file(fp, fname, "r") {
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    rstat.io_status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);

    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      if (hdr.etype != GKYL_DOUBLE)
        rstat.io_status = GKYL_ARRAY_RIO_DATA_MISMATCH;
    }

    struct gyrokinetic_output_meta meta =
      gyrokinetic_meta_from_mpack( &(struct gkyl_array_meta) {
          .meta = hdr.meta,
          .meta_sz = hdr.meta_size
        }
      );

    rstat.frame = meta.frame;
    rstat.stime = meta.stime;

    gkyl_grid_sub_array_header_release(&hdr);
  }
  
  return rstat;
}

//
// ............. Reading functions ............... //
// 
static void
gyrokinetic_app_geometry_read_and_copy(gkyl_gyrokinetic_app* app, struct gkyl_array *arr,
  struct gkyl_array *arr_host, char *varNm)
{
  cstr fileNm = cstr_from_fmt("%s-%s.gkyl", app->name, varNm);

  struct gkyl_app_restart_status rstat = header_from_file(app, fileNm.str);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, arr_host, fileNm.str);
    gkyl_array_copy(arr, arr_host);
  }
}

void
gkyl_gyrokinetic_app_read_geometry(gkyl_gyrokinetic_app* app)
{
  struct gkyl_array* arr_ho1 = mkarr(false,   app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array* arr_ho3 = mkarr(false, 3*app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array* arr_ho6 = mkarr(false, 6*app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array* arr_ho9 = mkarr(false, 9*app->confBasis.num_basis, app->local_ext.volume);

  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->mc2p        , arr_ho3, "mapc2p");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->bmag        , arr_ho1, "bmag");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->g_ij        , arr_ho6, "g_ij");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->dxdz        , arr_ho9, "dxdz");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->dzdx        , arr_ho9, "dzdx");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->normals     , arr_ho9, "normals");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->jacobgeo    , arr_ho1, "jacobgeo");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->jacobgeo_inv, arr_ho1, "jacobgeo_inv");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->gij         , arr_ho6, "gij");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->b_i         , arr_ho3, "b_i");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->bcart       , arr_ho3, "bcart");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->cmag        , arr_ho1, "cmag");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->jacobtot    , arr_ho1, "jacobtot");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->jacobtot_inv, arr_ho1, "jacobtot_inv");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->bmag_inv    , arr_ho1, "bmag_inv");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->bmag_inv_sq , arr_ho1, "bmag_inv_sq");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->gxxj        , arr_ho1, "gxxj");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->gxyj        , arr_ho1, "gxyj");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->gyyj        , arr_ho1, "gyyj");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->gxzj        , arr_ho1, "gxzj");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->eps2        , arr_ho1, "eps2");

  gkyl_array_release(arr_ho1);
  gkyl_array_release(arr_ho3);
  gkyl_array_release(arr_ho6);
  gkyl_array_release(arr_ho9);
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_app_from_file_field(gkyl_gyrokinetic_app *app, const char *fname)
{
  if (app->update_field != 1)
    return (struct gkyl_app_restart_status) {
      .io_status = GKYL_ARRAY_RIO_SUCCESS,
      .frame = 0,
      .stime = 0.0
    };

  struct gkyl_app_restart_status rstat = header_from_file(app, fname);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, app->field->phi_host, fname);
    if (app->use_gpu)
      gkyl_array_copy(app->field->phi_smooth, app->field->phi_host);
  }
  
  return rstat;
}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_app_from_file_species(gkyl_gyrokinetic_app *app, int sidx,
  const char *fname)
{
  struct gkyl_app_restart_status rstat = header_from_file(app, fname);

  struct gk_species *gk_s = &app->species[sidx];
  
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(gk_s->comm, &gk_s->grid, &gk_s->local, gk_s->f_host, fname);
    if (app->use_gpu)
      gkyl_array_copy(gk_s->f, gk_s->f_host);
    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      gk_species_source_calc(app, gk_s, &gk_s->src, 0.0);
    }
  }

  return rstat;
}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_app_from_file_neut_species(gkyl_gyrokinetic_app *app, int sidx,
  const char *fname)
{
  struct gkyl_app_restart_status rstat = header_from_file(app, fname);

  struct gk_neut_species *gk_ns = &app->neut_species[sidx];
  
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(gk_ns->comm, &gk_ns->grid, &gk_ns->local, gk_ns->f_host, fname);
    if (app->use_gpu)
      gkyl_array_copy(gk_ns->f, gk_ns->f_host);
    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      gk_neut_species_source_calc(app, gk_ns, &gk_ns->src, 0.0);
    }
  }

  return rstat;
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_app_from_frame_field(gkyl_gyrokinetic_app *app, int frame)
{
  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "field", frame);
  struct gkyl_app_restart_status rstat = gkyl_gyrokinetic_app_from_file_field(app, fileNm.str);
  app->field->is_first_energy_write_call = false; // Append to existing diagnostic.
  cstr_drop(&fileNm);
  
  return rstat;
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_app_from_frame_species(gkyl_gyrokinetic_app *app, int sidx, int frame)
{
  struct gk_species *gk_s = &app->species[sidx];

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, gk_s->info.name, frame);
  struct gkyl_app_restart_status rstat = gkyl_gyrokinetic_app_from_file_species(app, sidx, fileNm.str);
  cstr_drop(&fileNm);

  // Append to existing integrated diagnostics.
  gk_s->is_first_integ_write_call = false;
  if (app->enforce_positivity)
    gk_s->is_first_ps_integ_write_call = false;
  if (gk_s->rad.radiation_id == GKYL_GK_RADIATION)
    gk_s->rad.is_first_integ_write_call = false;
  if (gk_s->src.source_id)
    gk_s->src.is_first_integ_write_call = false;
  if (gk_s->bgk.collision_id == GKYL_BGK_COLLISIONS)
    gk_s->bgk.is_first_corr_status_write_call = false;

  return rstat;
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_app_from_frame_neut_species(gkyl_gyrokinetic_app *app, int sidx, int frame)
{
  struct gk_neut_species *gk_ns = &app->neut_species[sidx];

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, gk_ns->info.name, frame);
  struct gkyl_app_restart_status rstat = gkyl_gyrokinetic_app_from_file_neut_species(app, sidx, fileNm.str);
  gk_ns->is_first_integ_write_call = false; // append to existing diagnostic
  cstr_drop(&fileNm);
  
  // Append to existing integrated diagnostics.
  gk_ns->is_first_integ_write_call = false;
  if (gk_ns->src.source_id)
    gk_ns->src.is_first_integ_write_call = false;

  return rstat;
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_app_read_from_frame(gkyl_gyrokinetic_app *app, int frame)
{
  struct gkyl_app_restart_status rstat;
  for (int i=0; i<app->num_neut_species; i++) {
    int neut_frame = frame;
    if (app->neut_species[i].info.is_static) {
      neut_frame = 0;
    }
    rstat = gkyl_gyrokinetic_app_from_frame_neut_species(app, i, neut_frame);
  }
  for (int i=0; i<app->num_species; i++) {
    rstat = gkyl_gyrokinetic_app_from_frame_species(app, i, frame);
  }
  
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    // Compute the fields and apply BCs.
    struct gkyl_array *distf[app->num_species];
    struct gkyl_array *distf_neut[app->num_neut_species];
    for (int i=0; i<app->num_species; ++i) {
      distf[i] = app->species[i].f;
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      distf_neut[i] = app->neut_species[i].f;
    }
    if (app->update_field && app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
      for (int i=0; i<app->num_species; ++i) {
        struct gk_species *s = &app->species[i];

        // Compute advection speeds so we can compute the initial boundary flux.
        gkyl_dg_calc_gyrokinetic_vars_alpha_surf(s->calc_gk_vars, 
          &app->local, &s->local, &s->local_ext, app->field->phi_smooth,
          s->alpha_surf, s->sgn_alpha_surf, s->const_sgn_alpha);

        // Compute and store (in the ghost cell of of out) the boundary fluxes.
        // NOTE: this overwrites ghost cells that may be used for sourcing.
        gk_species_bflux_rhs(app, s, &s->bflux, distf[i], distf[i]);
      }
    }
    calc_field_and_apply_bc(app, rstat.stime, distf, distf_neut);
  }

  app->field->is_first_energy_write_call = false; // Append to existing diagnostic.

  return rstat;
}

// private function to handle variable argument list for printing
static void
v_gk_app_cout(const gkyl_gyrokinetic_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank, r = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp)
    vfprintf(fp, fmt, argp);
}

void
gkyl_gyrokinetic_app_cout(const gkyl_gyrokinetic_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_gk_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_gyrokinetic_app_release(gkyl_gyrokinetic_app* app)
{
  if (app->enforce_positivity) {
    gkyl_array_release(app->ps_delta_m0_ions);
  }

  gkyl_gk_geometry_release(app->gk_geom);

  gk_field_release(app, app->field);

  for (int i=0; i<app->num_species; ++i)
    gk_species_release(app, &app->species[i]);
  if (app->num_species > 0)
    gkyl_free(app->species);

  for (int i=0; i<app->num_neut_species; ++i)
    gk_neut_species_release(app, &app->neut_species[i]);
  if (app->num_neut_species > 0)
    gkyl_free(app->neut_species);

  gkyl_comm_release(app->comm);
  gkyl_rect_decomp_release(app->decomp);

  if (app->use_gpu) {
    gkyl_cu_free(app->basis_on_dev.basis);
    gkyl_cu_free(app->basis_on_dev.confBasis);
  }

  gkyl_free(app);
}
