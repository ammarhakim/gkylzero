#include <stdarg.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_basis.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>
#include <gkyl_null_comm.h>

#include <gkyl_vlasov_poisson_priv.h>
#include <gkyl_app_priv.h>

#include <mpack.h>

// returned gkyl_array_meta must be freed using vlasov_poisson_array_meta_release
static struct gkyl_array_meta*
vlasov_poisson_array_meta_new(struct vlasov_poisson_output_meta meta)
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
vlasov_poisson_array_meta_release(struct gkyl_array_meta *mt)
{
  if (!mt) return;
  MPACK_FREE(mt->meta);
  gkyl_free(mt);
}

static struct vlasov_poisson_output_meta
vlasov_poisson_meta_from_mpack(struct gkyl_array_meta *mt)
{
  struct vlasov_poisson_output_meta meta = { .frame = 0, .stime = 0.0 };

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

gkyl_vlasov_poisson_app*
gkyl_vlasov_poisson_app_new(struct gkyl_vp *vp)
{
  disable_denorm_float();

  assert(vp->num_species <= GKYL_MAX_SPECIES);

  gkyl_vlasov_poisson_app *app = gkyl_malloc(sizeof(gkyl_vlasov_poisson_app));

  int cdim = app->cdim = vp->cdim;
  int vdim = app->vdim = vp->vdim;
  int pdim = cdim+vdim;
  int poly_order = app->poly_order = vp->poly_order;
  int ns = app->num_species = vp->num_species;

  double cfl_frac = vp->cfl_frac == 0 ? 1.0 : vp->cfl_frac;
  app->cfl = cfl_frac;

#ifdef GKYL_HAVE_CUDA
  app->use_gpu = vp->use_gpu;
#else
  app->use_gpu = false; // Can't use GPUs if we don't have them!
#endif

  app->num_periodic_dir = vp->num_periodic_dir;
  for (int d=0; d<cdim; ++d)
    app->periodic_dirs[d] = vp->periodic_dirs[d];

  strcpy(app->name, vp->name);
  app->tcurr = 0.0; // reset on init

  if (app->use_gpu) {
    // Allocate device basis if we are using GPUs.
    app->basis_on_dev.basis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    app->basis_on_dev.confBasis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  } else {
    app->basis_on_dev.basis = &app->basis;
    app->basis_on_dev.confBasis = &app->confBasis;
  }

  // basis functions
  switch (vp->basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);
      if (poly_order == 1) {
        /* Force hybrid basis (p=2 in velocity space). */
        gkyl_cart_modal_hybrid(&app->basis, cdim, vdim);
        gkyl_cart_modal_serendip(&app->velBasis, vdim, 2);
      } else {
        gkyl_cart_modal_serendip(&app->basis, pdim, poly_order);
        gkyl_cart_modal_serendip(&app->velBasis, vdim, poly_order);
      } 

      if (app->use_gpu) {
        gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.confBasis, cdim, poly_order);
        if (poly_order == 1) {
          /* Force hybrid basis (p=2 in velocity space). */
          gkyl_cart_modal_hybrid_cu_dev(app->basis_on_dev.basis, cdim, vdim); 
        } else {
          gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.basis, pdim, poly_order);
        } 
      }
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      gkyl_cart_modal_tensor(&app->confBasis, cdim, poly_order);
      gkyl_cart_modal_tensor(&app->velBasis, vdim, poly_order);
      gkyl_cart_modal_tensor(&app->basis, pdim, poly_order);
      break;

    default:
      assert(false);
      break;
  }

  gkyl_rect_grid_init(&app->grid, cdim, vp->lower, vp->upper, vp->cells);

  int ghost[] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&app->grid, ghost, &app->global_ext, &app->global);

  if (vp->has_low_inp) {
    // Create local and local_ext from user-supplied local range.
    gkyl_create_ranges(&vp->low_inp.local_range, ghost, &app->local_ext, &app->local);
    
    if (vp->low_inp.comm)
      app->comm = gkyl_comm_acquire(vp->low_inp.comm);
    else {
      int cuts[3] = { 1, 1, 1 };
      struct gkyl_rect_decomp *rect_decomp =
        gkyl_rect_decomp_new_from_cuts(cdim, cuts, &app->global);
      
      app->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
          .decomp = rect_decomp,
          .use_gpu = app->use_gpu
        }
      );

      gkyl_rect_decomp_release(rect_decomp);
    }
  }
  else {
    // Global and local ranges are same, and so just copy.
    memcpy(&app->local, &app->global, sizeof(struct gkyl_range));
    memcpy(&app->local_ext, &app->global_ext, sizeof(struct gkyl_range));

    int cuts[3] = { 1, 1, 1 };
    struct gkyl_rect_decomp *rect_decomp =
      gkyl_rect_decomp_new_from_cuts(cdim, cuts, &app->global);
    
    app->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = rect_decomp,
        .use_gpu = app->use_gpu
      }
    );
    
    gkyl_rect_decomp_release(rect_decomp);
  }

  // Local skin and ghost ranges for configuration space fields.
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&app->lower_skin[dir], &app->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &app->local_ext, ghost); 
    gkyl_skin_ghost_ranges(&app->upper_skin[dir], &app->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &app->local_ext, ghost);
  }

  // Allocate space to store species objects.
  app->species = ns>0 ? gkyl_malloc(sizeof(struct vp_species[ns])) : 0;
  for (int i=0; i<ns; ++i)
    app->species[i] = (struct vp_species) { };

  // Set info for each species: this needs to be done here as we need
  // to access species name from vp_species_init
  for (int i=0; i<ns; ++i)
    app->species[i].info = vp->species[i];

  // Initialize the fields.
  app->has_field = !vp->skip_field; // Note inversion of truth.
  if (app->has_field)
    app->field = vp_field_new(vp, app);

  // Initialize each species
  for (int i=0; i<ns; ++i) 
    vp_species_init(vp, app, &app->species[i]);

  // initialize species wall emission terms: these rely
  // on other species which must be allocated in the previous step
  for (int i=0; i<ns; ++i) {
    if (app->species[i].emit_lo)
      vp_species_emission_cross_init(app, &app->species[i], &app->species[i].bc_emission_lo);
    if (app->species[i].emit_up)
      vp_species_emission_cross_init(app, &app->species[i], &app->species[i].bc_emission_up);
  }

  // Initialize each species cross-collisions terms: this has to be done here
  // as need pointers to colliding species' collision objects
  // allocated in vp_species_init.
  for (int i=0; i<ns; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      if (app->species[i].lbo.num_cross_collisions) {
        vp_species_lbo_cross_init(app, &app->species[i], &app->species[i].lbo);
      }
    }
  }

  // Initialize each species source terms
  for (int i=0; i<ns; ++i) {
    if (app->species[i].source_id) {
      vp_species_source_init(app, &app->species[i], &app->species[i].src);
    }
  }

  // Initialize stat object.
  app->stat = (struct gkyl_vlasov_poisson_stat) {
    .use_gpu = app->use_gpu,
    .stage_2_dt_diff = { DBL_MAX, 0.0 },
    .stage_3_dt_diff = { DBL_MAX, 0.0 },
  };

  return app;
}

// Compute fields.
static void
calc_field(gkyl_vlasov_poisson_app* app, double tcurr, const struct gkyl_array *fin[])
{
  if (app->has_field) {
    // Compute electrostatic potential from Poisson's equation.
    vp_field_accumulate_rho_c(app, app->field, fin);
  
    // Solve the field equation.
    vp_field_rhs(app, app->field);

    if (app->field->ext_em_evolve) {
      // Compute the external potentials.
      vp_field_calc_ext_em(app, app->field, tcurr);
    }

  }
}

// Compute fields and apply BCs.
static void
calc_field_and_apply_bc(gkyl_vlasov_poisson_app* app, double tcurr, struct gkyl_array *distf[])
{

  // Compute the field.
  calc_field(app, tcurr, (const struct gkyl_array **) distf);

  // Apply boundary conditions.
  for (int i=0; i<app->num_species; ++i) {
    vp_species_apply_bc(app, &app->species[i], distf[i], tcurr);
  }

}

struct vp_species *
vp_find_species(const gkyl_vlasov_poisson_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return &app->species[i];
  return 0;
}

int
vp_find_species_idx(const gkyl_vlasov_poisson_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return i;
  return -1;
}

void
gkyl_vlasov_poisson_app_apply_ic(gkyl_vlasov_poisson_app* app, double t0)
{
  app->tcurr = t0;
  for (int i=0; i<app->num_species; ++i)
    gkyl_vlasov_poisson_app_apply_ic_species(app, i, t0);

  // Compute the fields and apply BCs.
  struct gkyl_array *distf[app->num_species];
  for (int i=0; i<app->num_species; ++i) {
    distf[i] = app->species[i].f;
  }
  calc_field_and_apply_bc(app, t0, distf);
}

void
gkyl_vlasov_poisson_app_apply_ic_species(gkyl_vlasov_poisson_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  vp_species_apply_ic(app, &app->species[sidx], t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);
}

void
gkyl_vlasov_poisson_app_calc_mom(gkyl_vlasov_poisson_app* app)
{
  struct timespec wst = gkyl_wall_clock();
  for (int i=0; i<app->num_species; ++i) {
    struct vp_species *s = &app->species[i];

    for (int m=0; m<app->species[i].info.num_diag_moments; ++m) {
      vp_species_moment_calc(&s->moms[m], s->local, app->local, s->f);
      app->stat.nmom += 1;
    }
  }
  app->stat.mom_tm += gkyl_time_diff_now_sec(wst);
}

void
gkyl_vlasov_poisson_app_calc_integrated_mom(gkyl_vlasov_poisson_app* app, double tm)
{
  int vdim = app->vdim;
  double avals_global[2+vdim];

  struct timespec wst = gkyl_wall_clock();

  for (int i=0; i<app->num_species; ++i) {
    struct vp_species *s = &app->species[i];

    struct timespec wst = gkyl_wall_clock();

    vp_species_moment_calc(&s->integ_moms, s->local, app->local, s->f);
    // reduce to compute sum over whole domain, append to diagnostics
    gkyl_array_reduce_range(s->red_integ_diag, s->integ_moms.marr, GKYL_SUM, &app->local);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim,
      s->red_integ_diag, s->red_integ_diag_global);
    if (app->use_gpu) {
      gkyl_cu_memcpy(avals_global, s->red_integ_diag_global, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      memcpy(avals_global, s->red_integ_diag_global, sizeof(double[2+vdim]));
    }
    gkyl_dynvec_append(s->integ_diag, tm, avals_global);

    app->stat.mom_tm += gkyl_time_diff_now_sec(wst);
    app->stat.nmom += 1;
  }

  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_vlasov_poisson_app_calc_field_energy(gkyl_vlasov_poisson_app* app, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  if (app->has_field)
    vp_field_calc_energy(app, tm, app->field);
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_vlasov_poisson_app_write(gkyl_vlasov_poisson_app* app, double tm, int frame)
{
  app->stat.nio += 1;
  struct timespec wtm = gkyl_wall_clock();
  
  if (app->has_field)
    gkyl_vlasov_poisson_app_write_field(app, tm, frame);

  for (int i=0; i<app->num_species; ++i)
    gkyl_vlasov_poisson_app_write_species(app, i, tm, frame);

  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
}

void
gkyl_vlasov_poisson_app_write_field(gkyl_vlasov_poisson_app* app, double tm, int frame)
{
  if (app->has_field) {
    // Copy data from device to host before writing it out.
    if (app->use_gpu)
      gkyl_array_copy(app->field->phi_host, app->field->phi);
  
    struct gkyl_array_meta *mt = vlasov_poisson_array_meta_new( (struct vlasov_poisson_output_meta) {
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
  
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, app->field->phi_host, fileNm);
  
    vlasov_poisson_array_meta_release(mt);
  }
}

void
gkyl_vlasov_poisson_app_write_species(gkyl_vlasov_poisson_app* app, int sidx, double tm, int frame)
{
  struct gkyl_array_meta *mt = vlasov_poisson_array_meta_new( (struct vlasov_poisson_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }
  );

  struct vp_species *vps = &app->species[sidx];

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, vps->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, vps->info.name, frame);

  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(vps->f_host, vps->f);
  }

  gkyl_comm_array_write(vps->comm, &vps->grid, &vps->local, mt, vps->f_host, fileNm);

  if (vps->emit_lo)
    vp_species_emission_write(app, vps, &vps->bc_emission_lo, mt, frame);
  if (app->species[sidx].emit_up)
    vp_species_emission_write(app, vps, &vps->bc_emission_up, mt, frame);

  vlasov_poisson_array_meta_release(mt);
}

void
gkyl_vlasov_poisson_app_write_mom(gkyl_vlasov_poisson_app* app, double tm, int frame)
{
  struct gkyl_array_meta *mt = vlasov_poisson_array_meta_new( (struct vlasov_poisson_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  for (int i=0; i<app->num_species; ++i) {
    for (int m=0; m<app->species[i].info.num_diag_moments; ++m) {

      const char *fmt = "%s-%s_%s_%d.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, app->species[i].info.name,
        app->species[i].info.diag_moments[m], frame);
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[i].info.name,
        app->species[i].info.diag_moments[m], frame);

      if (app->use_gpu)
        gkyl_array_copy(app->species[i].moms[m].marr_host, app->species[i].moms[m].marr);

      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
        app->species[i].moms[m].marr_host, fileNm);
    }
  }

  vlasov_poisson_array_meta_release(mt);
}

void
gkyl_vlasov_poisson_app_write_integrated_mom(gkyl_vlasov_poisson_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out integrated diagnostic moments
      const char *fmt = "%s-%s_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, app->species[i].info.name, "integrated_moms");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[i].info.name, "integrated_moms");

      if (app->species[i].is_first_integ_write_call) {
        gkyl_dynvec_write(app->species[i].integ_diag, fileNm);
        app->species[i].is_first_integ_write_call = false;
      }
      else {
        gkyl_dynvec_awrite(app->species[i].integ_diag, fileNm);
      }
    }
    gkyl_dynvec_clear(app->species[i].integ_diag);
  }
}

void
gkyl_vlasov_poisson_app_write_field_energy(gkyl_vlasov_poisson_app* app)
{
  if (app->has_field) {
    // write out diagnostic moments
    const char *fmt = "%s-field-energy.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name);

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);

    if (rank == 0) {
      if (app->field->is_first_energy_write_call) {
        // write to a new file (this ensure previous output is removed)
        gkyl_dynvec_write(app->field->integ_energy, fileNm);
        app->field->is_first_energy_write_call = false;
      }
      else {
        // append to existing file
        gkyl_dynvec_awrite(app->field->integ_energy, fileNm);
      }
    }
    gkyl_dynvec_clear(app->field->integ_energy);
  }
}

void
gkyl_vlasov_poisson_app_write_lte_corr_status(gkyl_vlasov_poisson_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    struct vp_species *s = &app->species[i];

    if (s->collision_id == GKYL_BGK_COLLISIONS) {
       // write out diagnostic moments
      const char *fmt = "%s-%s-%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, app->species[i].info.name,
        "corr-lte-stat");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[i].info.name,
        "corr-lte-stat");

      int rank;
      gkyl_comm_get_rank(app->comm, &rank);

      if (rank == 0) {
        if (s->bgk.is_first_corr_status_write_call) {
          // write to a new file (this ensure previous output is removed)
          gkyl_dynvec_write(s->bgk.corr_stat, fileNm);
          s->bgk.is_first_corr_status_write_call = false;
        }
        else {
          // append to existing file
          gkyl_dynvec_awrite(s->bgk.corr_stat, fileNm);
        }
      }
      gkyl_dynvec_clear(s->bgk.corr_stat);
    }
  } 
}

// Take a forward Euler step with the suggested time-step dt. This may
// not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by
// stability. The actual time-step and dt_suggested are returned in
// the status object.
static void
forward_euler(gkyl_vlasov_poisson_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], struct gkyl_array *fout[],
  struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;

  double dtmin = DBL_MAX;

  // Compute necessary moments and boundary corrections for collisions.
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      vp_species_lbo_moms(app, &app->species[i], &app->species[i].lbo, fin[i]);
    }
    else if (app->species[i].collision_id == GKYL_BGK_COLLISIONS) {
      vp_species_bgk_moms(app, &app->species[i], 
        &app->species[i].bgk, fin[i]);
    }
  }

  // Compute necessary moments for cross-species collisions.
  // Needs to be done after self-collisions moments, so separate loop over species.
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      if (app->species[i].lbo.num_cross_collisions) {
        vp_species_lbo_cross_moms(app, &app->species[i], &app->species[i].lbo, fin[i]);
      }
    }
  }

  // Compute RHS of Vlasov-Poisson equation.
  for (int i=0; i<app->num_species; ++i) {
    struct vp_species *vps = &app->species[i];
    double dt1 = vp_species_rhs(app, vps, fin[i], fout[i]);
    dtmin = fmin(dtmin, dt1);
  }

  // Compute source term.
  // Done here as the RHS update for all species should be complete before
  // bflux calculation of the source species.
  for (int i=0; i<app->num_species; ++i)
    if (app->species[i].source_id)
      vp_species_source_rhs(app, &app->species[i], &app->species[i].src, fin, fout);

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // compute minimum time-step across all processors
  double dtmin_local = dtmin, dtmin_global;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dtmin_local, &dtmin_global);
  dtmin = dtmin_global;
  
  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // complete update of distribution function
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout[i], dta), 1.0, fin[i]);
  }
}

// Take time-step using the RK3 method. Also sets the status object
// which has the actual and suggested dts used. These can be different
// from the actual time-step.
static struct gkyl_update_status
rk3(gkyl_vlasov_poisson_app* app, double dt0)
{
  const struct gkyl_array *fin[app->num_species];
  struct gkyl_array *fout[app->num_species];
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

        forward_euler(app, tcurr, dt, fin, fout, &st);
	// Compute the fields and apply BCs.
        calc_field_and_apply_bc(app, tcurr, fout);

        dt = st.dt_actual;
        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
        }

        forward_euler(app, tcurr+dt, dt, fin, fout, &st);

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

        } else {
          for (int i=0; i<app->num_species; ++i)
            array_combine(app->species[i].f1,
              3.0/4.0, app->species[i].f, 1.0/4.0, app->species[i].fnew, &app->species[i].local_ext);

	  // Compute the fields and apply BCs.
          for (int i=0; i<app->num_species; ++i) {
            fout[i] = app->species[i].f1;
          }
          calc_field_and_apply_bc(app, tcurr, fout);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
        }

        forward_euler(app, tcurr+dt/2, dt, fin, fout, &st);

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

	  // Compute the fields and apply BCs
          for (int i=0; i<app->num_species; ++i) {
            fout[i] = app->species[i].f;
          }
          calc_field_and_apply_bc(app, tcurr, fout);

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
gkyl_vlasov_poisson_update(gkyl_vlasov_poisson_app* app, double dt)
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

struct gkyl_vlasov_poisson_stat
gkyl_vlasov_poisson_app_stat(gkyl_vlasov_poisson_app* app)
{
  vp_species_tm(app);
  vp_species_coll_tm(app);
  return app->stat;
}

void
gkyl_vlasov_poisson_app_species_ktm_rhs(gkyl_vlasov_poisson_app* app, int update_vol_term)
{
  for (int i=0; i<app->num_species; ++i) {

    struct vp_species *species = &app->species[i];

    const struct gkyl_array *fin = species->f;
    struct gkyl_array *rhs = species->f1;

    gkyl_array_clear(rhs, 0.0);
    gkyl_dg_updater_vlasov_poisson_advance(species->slvr, &species->local, 
      fin, species->cflrate, rhs); 
  }
}

static void
range_stat_write(gkyl_vlasov_poisson_app* app, const char *nm, const struct gkyl_range *r, FILE *fp)
{
  gkyl_vlasov_poisson_app_cout(app, fp, " %s_cells : [ ", nm);
  for (int i=0; i<r->ndim; ++i)
    gkyl_vlasov_poisson_app_cout(app, fp, " %d, ", gkyl_range_shape(r, i));
  gkyl_vlasov_poisson_app_cout(app, fp, " ],\n");
}

// ensure stats across processors are made consistent
static void
comm_reduce_app_stat(const gkyl_vlasov_poisson_app* app,
  const struct gkyl_vlasov_poisson_stat *local, struct gkyl_vlasov_poisson_stat *global)
{
  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);
  if (comm_sz == 1) {
    memcpy(global, local, sizeof(struct gkyl_vlasov_poisson_stat));
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
  gkyl_comm_allreduce(app->comm, GKYL_INT_64, GKYL_MAX, L_END, l_red, l_red_global);

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
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, D_END, d_red, d_red_global);
  
  global->total_tm = d_red_global[TOTAL_TM];
  global->init_species_tm = d_red_global[INIT_SPECIES_TM];
  global->species_coll_mom_tm = d_red_global[SPECIES_COLL_MOM_TM];
  global->species_coll_tm = d_red_global[SPECIES_COL_TM];
  global->species_rhs_tm = d_red_global[SPECIES_RHS_TM];
  global->field_rhs_tm = d_red_global[FIELD_RHS_TM];
  global->species_omega_cfl_tm = d_red_global[SPECIES_OMEGA_CFL_TM];
  global->mom_tm = d_red_global[MOM_TM];
  global->diag_tm = d_red_global[DIAG_TM];
  global->io_tm = d_red_global[IO_TM];
  global->species_bc_tm = d_red_global[SPECIES_BC_TM];

  // misc data needing reduction

  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_2_dt_diff,
    global->stage_2_dt_diff);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_3_dt_diff,
    global->stage_3_dt_diff);

  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, GKYL_MAX_SPECIES, local->species_lbo_coll_drag_tm,
    global->species_lbo_coll_drag_tm);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, GKYL_MAX_SPECIES, local->species_lbo_coll_diff_tm,
    global->species_lbo_coll_diff_tm);
}

void
gkyl_vlasov_poisson_app_stat_write(gkyl_vlasov_poisson_app* app)
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

  vp_species_coll_tm(app);
  vp_species_tm(app);

  struct gkyl_vlasov_poisson_stat stat = { };
  comm_reduce_app_stat(app, &app->stat, &stat);
  
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  if (rank == 0) fp = fopen(fileNm, "a");

  gkyl_vlasov_poisson_app_cout(app, fp, "{\n");

  if (strftime(buff, sizeof buff, "%c", &curr_tm))
    gkyl_vlasov_poisson_app_cout(app, fp, " date : %s,\n", buff);

  gkyl_vlasov_poisson_app_cout(app, fp, " use_gpu : %d,\n", stat.use_gpu);
  gkyl_vlasov_poisson_app_cout(app, fp, " num_ranks : %d,\n", num_ranks); 
  
  for (int s=0; s<app->num_species; ++s)
    range_stat_write(app, app->species[s].info.name, &app->species[s].global, fp);
  
  gkyl_vlasov_poisson_app_cout(app, fp, " nup : %ld,\n", stat.nup);
  gkyl_vlasov_poisson_app_cout(app, fp, " nfeuler : %ld,\n", stat.nfeuler);
  gkyl_vlasov_poisson_app_cout(app, fp, " nstage_2_fail : %ld,\n", stat.nstage_2_fail);
  gkyl_vlasov_poisson_app_cout(app, fp, " nstage_3_fail : %ld,\n", stat.nstage_3_fail);

  gkyl_vlasov_poisson_app_cout(app, fp, " stage_2_dt_diff : [ %lg, %lg ],\n",
    stat.stage_2_dt_diff[0], stat.stage_2_dt_diff[1]);
  gkyl_vlasov_poisson_app_cout(app, fp, " stage_3_dt_diff : [ %lg, %lg ],\n",
    stat.stage_3_dt_diff[0], stat.stage_3_dt_diff[1]);

  gkyl_vlasov_poisson_app_cout(app, fp, " total_tm : %lg,\n", stat.total_tm);
  gkyl_vlasov_poisson_app_cout(app, fp, " init_species_tm : %lg,\n", stat.init_species_tm);
  
  gkyl_vlasov_poisson_app_cout(app, fp, " species_rhs_tm : %lg,\n", stat.species_rhs_tm);
  if (app->has_field)
    gkyl_vlasov_poisson_app_cout(app, fp, " field_rhs_tm : %lg,\n", stat.field_rhs_tm);

  for (int s=0; s<app->num_species; ++s) {
    gkyl_vlasov_poisson_app_cout(app, fp, " species_coll_drag_tm[%d] : %lg,\n", s,
      stat.species_lbo_coll_drag_tm[s]);
    gkyl_vlasov_poisson_app_cout(app, fp, " species_coll_diff_tm[%d] : %lg,\n", s,
      stat.species_lbo_coll_diff_tm[s]);
  }

  gkyl_vlasov_poisson_app_cout(app, fp, " species_coll_mom_tm : %lg,\n", stat.species_coll_mom_tm);
  gkyl_vlasov_poisson_app_cout(app, fp, " species_coll_tm : %lg,\n", stat.species_coll_tm);

  gkyl_vlasov_poisson_app_cout(app, fp, " species_bc_tm : %lg,\n", stat.species_bc_tm);

  gkyl_vlasov_poisson_app_cout(app, fp, " nmom : %ld,\n", stat.nmom);
  gkyl_vlasov_poisson_app_cout(app, fp, " mom_tm : %lg\n", stat.mom_tm);

  gkyl_vlasov_poisson_app_cout(app, fp, " ndiag : %ld,\n", stat.ndiag);
  gkyl_vlasov_poisson_app_cout(app, fp, " diag_tm : %lg\n", stat.diag_tm);
  
  gkyl_vlasov_poisson_app_cout(app, fp, " nspecies_omega_cfl : %ld,\n", stat.nspecies_omega_cfl);
  gkyl_vlasov_poisson_app_cout(app, fp, " species_omega_cfl_tm : %lg\n", stat.species_omega_cfl_tm);

  gkyl_vlasov_poisson_app_cout(app, fp, " nio : %ld,\n", stat.nio);
  gkyl_vlasov_poisson_app_cout(app, fp, " io_tm : %lg\n", stat.io_tm);
  
  gkyl_vlasov_poisson_app_cout(app, fp, "}\n");

  if (rank == 0)
    fclose(fp);  

}

static struct gkyl_app_restart_status
header_from_file(gkyl_vlasov_poisson_app *app, const char *fname)
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

    struct vlasov_poisson_output_meta meta =
      vlasov_poisson_meta_from_mpack( &(struct gkyl_array_meta) {
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

struct gkyl_app_restart_status
gkyl_vlasov_poisson_app_from_file_field(gkyl_vlasov_poisson_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = header_from_file(app, fname);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, app->field->phi_host, fname);
    if (app->use_gpu)
      gkyl_array_copy(app->field->phi, app->field->phi_host);
  }

  return rstat;
}

struct gkyl_app_restart_status
gkyl_vlasov_poisson_app_from_file_species(gkyl_vlasov_poisson_app *app, int sidx,
  const char *fname)
{
  struct gkyl_app_restart_status rstat = header_from_file(app, fname);

  struct vp_species *vps = &app->species[sidx];

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(vps->comm, &vps->grid, &vps->local, vps->f_host, fname);
    if (app->use_gpu)
      gkyl_array_copy(vps->f, vps->f_host);
  }

  return rstat;
}

struct gkyl_app_restart_status
gkyl_vlasov_poisson_app_from_frame_field(gkyl_vlasov_poisson_app *app, int frame)
{
  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "field", frame);
  struct gkyl_app_restart_status rstat = gkyl_vlasov_poisson_app_from_file_field(app, fileNm.str);
  app->field->is_first_energy_write_call = false; // append to existing diagnostic
  cstr_drop(&fileNm);

  return rstat;
}

struct gkyl_app_restart_status
gkyl_vlasov_poisson_app_from_frame_species(gkyl_vlasov_poisson_app *app, int sidx, int frame)
{
  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, app->species[sidx].info.name, frame);
  struct gkyl_app_restart_status rstat = gkyl_vlasov_poisson_app_from_file_species(app, sidx, fileNm.str);
  app->species[sidx].is_first_integ_write_call = false; // append to existing diagnostic
  cstr_drop(&fileNm);

  vp_species_bflux_rhs(app, &app->species[sidx], &app->species[sidx].bflux, app->species[sidx].f,
    app->species[sidx].f);

  return rstat;
}

struct gkyl_app_restart_status
gkyl_vlasov_poisson_app_read_from_frame(gkyl_vlasov_poisson_app *app, int frame)
{
  struct gkyl_app_restart_status rstat;
  for (int i=0; i<app->num_species; i++) {
    rstat = gkyl_vlasov_poisson_app_from_frame_species(app, i, frame);
  }

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    // Compute the fields and apply BCs.
    struct gkyl_array *distf[app->num_species];
    for (int i=0; i<app->num_species; ++i) {
      distf[i] = app->species[i].f;
    }
    calc_field_and_apply_bc(app, rstat.stime, distf);

    if (app->field->has_ext_em) {
      // Compute the external fields.
      vp_field_calc_ext_em(app, app->field, rstat.stime);
    }
  }

  return rstat;
}

// private function to handle variable argument list for printing
static void
v_vlasov_poisson_app_cout(const gkyl_vlasov_poisson_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank, r = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp)
    vfprintf(fp, fmt, argp);
}

void
gkyl_vlasov_poisson_app_cout(const gkyl_vlasov_poisson_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_vlasov_poisson_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_vlasov_poisson_app_release(gkyl_vlasov_poisson_app* app)
{
  if (app->has_field)
    vp_field_release(app, app->field);

  for (int i=0; i<app->num_species; ++i)
    vp_species_release(app, &app->species[i]);
  if (app->num_species > 0)
    gkyl_free(app->species);

  gkyl_comm_release(app->comm);

  if (app->use_gpu) {
    gkyl_cu_free(app->basis_on_dev.basis);
    gkyl_cu_free(app->basis_on_dev.confBasis);
  }

  gkyl_free(app);
}
