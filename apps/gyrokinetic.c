#include <stdarg.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>
#include <gkyl_null_comm.h>

#include <gkyl_gyrokinetic_priv.h>

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
  app->use_gpu = gk->use_gpu;
#else
  app->use_gpu = false; // can't use GPUs if we don't have them!
#endif

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

  if (gk->has_low_inp) {
    // create local and local_ext from user-supplied local range
    gkyl_create_ranges(&gk->low_inp.local_range, ghost, &app->local_ext, &app->local);
    
    if (gk->low_inp.comm)
      app->comm = gkyl_comm_acquire(gk->low_inp.comm);
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
    // global and local ranges are same, and so just copy
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
  // local skin and ghost ranges for configuration space fields
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&app->lower_skin[dir], &app->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &app->local_ext, ghost); 
    gkyl_skin_ghost_ranges(&app->upper_skin[dir], &app->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &app->local_ext, ghost);
  }


  // Configuration space geometry initialization
  struct gkyl_rect_grid geo_grid;
  struct gkyl_range geo_local;
  struct gkyl_range geo_local_ext;
  struct gkyl_basis geo_basis;
  bool geo_3d_use_gpu = app->use_gpu;

  if(app->cdim < 3){
    geo_grid = agument_grid(app->grid, gk->geometry);
    gkyl_create_grid_ranges(&geo_grid, ghost, &geo_local_ext, &geo_local);
    geo_3d_use_gpu = false;
    switch (gk->basis_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
        gkyl_cart_modal_serendip(&geo_basis, 3, poly_order);
        break;
      default:
        assert(false);
        break;
    }
  }
  else{
    geo_grid = app->grid;
    geo_local = app->local;
    geo_local_ext = app->local_ext;
    geo_basis = app->confBasis;
  }

  struct gk_geometry* gk_geom_3d;
  switch (gk->geometry.geometry_id) {
    case GKYL_GEOMETRY_FROMFILE:
      gk_geom_3d = gkyl_gk_geometry_fromfile_new(&app->grid, &app->local, &app->local_ext, &app->confBasis, app->use_gpu);
      break;
    case GKYL_TOKAMAK:
      gk_geom_3d = gkyl_gk_geometry_tok_new(&geo_grid, &geo_local, &geo_local_ext, &geo_basis, 
          gk->geometry.tok_efit_info, gk->geometry.tok_grid_info, geo_3d_use_gpu);
      break;
    case GKYL_MIRROR:
      gk_geom_3d = gkyl_gk_geometry_mirror_new(&geo_grid, &geo_local, &geo_local_ext, &geo_basis, 
          gk->geometry.mirror_efit_info, gk->geometry.mirror_grid_info, geo_3d_use_gpu);
      break;
    case GKYL_MAPC2P:
      gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geo_grid, &geo_local, &geo_local_ext, &geo_basis, 
          gk->geometry.mapc2p, gk->geometry.c2p_ctx, gk->geometry.bmag_func,  gk->geometry.bmag_ctx, geo_3d_use_gpu);
      break;
  }

  // deflate geometry if necessary
  if(app->cdim < 3 && gk->geometry.geometry_id != GKYL_GEOMETRY_FROMFILE)
    app->gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &app->grid, &app->local, &app->local_ext, 
        &app->confBasis, app->use_gpu);
  else
    app->gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);

  gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry

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
    // initialize cross-species collisions (e.g, LBO or BGK)
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      if (app->species[i].lbo.num_cross_collisions) {
        gk_species_lbo_cross_init(app, &app->species[i], &app->species[i].lbo);
      }
    }
    else if (app->species[i].collision_id == GKYL_BGK_COLLISIONS) {
      if (app->species[i].bgk.num_cross_collisions) {
        gk_species_bgk_cross_init(app, &app->species[i], &app->species[i].bgk);
      }
    }
    // initialize cross-species reactions with plasma species (e.g., ionization, recombination, or charge exchange)
    if (app->species[i].has_reactions) {
      gk_species_react_cross_init(app, &app->species[i], &app->species[i].react);
    }
    // initialize cross-species reactions with neutral species (e.g., ionization, recombination, or charge exchange)
    if (app->species[i].has_neutral_reactions) {
      gk_species_react_cross_init(app, &app->species[i], &app->species[i].react_neut);
    }
    // initial radiation (e.g., bremmstrahlung model from cross-collisions of electrons with ions)
    if (app->species[i].radiation_id == GKYL_GK_RADIATION) {
      gk_species_radiation_init(app, &app->species[i], &app->species[i].rad);
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
    if (app->species[i].source_id) {
      gk_species_source_init(app, &app->species[i], &app->species[i].src);
    }
  }
  for (int i=0; i<neuts; ++i) {
    if (app->neut_species[i].source_id) {
      gk_neut_species_source_init(app, &app->neut_species[i], &app->neut_species[i].src);
    }
  }

  // initialize stat object
  app->stat = (struct gkyl_gyrokinetic_stat) {
    .use_gpu = app->use_gpu,
    .stage_2_dt_diff = { DBL_MAX, 0.0 },
    .stage_3_dt_diff = { DBL_MAX, 0.0 },
  };

  return app;
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
}

void
gkyl_gyrokinetic_app_apply_ic_species(gkyl_gyrokinetic_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  gk_species_apply_ic(app, &app->species[sidx], t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);

  gk_species_apply_bc(app, &app->species[sidx], app->species[sidx].f);
}

void
gkyl_gyrokinetic_app_apply_ic_neut_species(gkyl_gyrokinetic_app* app, int sidx, double t0)
{
  assert(sidx < app->num_neut_species);

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  gk_neut_species_apply_ic(app, &app->neut_species[sidx], t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);

  gk_neut_species_apply_bc(app, &app->neut_species[sidx], app->neut_species[sidx].f);
}

void
gkyl_gyrokinetic_app_calc_mom(gkyl_gyrokinetic_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];

    for (int m=0; m<app->species[i].info.num_diag_moments; ++m) {
      struct timespec wst = gkyl_wall_clock();
      gk_species_moment_calc(&s->moms[m], s->local, app->local, s->f);
      app->stat.mom_tm += gkyl_time_diff_now_sec(wst);
      app->stat.nmom += 1;
    }
  }
}

void
gkyl_gyrokinetic_app_calc_integrated_mom(gkyl_gyrokinetic_app* app, double tm)
{
  int vdim = app->vdim;
  double avals[2+vdim], avals_global[2+vdim];

  struct timespec wst = gkyl_wall_clock();

  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];

    struct timespec wst = gkyl_wall_clock();

    gk_species_moment_calc(&s->integ_moms, s->local, app->local, s->f);

    // Rescale integrated moment by inverse of Jacobian before integration over configuration space
    gkyl_dg_div_op_range(app->species[i].integ_moms.mem_geo, app->confBasis, 
      0, app->species[i].integ_moms.marr, 0, app->species[i].integ_moms.marr, 0, app->gk_geom->jacobgeo, &app->local); 

    // reduce to compute sum over whole domain, append to diagnostics
    if (app->use_gpu) {
      gkyl_array_reduce_range(s->red_integ_diag, s->integ_moms.marr, GKYL_SUM, &(app->local));
      gkyl_cu_memcpy(avals, s->red_integ_diag, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      gkyl_array_reduce_range(avals, s->integ_moms.marr_host, GKYL_SUM, &(app->local));
    }

    gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, avals, avals_global);
    gkyl_dynvec_append(s->integ_diag, tm, avals_global);

    app->stat.mom_tm += gkyl_time_diff_now_sec(wst);
    app->stat.nmom += 1;
  }

  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_gyrokinetic_app_calc_field_energy(gkyl_gyrokinetic_app* app, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  gk_field_calc_energy(app, tm, app->field);
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_gyrokinetic_app_write(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  app->stat.nio += 1;
  struct timespec wtm = gkyl_wall_clock();
  
  gkyl_gyrokinetic_app_write_field(app, tm, frame);

  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_app_write_species(app, i, tm, frame);
    if (app->species[i].source_id)
      if (app->species[i].src.write_source)
        gkyl_gyrokinetic_app_write_source_species(app, i, tm, frame);
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS)
      gkyl_gyrokinetic_app_write_coll_mom(app, i, tm, frame);
    if (app->species[i].radiation_id == GKYL_GK_RADIATION)
      gkyl_gyrokinetic_app_write_rad_drag(app, i, tm, frame);
  }

  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
}

void
gkyl_gyrokinetic_app_write_field(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  const char *fmt = "%s-field_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);

  // Compute electrostatic potential at desired output time
  const struct gkyl_array *fin[app->num_species];
  struct gkyl_array *fout[app->num_species];
  for (int i=0; i<app->num_species; ++i) {
    fin[i] = app->species[i].f;
    fout[i] = app->species[i].f1;
  }
  gk_field_accumulate_rho_c(app, app->field, fin);
  if (app->field->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
    gk_field_calc_ambi_pot_sheath_vals(app, app->field, fin, fout);
  }
  gk_field_rhs(app, app->field);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->field->phi_host, app->field->phi_smooth);
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, app->field->phi_host, fileNm);
  }
  else {
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, app->field->phi_smooth, fileNm);
  }
}

void
gkyl_gyrokinetic_app_write_species(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, app->species[sidx].info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[sidx].info.name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->species[sidx].f_host, app->species[sidx].f);
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, &app->species[sidx].local,
      app->species[sidx].f_host, fileNm);
  }
  else {
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, &app->species[sidx].local,
      app->species[sidx].f, fileNm);
  }
}

void
gkyl_gyrokinetic_app_write_source_species(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *s = &app->species[sidx];

  // Write out the source distribution function
  const char *fmt = "%s-%s_source_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, s->info.name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(s->src.source_host, s->src.source);
    gkyl_comm_array_write(s->comm, &s->grid, &s->local, s->src.source_host, fileNm);
  }
  else {
    gkyl_comm_array_write(s->comm, &s->grid, &s->local, s->src.source, fileNm);
  }

  // Write out density of the source
  const char *fmt_M0 = "%s-%s_source_M0_%d.gkyl";
  int sz_M0 = gkyl_calc_strlen(fmt_M0, app->name, s->info.name, frame);
  char fileNm_M0[sz_M0+1]; // ensures no buffer overflow
  snprintf(fileNm_M0, sizeof fileNm_M0, fmt_M0, app->name, s->info.name, frame); 

  gk_species_moment_calc(&s->m0, s->local, app->local, s->src.source); 
  // Rescale moment by inverse of Jacobian
  gkyl_dg_div_op_range(s->m0.mem_geo, app->confBasis, 
    0, s->m0.marr, 0, s->m0.marr, 0, app->gk_geom->jacobgeo, &app->local); 

  if (app->use_gpu)
    gkyl_array_copy(s->m0.marr_host, s->m0.marr);

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, s->m0.marr_host, fileNm_M0);   
}

void
gkyl_gyrokinetic_app_write_coll_mom(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *s = &app->species[sidx];

  // Construct the file handles for collision frequency and primitive moments
  const char *fmt = "%s-%s_nu_sum_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, s->info.name, frame);

  const char *fmt_prim = "%s-%s_prim_moms_%d.gkyl";
  int sz_prim = gkyl_calc_strlen(fmt_prim, app->name, s->info.name, frame);
  char fileNm_prim[sz_prim+1]; // ensures no buffer overflow
  snprintf(fileNm_prim, sizeof fileNm_prim, fmt_prim, app->name, s->info.name, frame);

  const char *fmt_nu_prim = "%s-%s_nu_prim_moms_%d.gkyl";
  int sz_nu_prim = gkyl_calc_strlen(fmt_nu_prim, app->name, s->info.name, frame);
  char fileNm_nu_prim[sz_nu_prim+1]; // ensures no buffer overflow
  snprintf(fileNm_nu_prim, sizeof fileNm_nu_prim, fmt_nu_prim, app->name, s->info.name, frame);

  // Compute primitive moments
  const struct gkyl_array *fin[app->num_species];
  gk_species_lbo_moms(app, s, &s->lbo, s->f);
  if (s->lbo.num_cross_collisions)
    gk_species_lbo_moms(app, s, &s->lbo, s->f);

  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(s->lbo.nu_sum_host, s->lbo.nu_sum);
    gkyl_array_copy(s->lbo.prim_moms_host, s->lbo.prim_moms);
    gkyl_array_copy(s->lbo.nu_prim_moms_host, s->lbo.nu_prim_moms);
  }

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, s->lbo.nu_sum_host, fileNm);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, s->lbo.prim_moms_host, fileNm_prim);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, s->lbo.nu_prim_moms_host, fileNm_nu_prim);
}

void
gkyl_gyrokinetic_app_write_rad_drag(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *s = &app->species[sidx];

  // Construct the file handles for vparallel and mu drag
  const char *fmt_vnu = "%s-%s_vnu_%d.gkyl";
  int sz_vnu = gkyl_calc_strlen(fmt_vnu, app->name, s->info.name, frame);
  char fileNm_vnu[sz_vnu+1]; // ensures no buffer overflow
  snprintf(fileNm_vnu, sizeof fileNm_vnu, fmt_vnu, app->name, s->info.name, frame);

  const char *fmt_vsqnu = "%s-%s_vsqnu_%d.gkyl";
  int sz_vsqnu = gkyl_calc_strlen(fmt_vsqnu, app->name, s->info.name, frame);
  char fileNm_vsqnu[sz_vsqnu+1]; // ensures no buffer overflow
  snprintf(fileNm_vsqnu, sizeof fileNm_vsqnu, fmt_vsqnu, app->name, s->info.name, frame);

  // Compute radiation drag coefficients
  const struct gkyl_array *fin[app->num_species];
  for (int i=0; i<app->num_species; ++i) 
    fin[i] = app->species[i].f;
  gk_species_radiation_moms(app, s, &s->rad, fin);

  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(s->rad.nvnu_sum_host, s->rad.nvnu_sum);
    gkyl_array_copy(s->rad.nvsqnu_sum_host, s->rad.nvsqnu_sum);
  }

  gkyl_comm_array_write(s->comm, &s->grid, &s->local, s->rad.nvnu_sum_host, fileNm_vnu);
  gkyl_comm_array_write(s->comm, &s->grid, &s->local, s->rad.nvsqnu_sum_host, fileNm_vsqnu);
}

void
gkyl_gyrokinetic_app_write_mom(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    for (int m=0; m<app->species[i].info.num_diag_moments; ++m) {

      const char *fmt = "%s-%s_%s_%d.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, app->species[i].info.name,
        app->species[i].info.diag_moments[m], frame);
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[i].info.name,
        app->species[i].info.diag_moments[m], frame);

      // Rescale moment by inverse of Jacobian
      gkyl_dg_div_op_range(app->species[i].moms[m].mem_geo, app->confBasis, 
        0, app->species[i].moms[m].marr, 0, app->species[i].moms[m].marr, 0, app->gk_geom->jacobgeo, &app->local);      

      if (app->use_gpu)
        gkyl_array_copy(app->species[i].moms[m].marr_host, app->species[i].moms[m].marr);

      gkyl_comm_array_write(app->comm, &app->grid, &app->local, app->species[i].moms[m].marr_host, fileNm);
    }
  }
}

void
gkyl_gyrokinetic_app_write_integrated_mom(gkyl_gyrokinetic_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out integrated diagnostic moments
      const char *fmt = "%s-%s-%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, app->species[i].info.name,
        "imom");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[i].info.name,
        "imom");

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
gkyl_gyrokinetic_app_write_field_energy(gkyl_gyrokinetic_app* app)
{
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

  if (app->update_field) {
    // Compute biased wall potential if present and time-dependent.
    // Note: biased wall potential use proj_on_basis 
    // so does copy to GPU every call if app->use_gpu = true.
    if (app->field->phi_wall_lo_evolve || app->field->phi_wall_up_evolve)
      gk_field_calc_phi_wall(app, app->field, tcurr);

    // compute electrostatic potential from gyrokinetic Poisson's equation
    gk_field_accumulate_rho_c(app, app->field, fin);
    gk_field_rhs(app, app->field);
  }

  // compute necessary moments and boundary corrections for collisions
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      gk_species_lbo_moms(app, &app->species[i], 
        &app->species[i].lbo, fin[i]);
    }
    else if (app->species[i].collision_id == GKYL_BGK_COLLISIONS) {
      gk_species_bgk_moms(app, &app->species[i], 
        &app->species[i].bgk, fin[i]);
    }
  }

  // compute necessary moments for cross-species collisions
  // needs to be done after self-collisions moments, so separate loop over species
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) { 
      if (app->species[i].lbo.num_cross_collisions) {
        gk_species_lbo_cross_moms(app, &app->species[i], 
          &app->species[i].lbo, fin[i]);        
      }
    }
    else if (app->species[i].collision_id == GKYL_BGK_COLLISIONS) {
      if (app->species[i].bgk.num_cross_collisions) {
        gk_species_bgk_cross_moms(app, &app->species[i], 
          &app->species[i].bgk, fin[i]);        
      }
    }
    // compute necessary reaction rates (e.g., ionization, recombination, or charge exchange)
    if (app->species[i].has_reactions) {
      gk_species_react_cross_moms(app, &app->species[i], 
        &app->species[i].react, fin[i], fin, fin_neut);
    }
    if (app->species[i].has_neutral_reactions) {
      gk_species_react_cross_moms(app, &app->species[i], 
        &app->species[i].react_neut, fin[i], fin, fin_neut);
    }
    // compute necessary drag coefficients for radiation operator
    if (app->species[i].radiation_id == GKYL_GK_RADIATION) {
      gk_species_radiation_moms(app, &app->species[i], 
        &app->species[i].rad, fin);
    }
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    // compute necessary reaction rates (e.g., ionization, recombination, or charge exchange)
    if (app->neut_species[i].has_neutral_reactions) {
      gk_neut_species_react_cross_moms(app, &app->neut_species[i], 
        &app->neut_species[i].react_neut, fin[i], fin, fin_neut);
    }
  }

  // compute RHS of Gyrokinetic equation
  for (int i=0; i<app->num_species; ++i) {
    double dt1 = gk_species_rhs(app, &app->species[i], fin[i], fout[i]);
    dtmin = fmin(dtmin, dt1);
  }

  // compute RHS of Neutrals 
  for (int i=0; i<app->num_neut_species; ++i) {
    double dt1 = gk_neut_species_rhs(app, &app->neut_species[i], fin_neut[i], fout_neut[i]);
    dtmin = fmin(dtmin, dt1);
  }

  // compute plasma source term
  // done here as the RHS update for all species should be complete before
  // in case we are using boundary fluxes as a component of our source function
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].source_id) {
      gk_species_source_rhs(app, &app->species[i], 
        &app->species[i].src, fin[i], fout[i]);
    }
  }

  // compute neutral source term
  // done here as the RHS update for all species should be complete before
  // in case we are using boundary fluxes as a component of our source function
  for (int i=0; i<app->num_neut_species; ++i) {
    if (app->neut_species[i].source_id) {
      gk_neut_species_source_rhs(app, &app->neut_species[i], 
        &app->neut_species[i].src, fin_neut[i], fout_neut[i]);
    }
  }

  if (app->update_field) {
    // Compute ambipolar potential sheath values if using adiabatic electrons
    // done here as the RHS update for all species should be complete before
    // boundary fluxes are computed (ion fluxes needed for sheath values) 
    // and these boundary fluxes are stored temporarily in ghost cells of RHS
    if (app->field->gkfield_id == GKYL_GK_FIELD_ADIABATIC)
      gk_field_calc_ambi_pot_sheath_vals(app, app->field, fin, fout);
  }

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // compute minimum time-step across all processors
  double dtmin_local = dtmin, dtmin_global;
  gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dtmin_local, &dtmin_global);
  dtmin = dtmin_global;
  
  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // complete update of distribution functions and apply boundary conditions
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout[i], dta), 1.0, fin[i]);
    gk_species_apply_bc(app, &app->species[i], fout[i]);
  }

  // complete update of neutral distribution functions and apply boundary conditions
  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout_neut[i], dta), 1.0, fin_neut[i]);
    gk_neut_species_apply_bc(app, &app->neut_species[i], fout_neut[i]);
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
          fout_neut[i] = app->neut_species[i].f1;
        }
        forward_euler(app, tcurr, dt, fin, fout, fin_neut, fout_neut, &st);
        dt = st.dt_actual;
        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
        }
        for (int i=0; i<app->num_neut_species; ++i) {
          fin_neut[i] = app->neut_species[i].f1;
          fout_neut[i] = app->neut_species[i].fnew;
        }
        forward_euler(app, tcurr+dt, dt, fin, fout, fin_neut, fout_neut, &st);
        if (st.dt_actual < dt) {

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
            array_combine(app->neut_species[i].f1,
              3.0/4.0, app->neut_species[i].f, 1.0/4.0, app->neut_species[i].fnew, &app->neut_species[i].local_ext);
          }
          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
        }
        for (int i=0; i<app->num_neut_species; ++i) {
          fin_neut[i] = app->neut_species[i].f1;
          fout_neut[i] = app->neut_species[i].fnew;
        }
        forward_euler(app, tcurr+dt/2, dt, fin, fout, fin_neut, fout_neut, &st);
        if (st.dt_actual < dt) {
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
            array_combine(app->neut_species[i].f1,
              1.0/3.0, app->neut_species[i].f, 2.0/3.0, app->neut_species[i].fnew, &app->neut_species[i].local_ext);
            gkyl_array_copy_range(app->neut_species[i].f, app->neut_species[i].f1, &app->neut_species[i].local_ext);
          }
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
  gkyl_comm_all_reduce(app->comm, GKYL_INT_64, GKYL_MAX, L_END, l_red, l_red_global);

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
  gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_MAX, D_END, d_red, d_red_global);
  
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

  gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_2_dt_diff,
    global->stage_2_dt_diff);
  gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_3_dt_diff,
    global->stage_3_dt_diff);

  gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_MAX, GKYL_MAX_SPECIES, local->species_lbo_coll_drag_tm,
    global->species_lbo_coll_drag_tm);
  gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_MAX, GKYL_MAX_SPECIES, local->species_lbo_coll_diff_tm,
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

// private function to handle variable argument list for printing
static void
v_gk_app_cout(const gkyl_gyrokinetic_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank, r = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0)
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
  gkyl_gk_geometry_release(app->gk_geom);

  gk_field_release(app, app->field);

  for (int i=0; i<app->num_species; ++i)
    gk_species_release(app, &app->species[i]);
  if (app->num_species > 0)
    gkyl_free(app->species);

  gkyl_comm_release(app->comm);

  if (app->use_gpu) {
    gkyl_cu_free(app->basis_on_dev.basis);
    gkyl_cu_free(app->basis_on_dev.confBasis);
  }

  gkyl_free(app);
}
