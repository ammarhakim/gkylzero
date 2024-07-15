#include <stdarg.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_comm_io.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>
#include <gkyl_null_comm.h>
#include <gkyl_pkpm_priv.h>

gkyl_pkpm_app*
gkyl_pkpm_app_new(struct gkyl_pkpm *pkpm)
{
  disable_denorm_float();

  assert(pkpm->num_species <= GKYL_MAX_SPECIES);

  gkyl_pkpm_app *app = gkyl_malloc(sizeof(gkyl_pkpm_app));

  int cdim = app->cdim = pkpm->cdim;
  int vdim = app->vdim = pkpm->vdim;
  int pdim = cdim+vdim;
  int poly_order = app->poly_order = pkpm->poly_order;
  int ns = app->num_species = pkpm->num_species;

  double cfl_frac = pkpm->cfl_frac == 0 ? 1.0 : pkpm->cfl_frac;
  app->cfl = cfl_frac;

#ifdef GKYL_HAVE_CUDA
  app->use_gpu = pkpm->use_gpu;
#else
  app->use_gpu = false; // can't use GPUs if we don't have them!
#endif

  app->num_periodic_dir = pkpm->num_periodic_dir;
  for (int d=0; d<cdim; ++d)
    app->periodic_dirs[d] = pkpm->periodic_dirs[d];

  strcpy(app->name, pkpm->name);
  app->tcurr = 0.0; // reset on init

  if (app->use_gpu) {
    // allocate device basis if we are using GPUs
    app->basis_on_dev.basis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    app->basis_on_dev.confBasis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  }
  else {
    app->basis_on_dev.basis = &app->basis;
    app->basis_on_dev.confBasis = &app->confBasis;
  }

  // basis functions
  switch (pkpm->basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);
      if (poly_order > 1) {
        gkyl_cart_modal_serendip(&app->basis, pdim, poly_order);
        if (vdim > 0)
          gkyl_cart_modal_serendip(&app->velBasis, vdim, poly_order);
      } else if (poly_order == 1) {
        /* Force hybrid basis (p=2 in velocity space). */
        gkyl_cart_modal_hybrid(&app->basis, cdim, vdim);
        if (vdim > 0)
          gkyl_cart_modal_serendip(&app->velBasis, vdim, 2);
      }

      if (app->use_gpu) {
        gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.confBasis, cdim, poly_order);
        if (poly_order > 1) {
          gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.basis, pdim, poly_order);
        } else if (poly_order == 1) {
          /* Force hybrid basis (p=2 in velocity space). */
          gkyl_cart_modal_hybrid_cu_dev(app->basis_on_dev.basis, cdim, vdim);
        }
      }
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      gkyl_cart_modal_tensor(&app->basis, pdim, poly_order);
      gkyl_cart_modal_tensor(&app->confBasis, cdim, poly_order);
      if (vdim > 0)
        gkyl_cart_modal_tensor(&app->velBasis, vdim, poly_order);
      break;

    default:
      assert(false);
      break;
  }

  gkyl_rect_grid_init(&app->grid, cdim, pkpm->lower, pkpm->upper, pkpm->cells);

  int ghost[] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&app->grid, ghost, &app->global_ext, &app->global);
  if (pkpm->has_low_inp) {
    // create local and local_ext from user-supplied local range
    gkyl_create_ranges(&pkpm->low_inp.local_range, ghost, &app->local_ext, &app->local);
    
    if (pkpm->low_inp.comm)
      app->comm = gkyl_comm_acquire(pkpm->low_inp.comm);
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

  app->has_field = !pkpm->skip_field; // note inversion of truth value
  if (app->has_field)
    app->field = pkpm_field_new(pkpm, app);

  // allocate space to store species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct pkpm_species[ns])) : 0;

  // set info for each species: this needs to be done here as we need
  // to access species name from pkpm_species_init
  for (int i=0; i<ns; ++i)
    app->species[i].info = pkpm->species[i];

  // initialize each species
  for (int i=0; i<ns; ++i) 
    pkpm_species_init(pkpm, app, &app->species[i]);

  // initialize each species cross-species terms: this has to be done here
  // as need pointers to colliding species' collision objects
  // allocated in the previous step
  for (int i=0; i<ns; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS
      && app->species[i].lbo.num_cross_collisions) {
      pkpm_species_lbo_cross_init(app, &app->species[i], &app->species[i].lbo);
    }
  }
  
  // initialize stat object
  app->stat = (struct gkyl_pkpm_stat) {
    .use_gpu = app->use_gpu,
    .stage_2_dt_diff = { DBL_MAX, 0.0 },
    .stage_3_dt_diff = { DBL_MAX, 0.0 },
  };

  return app;
}

struct pkpm_species *
pkpm_find_species(const gkyl_pkpm_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return &app->species[i];
  return 0;
}

int
pkpm_find_species_idx(const gkyl_pkpm_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return i;
  return -1;
}

void
gkyl_pkpm_app_apply_ic(gkyl_pkpm_app* app, double t0)
{
  app->tcurr = t0;
  if (app->has_field) 
    gkyl_pkpm_app_apply_ic_field(app, t0);
  for (int i=0; i<app->num_species; ++i)
    gkyl_pkpm_app_apply_ic_species(app, i, t0);
}

void
gkyl_pkpm_app_apply_ic_field(gkyl_pkpm_app* app, double t0)
{
  app->tcurr = t0;

  struct timespec wtm = gkyl_wall_clock();
  pkpm_field_apply_ic(app, app->field, t0);
  app->stat.init_field_tm += gkyl_time_diff_now_sec(wtm);

  pkpm_field_apply_bc(app, app->field, app->field->em);
  // bvar is computed over extended range, so calculate it after
  // we apply BCs in the initialization step. Note that the apply_bc call
  // may not apply BCs in the corner cells and the corner values will just
  // be what comes from initializing the EM field over the extended range
  pkpm_field_calc_bvar(app, app->field, app->field->em); 
}

void
gkyl_pkpm_app_apply_ic_species(gkyl_pkpm_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  pkpm_species_apply_ic(app, &app->species[sidx], t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);

  pkpm_species_apply_bc(app, &app->species[sidx], app->species[sidx].f, app->species[sidx].fluid);
}

void
gkyl_pkpm_app_calc_integrated_mom(gkyl_pkpm_app* app, double tm)
{
  double avals[9], avals_global[9];

  struct timespec wst = gkyl_wall_clock();

  for (int i=0; i<app->num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    gkyl_array_clear(s->integ_pkpm_mom, 0.0);

    pkpm_species_calc_pkpm_vars(app, s, s->f, s->fluid);
    gkyl_dg_calc_pkpm_integrated_vars(s->calc_pkpm_vars, &app->local, 
      s->pkpm_moms.marr, s->fluid, 
      s->pkpm_prim, s->integ_pkpm_mom);
    gkyl_array_scale_range(s->integ_pkpm_mom, app->grid.cellVolume, &(app->local));
    if (app->use_gpu) {
      gkyl_array_reduce_range(s->red_integ_diag, s->integ_pkpm_mom, GKYL_SUM, &(app->local));
      gkyl_cu_memcpy(avals, s->red_integ_diag, sizeof(double[9]), GKYL_CU_MEMCPY_D2H);
    }
    else { 
      gkyl_array_reduce_range(avals, s->integ_pkpm_mom, GKYL_SUM, &(app->local));
    }

    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 9, avals, avals_global);
    gkyl_dynvec_append(s->integ_diag, tm, avals_global);
  }

  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += app->num_species;
}

void
gkyl_pkpm_app_calc_integrated_L2_f(gkyl_pkpm_app* app, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  for (int i=0; i<app->num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    pkpm_species_calc_L2(app, tm, s);
  }
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += app->num_species;
}

void
gkyl_pkpm_app_calc_field_energy(gkyl_pkpm_app* app, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  pkpm_field_calc_energy(app, tm, app->field);
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_pkpm_app_write(gkyl_pkpm_app* app, double tm, int frame)
{
  app->stat.nio += 1;
  struct timespec wtm = gkyl_wall_clock();
  
  if (app->has_field)
    gkyl_pkpm_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i) {
    gkyl_pkpm_app_write_species(app, i, tm, frame);
    gkyl_pkpm_app_write_mom(app, i, tm, frame);
  }

  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
}

void
gkyl_pkpm_app_write_field(gkyl_pkpm_app* app, double tm, int frame)
{
  const char *fmt = "%s-field_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->field->em_host, app->field->em);
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, 0, app->field->em_host, fileNm);
  }
  else {
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, 0, app->field->em, fileNm);
  }
}

void
gkyl_pkpm_app_write_species(gkyl_pkpm_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, app->species[sidx].info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[sidx].info.name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->species[sidx].f_host, app->species[sidx].f);
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, &app->species[sidx].local,
      0, app->species[sidx].f_host, fileNm);
  }
  else {
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, &app->species[sidx].local,
      0, app->species[sidx].f, fileNm);
  }
}

void
gkyl_pkpm_app_write_mom(gkyl_pkpm_app* app, int sidx, double tm, int frame)
{
  struct pkpm_species *s = &app->species[sidx];

  // Construct the file handles for the three quantities (PKPM moments, PKPM fluid variables, PKPM update variables)
  const char *fmt = "%s-%s_pkpm_moms_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, s->info.name, frame);

  const char *fmt_fluid = "%s-%s_pkpm_fluid_%d.gkyl";
  int sz_fluid = gkyl_calc_strlen(fmt_fluid, app->name, s->info.name, frame);
  char fileNm_fluid[sz_fluid+1]; // ensures no buffer overflow
  snprintf(fileNm_fluid, sizeof fileNm_fluid, fmt_fluid, app->name, s->info.name, frame);

  const char *fmt_pkpm_vars = "%s-%s_pkpm_vars_%d.gkyl";
  int sz_pkpm_vars = gkyl_calc_strlen(fmt_pkpm_vars, app->name, s->info.name, frame);
  char fileNm_pkpm_vars[sz_pkpm_vars+1]; // ensures no buffer overflow
  snprintf(fileNm_pkpm_vars, sizeof fileNm_pkpm_vars, fmt_pkpm_vars, app->name, s->info.name, frame);

  // Compute the PKPM variables including moments and primitive variables 
  // and construct arrays for writing out fluid and other pkpm variables.
  pkpm_species_moment_calc(&s->pkpm_moms_diag, s->local, app->local, s->f);
  pkpm_species_calc_pkpm_vars(app, s, s->f, s->fluid);
  pkpm_species_calc_pkpm_update_vars(app, s, s->f); 
  gkyl_dg_calc_pkpm_vars_io(s->calc_pkpm_vars, &app->local, 
    s->pkpm_moms.marr, s->fluid, 
    s->pkpm_p_ij, s->pkpm_prim, 
    s->pkpm_accel, s->fluid_io, s->pkpm_vars_io);

  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(s->pkpm_moms_diag.marr_host, s->pkpm_moms_diag.marr);
    gkyl_array_copy(s->fluid_io_host, s->fluid_io);
    gkyl_array_copy(s->pkpm_vars_io_host, s->pkpm_vars_io);
  }

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 0, s->pkpm_moms_diag.marr_host, fileNm);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 0, s->fluid_io_host, fileNm_fluid);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 0, s->pkpm_vars_io_host, fileNm_pkpm_vars);
}

void
gkyl_pkpm_app_write_integrated_mom(gkyl_pkpm_app *app)
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
gkyl_pkpm_app_write_integrated_L2_f(gkyl_pkpm_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out integrated L^2
      const char *fmt = "%s-%s-%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, app->species[i].info.name,
        "L2");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[i].info.name,
        "L2");

      if (app->species[i].is_first_integ_L2_write_call) {
        // write to a new file (this ensure previous output is removed)
        gkyl_dynvec_write(app->species[i].integ_L2_f, fileNm);
        app->species[i].is_first_integ_L2_write_call = false;
      }
      else {
        // append to existing file
        gkyl_dynvec_awrite(app->species[i].integ_L2_f, fileNm);
      }
    }
    gkyl_dynvec_clear(app->species[i].integ_L2_f);
  }
}

void
gkyl_pkpm_app_write_field_energy(gkyl_pkpm_app* app)
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
forward_euler(gkyl_pkpm_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[], const struct gkyl_array *emin,
  struct gkyl_array *fout[], struct gkyl_array *fluidout[], struct gkyl_array *emout, struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;

  double dtmin = DBL_MAX;

  // Compute external EM field or applied currents if present.
  // Note: external EM field and  applied currents use proj_on_basis 
  // so does copy to GPU every call if app->use_gpu = true.
  if (app->has_field) {
    if (app->field->ext_em_evolve)
      pkpm_field_calc_ext_em(app, app->field, tcurr);
    if (app->field->app_current_evolve)
      pkpm_field_calc_app_current(app, app->field, tcurr); 
    // Also compute magnetic field unit vector and tensor,
    // div(b) and div(b) at quadrature points for update
    pkpm_field_calc_bvar(app, app->field, emin);  
  }

  // Two separate loops over number of species to compute needed auxiliary quantities.
  for (int i=0; i<app->num_species; ++i) {
    // Compute parallel-kinetic-perpendicular moment (pkpm) variables.
    // These are the coupling moments [rho, p_par, p_perp], the self-consistent
    // pressure force (div(p_par b_hat)), and the primitive variables
    pkpm_species_calc_pkpm_vars(app, &app->species[i], fin[i], fluidin[i]);
    // compute necessary moments and boundary corrections for collisions
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      pkpm_species_lbo_moms(app, &app->species[i], &app->species[i].lbo, fin[i]);
    }
  }
  for (int i=0; i<app->num_species; ++i) {
    // compute necessary moments for cross-species collisions
    // needs to be done after self-collisions moments, so separate loop over species
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS
      && app->species[i].lbo.num_cross_collisions) {
      pkpm_species_lbo_cross_moms(app, &app->species[i], &app->species[i].lbo, fin[i]);
    }
    // Finish computing parallel-kinetic-perpendicular moment (pkpm) variables.
    // These are the update variables including the acceleration variables in the kinetic
    // equation and the source distribution functions for Laguerre couplings.
    // Needs to be done after all collisional moment computations for collisional sources
    // in Laguerre couplings.
    pkpm_species_calc_pkpm_update_vars(app, &app->species[i], fin[i]); 
  }

  // compute RHS of pkpm equations
  for (int i=0; i<app->num_species; ++i) {
    double dt1 = pkpm_species_rhs(app, &app->species[i], fin[i], fluidin[i], emin, fout[i], fluidout[i]);
    dtmin = fmin(dtmin, dt1);
  }

  // compute RHS of Maxwell equations
  if (app->has_field) {
    double dt1 = pkpm_field_rhs(app, app->field, emin, emout);
    dtmin = fmin(dtmin, dt1);
  }

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // compute minimum time-step across all processors
  double dtmin_local = dtmin, dtmin_global;
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dtmin_local, &dtmin_global);
  dtmin = dtmin_global;
  
  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // complete update of distribution function
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout[i], dta), 1.0, fin[i]);
    gkyl_array_accumulate(gkyl_array_scale(fluidout[i], dta), 1.0, fluidin[i]);
    pkpm_species_apply_bc(app, &app->species[i], fout[i], fluidout[i]);
  }

  if (app->has_field) {
    struct timespec wst = gkyl_wall_clock();

    // (can't accumulate current when field is static)
    if (!app->field->info.is_static) {
      // accumulate current contribution from kinetic species to electric field terms
      pkpm_field_accumulate_current(app, fluidin, emout);
      app->stat.current_tm += gkyl_time_diff_now_sec(wst);
    }

    // complete update of field (even when field is static, it is
    // safest to do this accumulate as it ensure emout = emin)
    gkyl_array_accumulate(gkyl_array_scale(emout, dta), 1.0, emin);

    pkpm_field_apply_bc(app, app->field, emout);
  }
}

// Take time-step using the RK3 method. Also sets the status object
// which has the actual and suggested dts used. These can be different
// from the actual time-step.
static struct gkyl_update_status
rk3(gkyl_pkpm_app* app, double dt0)
{
  const struct gkyl_array *fin[app->num_species];
  struct gkyl_array *fout[app->num_species];
  const struct gkyl_array *fluidin[app->num_species];
  struct gkyl_array *fluidout[app->num_species];
  struct gkyl_update_status st = { .success = true };

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f;
          fluidin[i] = app->species[i].fluid;
          fout[i] = app->species[i].f1;
          fluidout[i] = app->species[i].fluid1;
        }
        forward_euler(app, tcurr, dt, fin, fluidin, app->has_field ? app->field->em : 0,
          fout, fluidout, app->has_field ? app->field->em1 : 0,
          &st
        );
        dt = st.dt_actual;
        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fluidin[i] = app->species[i].fluid1;
          fout[i] = app->species[i].fnew;
          fluidout[i] = app->species[i].fluidnew;
        }
        forward_euler(app, tcurr+dt, dt, fin, fluidin, app->has_field ? app->field->em1 : 0,
          fout, fluidout, app->has_field ? app->field->emnew : 0,
          &st
        );
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

        } else {
          for (int i=0; i<app->num_species; ++i)
            array_combine(app->species[i].f1,
              3.0/4.0, app->species[i].f, 1.0/4.0, app->species[i].fnew, &app->species[i].local_ext);
          for (int i=0; i<app->num_species; ++i)
            array_combine(app->species[i].fluid1,
              3.0/4.0, app->species[i].fluid, 1.0/4.0, app->species[i].fluidnew, &app->local_ext);
          if (app->has_field)
            array_combine(app->field->em1,
              3.0/4.0, app->field->em, 1.0/4.0, app->field->emnew, &app->local_ext);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fluidin[i] = app->species[i].fluid1;
          fout[i] = app->species[i].fnew;
          fluidout[i] = app->species[i].fluidnew;
        }
        forward_euler(app, tcurr+dt/2, dt, fin, fluidin, app->has_field ? app->field->em1 : 0,
          fout, fluidout, app->has_field ? app->field->emnew : 0,
          &st
        );
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
          for (int i=0; i<app->num_species; ++i) {
            array_combine(app->species[i].fluid1,
              1.0/3.0, app->species[i].fluid, 2.0/3.0, app->species[i].fluidnew, &app->local_ext);
            gkyl_array_copy_range(app->species[i].fluid, app->species[i].fluid1, &app->local_ext);
          }
          if (app->has_field) {
            array_combine(app->field->em1,
              1.0/3.0, app->field->em, 2.0/3.0, app->field->emnew, &app->local_ext);
            gkyl_array_copy_range(app->field->em, app->field->em1, &app->local_ext);
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
gkyl_pkpm_update(gkyl_pkpm_app* app, double dt)
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

struct gkyl_pkpm_stat
gkyl_pkpm_app_stat(gkyl_pkpm_app* app)
{
  pkpm_species_tm(app);
  pkpm_species_coll_tm(app);
  return app->stat;
}

static void
range_stat_write(gkyl_pkpm_app* app, const char *nm, const struct gkyl_range *r, FILE *fp)
{
  gkyl_pkpm_app_cout(app, fp, " %s_cells : [ ", nm);
  for (int i=0; i<r->ndim; ++i)
    gkyl_pkpm_app_cout(app, fp, " %d, ", gkyl_range_shape(r, i));
  gkyl_pkpm_app_cout(app, fp, " ],\n");
}

// ensure stats across processors are made consistent
static void
comm_reduce_app_stat(const gkyl_pkpm_app* app,
  const struct gkyl_pkpm_stat *local, struct gkyl_pkpm_stat *global)
{
  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);
  if (comm_sz == 1) {
    memcpy(global, local, sizeof(struct gkyl_pkpm_stat));
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
    TOTAL_TM, INIT_SPECIES_TM, INIT_FLUID_SPECIES_TM, INIT_FIELD_TM,
    SPECIES_RHS_TM, FLUID_SPECIES_RHS_TM, SPECIES_COLL_MOM_TM,
    SPECIES_COL_TM, SPECIES_PKPM_VARS_TM, FIELD_RHS_TM, FIELD_EM_VARS_TM, CURRENT_TM,
    SPECIES_OMEGA_CFL_TM, FIELD_OMEGA_CFL_TM, DIAG_TM, IO_TM,
    SPECIES_BC_TM, FLUID_SPECIES_BC_TM, FIELD_BC_TM,
    D_END
  };

  double d_red[D_END] = {
    [TOTAL_TM] = local->total_tm,
    [INIT_SPECIES_TM] = local->init_species_tm,
    [INIT_FLUID_SPECIES_TM] = local->init_fluid_species_tm,
    [INIT_FIELD_TM] = local->field_rhs_tm,
    [SPECIES_RHS_TM] = local->species_rhs_tm,
    [FLUID_SPECIES_RHS_TM] = local->fluid_species_rhs_tm,
    [SPECIES_COLL_MOM_TM] = local->species_coll_mom_tm,
    [SPECIES_COL_TM] = local->species_coll_tm,
    [SPECIES_PKPM_VARS_TM] = local->species_pkpm_vars_tm,
    [FIELD_RHS_TM] = local->field_rhs_tm,
    [FIELD_EM_VARS_TM] = local->field_em_vars_tm,
    [CURRENT_TM] = local->current_tm,
    [SPECIES_OMEGA_CFL_TM] = local->species_omega_cfl_tm,
    [FIELD_OMEGA_CFL_TM] = local->field_omega_cfl_tm,
    [DIAG_TM] = local->diag_tm,
    [IO_TM] = local->io_tm,
    [SPECIES_BC_TM] = local->species_bc_tm,
    [FLUID_SPECIES_BC_TM] = local->fluid_species_bc_tm,
    [FIELD_BC_TM] = local->field_bc_tm
  };

  double d_red_global[D_END];
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, D_END, d_red, d_red_global);
  
  global->total_tm = d_red_global[TOTAL_TM];
  global->init_species_tm = d_red_global[INIT_SPECIES_TM];
  global->init_fluid_species_tm = d_red_global[INIT_FLUID_SPECIES_TM];
  global->field_rhs_tm = d_red_global[INIT_FIELD_TM];
  global->species_rhs_tm = d_red_global[SPECIES_RHS_TM];
  global->fluid_species_rhs_tm = d_red_global[FLUID_SPECIES_RHS_TM];
  global->species_coll_mom_tm = d_red_global[SPECIES_COLL_MOM_TM];
  global->species_coll_tm = d_red_global[SPECIES_COL_TM];
  global->species_pkpm_vars_tm = d_red_global[SPECIES_PKPM_VARS_TM];
  global->field_rhs_tm = d_red_global[FIELD_RHS_TM];
  global->field_em_vars_tm = d_red_global[FIELD_EM_VARS_TM];
  global->current_tm = d_red_global[CURRENT_TM];
  global->species_omega_cfl_tm = d_red_global[SPECIES_OMEGA_CFL_TM];
  global->field_omega_cfl_tm = d_red_global[FIELD_OMEGA_CFL_TM];
  global->diag_tm = d_red_global[DIAG_TM];
  global->io_tm = d_red_global[IO_TM];
  global->species_bc_tm = d_red_global[SPECIES_BC_TM];
  global->fluid_species_bc_tm = d_red_global[FLUID_SPECIES_BC_TM];
  global->field_bc_tm = d_red_global[FIELD_BC_TM];

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
gkyl_pkpm_app_stat_write(gkyl_pkpm_app* app)
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

  pkpm_species_coll_tm(app);
  pkpm_species_tm(app);

  struct gkyl_pkpm_stat stat = { };
  comm_reduce_app_stat(app, &app->stat, &stat);
  
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  if (rank == 0) fp = fopen(fileNm, "a");

  gkyl_pkpm_app_cout(app, fp, "{\n");

  if (strftime(buff, sizeof buff, "%c", &curr_tm))
    gkyl_pkpm_app_cout(app, fp, " date : %s,\n", buff);

  gkyl_pkpm_app_cout(app, fp, " use_gpu : %d,\n", stat.use_gpu);
  gkyl_pkpm_app_cout(app, fp, " num_ranks : %d,\n", num_ranks); 
  
  for (int s=0; s<app->num_species; ++s)
    range_stat_write(app, app->species[s].info.name, &app->species[s].global, fp);
  
  gkyl_pkpm_app_cout(app, fp, " nup : %ld,\n", stat.nup);
  gkyl_pkpm_app_cout(app, fp, " nfeuler : %ld,\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, fp, " nstage_2_fail : %ld,\n", stat.nstage_2_fail);
  gkyl_pkpm_app_cout(app, fp, " nstage_3_fail : %ld,\n", stat.nstage_3_fail);

  gkyl_pkpm_app_cout(app, fp, " stage_2_dt_diff : [ %lg, %lg ],\n",
    stat.stage_2_dt_diff[0], stat.stage_2_dt_diff[1]);
  gkyl_pkpm_app_cout(app, fp, " stage_3_dt_diff : [ %lg, %lg ],\n",
    stat.stage_3_dt_diff[0], stat.stage_3_dt_diff[1]);

  gkyl_pkpm_app_cout(app, fp, " total_tm : %lg,\n", stat.total_tm);
  gkyl_pkpm_app_cout(app, fp, " init_species_tm : %lg,\n", stat.init_species_tm);
  if (app->has_field)
    gkyl_pkpm_app_cout(app, fp, " init_field_tm : %lg,\n", stat.init_field_tm);
  
  gkyl_pkpm_app_cout(app, fp, " species_rhs_tm : %lg,\n", stat.species_rhs_tm);

  for (int s=0; s<app->num_species; ++s) {
    gkyl_pkpm_app_cout(app, fp, " species_coll_drag_tm[%d] : %lg,\n", s,
      stat.species_lbo_coll_drag_tm[s]);
    gkyl_pkpm_app_cout(app, fp, " species_coll_diff_tm[%d] : %lg,\n", s,
      stat.species_lbo_coll_diff_tm[s]);
  }

  gkyl_pkpm_app_cout(app, fp, " species_coll_mom_tm : %lg,\n", stat.species_coll_mom_tm);
  gkyl_pkpm_app_cout(app, fp, " species_coll_tm : %lg,\n", stat.species_coll_tm);

  gkyl_pkpm_app_cout(app, fp, " species_pkpm_vars_tm : %lg,\n", stat.species_pkpm_vars_tm);

  gkyl_pkpm_app_cout(app, fp, " species_bc_tm : %lg,\n", stat.species_bc_tm);
  
  gkyl_pkpm_app_cout(app, fp, " fluid_species_rhs_tm : %lg,\n", stat.fluid_species_rhs_tm);

  gkyl_pkpm_app_cout(app, fp, " fluid_species_bc_tm : %lg,\n", stat.fluid_species_bc_tm);

  if (app->has_field) {
    gkyl_pkpm_app_cout(app, fp, " field_rhs_tm : %lg,\n", stat.field_rhs_tm);
    gkyl_pkpm_app_cout(app, fp, " field_bc_tm : %lg,\n", stat.field_bc_tm);
    
    gkyl_pkpm_app_cout(app, fp, " field_em_vars_tm : %lg,\n", stat.field_em_vars_tm);
    gkyl_pkpm_app_cout(app, fp, " current_tm : %lg,\n", stat.current_tm);
  }

  gkyl_pkpm_app_cout(app, fp, " ndiag : %ld,\n", stat.ndiag);
  gkyl_pkpm_app_cout(app, fp, " diag_tm : %lg\n", stat.diag_tm);
  
  gkyl_pkpm_app_cout(app, fp, " nspecies_omega_cfl : %ld,\n", stat.nspecies_omega_cfl);
  gkyl_pkpm_app_cout(app, fp, " species_omega_cfl_tm : %lg\n", stat.species_omega_cfl_tm);

  gkyl_pkpm_app_cout(app, fp, " nfield_omega_cfl : %ld,\n", stat.nfield_omega_cfl);
  gkyl_pkpm_app_cout(app, fp, " field_omega_cfl_tm : %lg\n", stat.field_omega_cfl_tm);

  gkyl_pkpm_app_cout(app, fp, " nio : %ld,\n", stat.nio);
  gkyl_pkpm_app_cout(app, fp, " io_tm : %lg\n", stat.io_tm);
  
  gkyl_pkpm_app_cout(app, fp, "}\n");

  if (rank == 0)
    fclose(fp);  

}

// private function to handle variable argument list for printing
static void
v_pkpm_app_cout(const gkyl_pkpm_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank, r = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp)
    vfprintf(fp, fmt, argp);
}

void
gkyl_pkpm_app_cout(const gkyl_pkpm_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_pkpm_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_pkpm_app_release(gkyl_pkpm_app* app)
{
  for (int i=0; i<app->num_species; ++i)
    pkpm_species_release(app, &app->species[i]);
  if (app->num_species > 0)
    gkyl_free(app->species);
  if (app->has_field)
    pkpm_field_release(app, app->field);

  gkyl_comm_release(app->comm);

  if (app->use_gpu) {
    gkyl_cu_free(app->basis_on_dev.basis);
    gkyl_cu_free(app->basis_on_dev.confBasis);
  }

  gkyl_free(app);
}
