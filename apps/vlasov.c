#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>

#include <gkyl_vlasov_priv.h>

gkyl_vlasov_app*
gkyl_vlasov_app_new(struct gkyl_vm *vm)
{
  disable_denorm_float();
  
  assert(vm->num_species <= GKYL_MAX_SPECIES);
  
  gkyl_vlasov_app *app = gkyl_malloc(sizeof(gkyl_vlasov_app));
  
  int cdim = app->cdim = vm->cdim;
  int vdim = app->vdim = vm->vdim;
  int pdim = cdim+vdim;
  int poly_order = app->poly_order = vm->poly_order;
  int ns = app->num_species = vm->num_species;
  int nsf = app->num_fluid_species = vm->num_fluid_species;

  double cfl_frac = vm->cfl_frac == 0 ? 1.0 : vm->cfl_frac;
  app->cfl = cfl_frac/(2*poly_order+1);

#ifdef GKYL_HAVE_CUDA
  app->use_gpu = vm->use_gpu;
#else
  app->use_gpu = false; // can't use GPUs if we don't have them!
#endif
  
  app->num_periodic_dir = vm->num_periodic_dir;
  for (int d=0; d<cdim; ++d)
    app->periodic_dirs[d] = vm->periodic_dirs[d];

  strcpy(app->name, vm->name);
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
  switch (vm->basis_type) {
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

  gkyl_rect_grid_init(&app->grid, cdim, vm->lower, vm->upper, vm->cells);

  int ghost[] = { 1, 1, 1 };  
  gkyl_create_grid_ranges(&app->grid, ghost, &app->local_ext, &app->local);
  skin_ghost_ranges_init(&app->skin_ghost, &app->local_ext, ghost);

  app->has_field = !vm->skip_field; // note inversion of truth value

  if (app->has_field)
    app->field = vm_field_new(vm, app);

  // allocate space to store species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct vm_species[ns])) : 0;
  
  // set info for each species: this needs to be done here as we need
  // to access species name from vm_species_init
  for (int i=0; i<ns; ++i)
    app->species[i].info = vm->species[i];

  // allocate space to store fluid species objects
  app->fluid_species = nsf>0 ? gkyl_malloc(sizeof(struct vm_fluid_species[nsf])) : 0;

  // set info for each fluid species: this needs to be done here as we
  // need to access species name from vm_fluid_species_init
  for (int i=0; i<nsf; ++i)
    app->fluid_species[i].info = vm->fluid_species[i];  
  
  // initialize each species
  for (int i=0; i<ns; ++i)
    vm_species_init(vm, app, &app->species[i]);

  // initialize each species cross-species terms: this has to be done here
  // as need pointers to colliding species' collision objects
  // allocated in the previous step
  for (int i=0; i<ns; ++i)
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS
      && app->species[i].lbo.num_cross_collisions) {
      vm_species_lbo_cross_init(app, &app->species[i], &app->species[i].lbo);
    }

  // initialize each species source terms: this has to be done here
  // as they may initialize a bflux updater for their source species
  for (int i=0; i<ns; ++i)
    if (app->species[i].source_id)
      vm_species_source_init(app, &app->species[i], &app->species[i].src);
  
  // initialize each fluid species
  for (int i=0; i<nsf; ++i)
    vm_fluid_species_init(vm, app, &app->fluid_species[i]);
  
  // initialize stat object
  app->stat = (struct gkyl_vlasov_stat) {
    .use_gpu = app->use_gpu,
    .stage_2_dt_diff = { DBL_MAX, 0.0 },
    .stage_3_dt_diff = { DBL_MAX, 0.0 },
  };
  
  return app;
}

struct vm_species *
vm_find_species(const gkyl_vlasov_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return &app->species[i];
  return 0;
}

int
vm_find_species_idx(const gkyl_vlasov_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return i;
  return -1;
}

struct vm_fluid_species *
vm_find_fluid_species(const gkyl_vlasov_app *app, const char *nm)
{
  for (int i=0; i<app->num_fluid_species; ++i)
    if (strcmp(nm, app->fluid_species[i].info.name) == 0)
      return &app->fluid_species[i];
  return 0;
}

int
vm_find_fluid_species_idx(const gkyl_vlasov_app *app, const char *nm)
{
  for (int i=0; i<app->num_fluid_species; ++i)
    if (strcmp(nm, app->fluid_species[i].info.name) == 0)
      return i;
  return -1;  
}

void
gkyl_vlasov_app_apply_ic(gkyl_vlasov_app* app, double t0)
{
  app->tcurr = t0;
  if (app->has_field)
    gkyl_vlasov_app_apply_ic_field(app, t0);
  for (int i=0; i<app->num_species; ++i)
    gkyl_vlasov_app_apply_ic_species(app, i, t0);
  for (int i=0; i<app->num_fluid_species; ++i)
    gkyl_vlasov_app_apply_ic_fluid_species(app, i, t0);
}

void
gkyl_vlasov_app_apply_ic_field(gkyl_vlasov_app* app, double t0)
{
  app->tcurr = t0;

  struct timespec wtm = gkyl_wall_clock();
  vm_field_apply_ic(app, app->field, t0);
  app->stat.init_field_tm += gkyl_time_diff_now_sec(wtm);
    
  vm_field_apply_bc(app, app->field, app->field->em);  
}

void
gkyl_vlasov_app_apply_ic_species(gkyl_vlasov_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  vm_species_apply_ic(app, &app->species[sidx], t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);

  vm_species_apply_bc(app, &app->species[sidx], app->species[sidx].f);
}

void
gkyl_vlasov_app_apply_ic_fluid_species(gkyl_vlasov_app* app, int sidx, double t0)
{
  assert(sidx < app->num_fluid_species);

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  vm_fluid_species_apply_ic(app, &app->fluid_species[sidx], t0);
  app->stat.init_fluid_species_tm += gkyl_time_diff_now_sec(wtm);

  vm_fluid_species_apply_bc(app, &app->fluid_species[sidx], app->fluid_species[sidx].fluid);
}

void
gkyl_vlasov_app_calc_mom(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *s = &app->species[i];
    
    for (int m=0; m<app->species[i].info.num_diag_moments; ++m) {
      struct timespec wst = gkyl_wall_clock();
      vm_species_moment_calc(&s->moms[m], s->local, app->local, s->f);
      app->stat.mom_tm += gkyl_time_diff_now_sec(wst);
      app->stat.nmom += 1;
    }
  }
}

void
gkyl_vlasov_app_calc_integrated_mom(gkyl_vlasov_app* app, double tm)
{
  double avals[2+GKYL_MAX_DIM];

  struct timespec wst = gkyl_wall_clock();
  
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *s = &app->species[i];
    
    struct timespec wst = gkyl_wall_clock();
    vm_species_moment_calc(&s->integ_moms, s->local, app->local, s->f);

    // reduce to compute sum over whole domain, append to diagnostics
    if (app->use_gpu) {
      gkyl_array_reduce_range(s->red_integ_diag, s->integ_moms.marr, GKYL_SUM, app->local);
      gkyl_cu_memcpy(avals, s->red_integ_diag, sizeof(double[2+GKYL_MAX_DIM]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      gkyl_array_reduce_range(avals, s->integ_moms.marr_host, GKYL_SUM, app->local);
    }
    gkyl_dynvec_append(s->integ_diag, tm, avals);
    
    app->stat.mom_tm += gkyl_time_diff_now_sec(wst);
    app->stat.nmom += 1;
  }

  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;  
}

void
gkyl_vlasov_app_calc_field_energy(gkyl_vlasov_app* app, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  vm_field_calc_energy(app, tm, app->field, app->field->em);
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_vlasov_app_write(gkyl_vlasov_app* app, double tm, int frame)
{
  if (app->has_field)
    gkyl_vlasov_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i) {
    gkyl_vlasov_app_write_species(app, i, tm, frame);
    if (app->species[i].field_id == GKYL_FIELD_SR_E_B && frame == 0)
      gkyl_vlasov_app_write_species_gamma(app, i, tm, frame);
    if (app->species[i].has_mirror_force && frame == 0) {
      gkyl_vlasov_app_write_magB(app, i, tm, frame);
      gkyl_vlasov_app_write_gradB(app, i, tm, frame);
    }
  }
  for (int i=0; i<app->num_fluid_species; ++i)
    gkyl_vlasov_app_write_fluid_species(app, i, tm, frame);
}
 
void
gkyl_vlasov_app_write_field(gkyl_vlasov_app* app, double tm, int frame)
{
  const char *fmt = "%s-field_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out    
    gkyl_array_copy(app->field->em_host, app->field->em);
    gkyl_grid_sub_array_write(&app->grid, &app->local, app->field->em_host, fileNm);
  }
  else {
    gkyl_grid_sub_array_write(&app->grid, &app->local, app->field->em, fileNm);
  }
}

void
gkyl_vlasov_app_write_species(gkyl_vlasov_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, app->species[sidx].info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[sidx].info.name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->species[sidx].f_host, app->species[sidx].f);
    gkyl_grid_sub_array_write(&app->species[sidx].grid, &app->species[sidx].local,
      app->species[sidx].f_host, fileNm);
  }
  else {
    gkyl_grid_sub_array_write(&app->species[sidx].grid, &app->species[sidx].local,
      app->species[sidx].f, fileNm);
  }
}

void
gkyl_vlasov_app_write_species_gamma(gkyl_vlasov_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_p_over_gamma_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, app->species[sidx].info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[sidx].info.name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->species[sidx].p_over_gamma_host, app->species[sidx].p_over_gamma);
    gkyl_grid_sub_array_write(&app->species[sidx].grid_vel, &app->species[sidx].local_vel,
      app->species[sidx].p_over_gamma_host, fileNm);
  }
  else {
    gkyl_grid_sub_array_write(&app->species[sidx].grid_vel, &app->species[sidx].local_vel,
      app->species[sidx].p_over_gamma, fileNm);
  }  
}

void
gkyl_vlasov_app_write_magB(gkyl_vlasov_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_magB_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, app->species[sidx].info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[sidx].info.name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->species[sidx].magB_host, app->species[sidx].magB);
    gkyl_grid_sub_array_write(&app->grid, &app->local,
      app->species[sidx].magB_host, fileNm);
  }
  else {
    gkyl_grid_sub_array_write(&app->grid, &app->local,
      app->species[sidx].magB, fileNm);
  }  
}

void
gkyl_vlasov_app_write_gradB(gkyl_vlasov_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_gradB_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, app->species[sidx].info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[sidx].info.name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->species[sidx].gradB_host, app->species[sidx].gradB);
    gkyl_grid_sub_array_write(&app->grid, &app->local,
      app->species[sidx].gradB_host, fileNm);
  }
  else {
    gkyl_grid_sub_array_write(&app->grid, &app->local,
      app->species[sidx].gradB, fileNm);
  }  
}

void
gkyl_vlasov_app_write_fluid_species(gkyl_vlasov_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, app->fluid_species[sidx].info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->fluid_species[sidx].info.name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->fluid_species[sidx].fluid_host, app->fluid_species[sidx].fluid);
    gkyl_grid_sub_array_write(&app->grid, &app->local,
      app->fluid_species[sidx].fluid_host, fileNm);
  }
  else {
    gkyl_grid_sub_array_write(&app->grid, &app->local,
      app->fluid_species[sidx].fluid, fileNm);
  }
}

void
gkyl_vlasov_app_write_mom(gkyl_vlasov_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    for (int m=0; m<app->species[i].info.num_diag_moments; ++m) {
      
      const char *fmt = "%s-%s-%s_%d.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, app->species[i].info.name,
        app->species[i].info.diag_moments[m], frame);
      char fileNm[sz+1]; // ensures no buffer overflow  
      snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[i].info.name,
        app->species[i].info.diag_moments[m], frame);

      if (app->use_gpu)
        gkyl_array_copy(app->species[i].moms[m].marr_host, app->species[i].moms[m].marr);
      
      gkyl_grid_sub_array_write(&app->grid, &app->local, app->species[i].moms[m].marr_host, fileNm);
    }
  }
}

void
gkyl_vlasov_app_write_integrated_mom(gkyl_vlasov_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    // write out diagnostic moments
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
    gkyl_dynvec_clear(app->species[i].integ_diag);
  }
}

void
gkyl_vlasov_app_write_field_energy(gkyl_vlasov_app* app)
{
  // write out diagnostic moments
  const char *fmt = "%s-field-energy.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name);

  if (app->field->is_first_energy_write_call) {
    // write to a new file (this ensure previous output is removed)
    gkyl_dynvec_write(app->field->integ_energy, fileNm);
    app->field->is_first_energy_write_call = false;
  }
  else {
    // append to existing file
    gkyl_dynvec_awrite(app->field->integ_energy, fileNm);
  }
  gkyl_dynvec_clear(app->field->integ_energy);
}

// Take a forward Euler step with the suggested time-step dt. This may
// not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by
// stability. The actual time-step and dt_suggested are returned in
// the status object.
static void
forward_euler(gkyl_vlasov_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[], const struct gkyl_array *emin,
  struct gkyl_array *fout[], struct gkyl_array *fluidout[], struct gkyl_array *emout, struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;
  
  double dtmin = DBL_MAX;

  // compute necessary moments and boundary corrections for collisions
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      // Pass full fluidin array to be agnostic to species ordering
      // vm_species_lbo_moms will find the particular fluid species being collided with if it exists
      vm_species_lbo_moms(app, &app->species[i], &app->species[i].lbo, fin[i], app->species[i].collides_with_fluid, fluidin);
    }
  }

  // compute necessary moments for cross-species collisions
  // needs to be done after self-collisions moments, so separate loop over species
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS
      && app->species[i].lbo.num_cross_collisions) {
      // Pass full fluidin array to be agnostic to species ordering
      // vm_species_lbo_moms will find the particular fluid species being collided with if it exists
      vm_species_lbo_cross_moms(app, &app->species[i], &app->species[i].lbo, fin[i], app->species[i].collides_with_fluid, fluidin);
    }
  }
  
  // compute RHS of Vlasov equations
  for (int i=0; i<app->num_species; ++i) {
    double dt1 = vm_species_rhs(app, &app->species[i], fin[i], emin, fout[i], fluidin);
    dtmin = fmin(dtmin, dt1);
  }
  for (int i=0; i<app->num_fluid_species; ++i) {
    double dt1 = vm_fluid_species_rhs(app, &app->fluid_species[i], fluidin[i], fluidout[i]);
    dtmin = fmin(dtmin, dt1);
  }
  // compute source term
  // done here as the RHS update for all species should be complete before
  // bflux calculation of the source species
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].source_id) {
      vm_species_source_rhs(app, &app->species[i], &app->species[i].src, fin, fout);
    }
  }
  // compute RHS of Maxwell equations
  if (app->has_field) {
    double dt1 = vm_field_rhs(app, app->field, emin, emout);
    dtmin = fmin(dtmin, dt1);
  }

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // complete update of distribution function
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate_range(gkyl_array_scale_range(fout[i], dta, app->species[i].local),
      1.0, fin[i], app->species[i].local);
    vm_species_apply_bc(app, &app->species[i], fout[i]);
  }

  // complete update of fluid species
  for (int i=0; i<app->num_fluid_species; ++i) {
    gkyl_array_accumulate_range(gkyl_array_scale_range(fluidout[i], dta, app->local),
      1.0, fluidin[i], app->local);
    vm_fluid_species_apply_bc(app, &app->fluid_species[i], fluidout[i]);
  }

  if (app->has_field) {
    struct timespec wst = gkyl_wall_clock();

    // (can't accumulate current when field is static)
    if (!app->field->info.is_static) {
      // accumulate current contribution from kinetic species to electric field terms
      for (int i=0; i<app->num_species; ++i) {
        struct vm_species *s = &app->species[i];
        vm_species_moment_calc(&s->m1i, s->local, app->local, fin[i]);
    
        double qbyeps = s->info.charge/app->field->info.epsilon0;
        if (s->has_mirror_force) {
	  // if has_mirror_force, m1i = m1i/B = J*m1i
	  // need to multiply by magB to get the right current density
          gkyl_dg_mul_op_range(app->confBasis, 0, s->m1i_no_J, 0,
            s->m1i.marr, 0, s->magB, &app->local);
	  gkyl_array_accumulate_range(emout, -qbyeps, s->m1i_no_J, app->local);
        } else {
          gkyl_array_accumulate_range(emout, -qbyeps, s->m1i.marr, app->local);
	}
      }
      // accumulate current contribution from fluid species to electric field terms
      for (int i=0; i<app->num_fluid_species; ++i) {
        struct vm_fluid_species *f = &app->fluid_species[i];

        // when charge is specified in fluid species, compute qbyeps, otherwise assume zero
        double qbyeps = f->info.charge ? f->info.charge/app->field->info.epsilon0 : 0.0;
        gkyl_array_accumulate_range(emout, -qbyeps, f->u, app->local);
      }
      app->stat.current_tm += gkyl_time_diff_now_sec(wst);
    }
  
    // complete update of field (even when field is static, it is
    // safest to do this accumulate as it ensure emout = emin)
    gkyl_array_accumulate_range(gkyl_array_scale_range(emout, dta, app->local),
      1.0, emin, app->local);
    
    vm_field_apply_bc(app, app->field, emout);
  }
}

// Take time-step using the RK3 method. Also sets the status object
// which has the actual and suggested dts used. These can be different
// from the actual time-step.
static struct gkyl_update_status
rk3(gkyl_vlasov_app* app, double dt0)
{
  const struct gkyl_array *fin[app->num_species];
  struct gkyl_array *fout[app->num_species];
  const struct gkyl_array *fluidin[app->num_fluid_species];
  struct gkyl_array *fluidout[app->num_fluid_species];
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
        for (int i=0; i<app->num_fluid_species; ++i) {
          fluidin[i] = app->fluid_species[i].fluid;
          fluidout[i] = app->fluid_species[i].fluid1;
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
          fout[i] = app->species[i].fnew;
        }
        for (int i=0; i<app->num_fluid_species; ++i) {
          fluidin[i] = app->fluid_species[i].fluid1;
          fluidout[i] = app->fluid_species[i].fluidnew;
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
              3.0/4.0, app->species[i].f, 1.0/4.0, app->species[i].fnew, app->species[i].local_ext);
          for (int i=0; i<app->num_fluid_species; ++i)
            array_combine(app->fluid_species[i].fluid1,
              3.0/4.0, app->fluid_species[i].fluid, 1.0/4.0, app->fluid_species[i].fluidnew, app->local_ext);
          if (app->has_field)
            array_combine(app->field->em1,
              3.0/4.0, app->field->em, 1.0/4.0, app->field->emnew, app->local_ext);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
        }
        for (int i=0; i<app->num_fluid_species; ++i) {
          fluidin[i] = app->fluid_species[i].fluid1;
          fluidout[i] = app->fluid_species[i].fluidnew;
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
              1.0/3.0, app->species[i].f, 2.0/3.0, app->species[i].fnew, app->species[i].local_ext);
            gkyl_array_copy_range(app->species[i].f, app->species[i].f1, app->species[i].local_ext);
          }
          for (int i=0; i<app->num_fluid_species; ++i) {
            array_combine(app->fluid_species[i].fluid1,
              1.0/3.0, app->fluid_species[i].fluid, 2.0/3.0, app->fluid_species[i].fluidnew, app->local_ext);
            gkyl_array_copy_range(app->fluid_species[i].fluid, app->fluid_species[i].fluid1, app->local_ext);
          }
          if (app->has_field) {
            array_combine(app->field->em1,
              1.0/3.0, app->field->em, 2.0/3.0, app->field->emnew, app->local_ext);
            gkyl_array_copy_range(app->field->em, app->field->em1, app->local_ext);
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
gkyl_vlasov_update(gkyl_vlasov_app* app, double dt)
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

struct gkyl_vlasov_stat
gkyl_vlasov_app_stat(gkyl_vlasov_app* app)
{
  vm_species_tm(app);
  vm_species_coll_tm(app);
  return app->stat;
}

void
gkyl_vlasov_app_species_ktm_rhs(gkyl_vlasov_app* app, int update_vol_term)
{
  for (int i=0; i<app->num_species; ++i) {
    
    struct vm_species *species = &app->species[i];
    
    const struct gkyl_array *qmem = species->qmem;
    const struct gkyl_array *p_over_gamma = species->p_over_gamma;

    const struct gkyl_array *fin = species->f;
    struct gkyl_array *rhs = species->f1;

    // if (app->use_gpu)
    //   gkyl_hyper_dg_set_update_vol_cu(species->slvr, update_vol_term);
    // else
    //   gkyl_hyper_dg_set_update_vol(species->slvr, update_vol_term);
    gkyl_array_clear_range(rhs, 0.0, species->local);
    if (app->use_gpu)
      gkyl_dg_updater_vlasov_advance_cu(species->slvr, species->field_id, &species->local, 
        species->qmem, species->p_over_gamma, 
        fin, species->cflrate, rhs);
    else
      gkyl_dg_updater_vlasov_advance(species->slvr, species->field_id, &species->local, 
        species->qmem, species->p_over_gamma, 
        fin, species->cflrate, rhs);
  }
}

static void
range_stat_write(const char *nm, const struct gkyl_range *r, FILE *fp)
{
  fprintf(fp, " %s_cells : [ ", nm);
  for (int i=0; i<r->ndim; ++i)
    fprintf(fp, " %d, ", gkyl_range_shape(r, i));
  fprintf(fp, " ],\n");
}

void
gkyl_vlasov_app_stat_write(gkyl_vlasov_app* app)
{
  const char *fmt = "%s-%s";
  int sz = gkyl_calc_strlen(fmt, app->name, "stat.json");
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, "stat.json");
  
  char buff[70];
  time_t t = time(NULL);
  struct tm curr_tm = *localtime(&t);

  vm_species_coll_tm(app);
  vm_species_tm(app);

  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  with_file (fp, fileNm, "a") {
    fprintf(fp, "{\n");

    if (strftime(buff, sizeof buff, "%c", &curr_tm))
      fprintf(fp, " date : %s,\n", buff);

    fprintf(fp, " use_gpu : %d,\n", app->stat.use_gpu);

    for (int s=0; s<app->num_species; ++s)
      range_stat_write(app->species[s].info.name, &app->species[s].local, fp);
    
    fprintf(fp, " nup : %ld,\n", app->stat.nup);
    fprintf(fp, " nfeuler : %ld,\n", app->stat.nfeuler);
    fprintf(fp, " nstage_2_fail : %ld,\n", app->stat.nstage_2_fail);
    fprintf(fp, " nstage_3_fail : %ld,\n", app->stat.nstage_3_fail);

    fprintf(fp, " stage_2_dt_diff : [ %lg, %lg ],\n",
      app->stat.stage_2_dt_diff[0], app->stat.stage_2_dt_diff[1]);
    fprintf(fp, " stage_3_dt_diff : [ %lg, %lg ],\n",
      app->stat.stage_3_dt_diff[0], app->stat.stage_3_dt_diff[1]);
    
    fprintf(fp, " total_tm : %lg,\n", app->stat.total_tm);
    fprintf(fp, " init_species_tm : %lg,\n", app->stat.init_species_tm);
    if (app->has_field)
      fprintf(fp, " init_field_tm : %lg,\n", app->stat.init_field_tm);
    
    fprintf(fp, " species_rhs_tm : %lg,\n", app->stat.species_rhs_tm);

    for (int s=0; s<app->num_species; ++s) {
      fprintf(fp, " species_coll_drag_tm[%d] : %lg,\n", s,
        app->stat.species_lbo_coll_drag_tm[s]);
      fprintf(fp, " species_coll_diff_tm[%d] : %lg,\n", s,
        app->stat.species_lbo_coll_diff_tm[s]);
    }

    fprintf(fp, " species_coll_mom_tm : %lg,\n", app->stat.species_coll_mom_tm);
    fprintf(fp, " species_coll_tm : %lg,\n", app->stat.species_coll_tm);

    fprintf(fp, " fluid_species_rhs_tm : %lg,\n", app->stat.fluid_species_rhs_tm);
    
    if (app->has_field) {
      fprintf(fp, " field_rhs_tm : %lg,\n", app->stat.field_rhs_tm);
      fprintf(fp, " current_tm : %lg,\n", app->stat.current_tm);
    }

    fprintf(fp, " nmom : %ld,\n", app->stat.nmom);
    fprintf(fp, " mom_tm : %lg\n", app->stat.mom_tm);

    fprintf(fp, " ndiag : %ld,\n", app->stat.ndiag);
    fprintf(fp, " diag_tm : %lg\n", app->stat.diag_tm);

    fprintf(fp, " nspecies_omega_cfl : %ld,\n", app->stat.nspecies_omega_cfl);
    fprintf(fp, " species_omega_cfl_tm : %lg\n", app->stat.species_omega_cfl_tm);

    fprintf(fp, " nfield_omega_cfl : %ld,\n", app->stat.nfield_omega_cfl);
    fprintf(fp, " field_omega_cfl_tm : %lg\n", app->stat.field_omega_cfl_tm);
  
    fprintf(fp, "}\n");
  }
}

void
gkyl_vlasov_app_release(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i)
    vm_species_release(app, &app->species[i]);
  for (int i=0; i<app->num_fluid_species; ++i)
    vm_fluid_species_release(app, &app->fluid_species[i]);
  if (app->num_species > 0)
    gkyl_free(app->species);
  if (app->num_fluid_species > 0)
    gkyl_free(app->fluid_species);
  if (app->has_field)
    vm_field_release(app, app->field);

  if (app->use_gpu) {
    gkyl_cu_free(app->basis_on_dev.basis);
    gkyl_cu_free(app->basis_on_dev.confBasis);
  }

  gkyl_free(app);
}
