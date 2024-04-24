#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_priv.h>

#include <assert.h>
#include <float.h>
#include <time.h>

// function to evaluate external electromagnetic field (this is needed 
// as external electromagnetic field function provided by the user
// returns 6 components, while the Vlasov solver
// expects 8 components to match the EM field)
static void
eval_ext_em(double t, const double *xn, double *ext_em_out, void *ctx)
{
  struct vm_eval_ext_em_ctx *ext_em_ctx = ctx;
  double ext_em[6]; // output external EM field
  ext_em_ctx->ext_em_func(t, xn, ext_em, ext_em_ctx->ext_em_ctx);
  
  for (int i=0; i<6; ++i) ext_em_out[i] = ext_em[i];
  for (int i=6; i<8; ++i) ext_em_out[i] = 0.0;
}

// function to evaluate applied current (this is needed as 
// applied current function provided by the user returns 3 components,
// while the EM solver expects 8 components to match the EM field)
static void
eval_app_current(double t, const double *xn, double *app_current_out, void *ctx)
{
  struct vm_eval_app_current_ctx *app_current_ctx = ctx;
  double app_current[3]; // output applied current
  app_current_ctx->app_current_func(t, xn, app_current, app_current_ctx->app_current_ctx);
  
  for (int i=0; i<3; ++i) app_current_out[i] = app_current[i];
  for (int i=3; i<8; ++i) app_current_out[i] = 0.0;
}

// initialize field object
struct vm_field* 
vm_field_new(struct gkyl_vm *vm, struct gkyl_vlasov_app *app)
{
  struct vm_field *f = gkyl_malloc(sizeof(struct vm_field));

  f->info = vm->field;

  // allocate EM arrays
  f->em = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  f->em1 = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  f->emnew = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  f->em_energy = mkarr(app->use_gpu, 6, app->local_ext.volume);

  // allocate a total field variable for methods which require ext_em + em such as b_hat calculation
  f->tot_em = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);

  f->em_host = f->em;  
  if (app->use_gpu) {
    f->em_host = mkarr(false, 8*app->confBasis.num_basis, app->local_ext.volume);
    f->em_energy_red = gkyl_cu_malloc(sizeof(double[6]));
  }

  f->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 6);
  f->is_first_energy_write_call = true;

  f->has_ext_em = false;
  f->ext_em_evolve = false;
  // setup external electromagnetic field
  if (f->info.ext_em) {
    f->has_ext_em = true;
    if (f->info.ext_em_evolve)
      f->ext_em_evolve = f->info.ext_em_evolve;
    // we need to ensure external electromagnetic field has same shape as EM
    // field as it will get added to qmem
    f->ext_em = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);

    f->ext_em_host = f->ext_em;
    if (app->use_gpu)
      f->ext_em_host = mkarr(false, 8*app->confBasis.num_basis, app->local_ext.volume);

    f->ext_em_ctx = (struct vm_eval_ext_em_ctx) {
      .ext_em_func = f->info.ext_em, .ext_em_ctx = f->info.ext_em_ctx
    };
    f->ext_em_proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis, app->confBasis.poly_order+1,
      8, eval_ext_em, &f->ext_em_ctx);
  }

  f->has_app_current = false;
  f->app_current_evolve = false;
  // setup external currents
  if (f->info.app_current) {
    f->has_app_current = true;
    if (f->info.app_current_evolve)
      f->app_current_evolve = f->info.app_current_evolve;
    // we need to ensure external electromagnetic field has same shape as EM
    // field as it will get added to qmem
    f->app_current = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);

    f->app_current_host = f->app_current;
    if (app->use_gpu)
      f->app_current_host = mkarr(false, 8*app->confBasis.num_basis, app->local_ext.volume);

    f->app_current_ctx = (struct vm_eval_app_current_ctx) {
      .app_current_func = f->info.app_current, .app_current_ctx = f->info.app_current_ctx
    };
    f->app_current_proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis, app->confBasis.poly_order+1,
      8, eval_app_current, &f->app_current_ctx);
  }

  // allocate cflrate (scalar array)
  f->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu)
    f->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
  else
    f->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  // equation object
  double c = 1/sqrt(f->info.epsilon0*f->info.mu0);
  double ef = f->info.elcErrorSpeedFactor, mf = f->info.mgnErrorSpeedFactor;

  struct gkyl_dg_eqn *eqn;
  eqn = gkyl_dg_maxwell_new(&app->confBasis, c, ef, mf, app->use_gpu);

  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_flags[2*GKYL_MAX_DIM] = {0, 0, 0, 0, 0, 0};

  // Maxwell solver
  f->slvr = gkyl_hyper_dg_new(&app->grid, &app->confBasis, eqn,
    app->cdim, up_dirs, zero_flux_flags, 1, app->use_gpu);

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int dir=0; dir<app->cdim; ++dir) {
    f->lower_bc[dir] = f->upper_bc[dir] = GKYL_FIELD_COPY;
    if (is_np[dir]) {
      const enum gkyl_field_bc_type *bc;
      if (dir == 0)
        bc = f->info.bcx;
      else if (dir == 1)
        bc = f->info.bcy;
      else
        bc = f->info.bcz;

      f->lower_bc[dir] = bc[0];
      f->upper_bc[dir] = bc[1];
    }
  }

  // allocate buffer for applying BCs 
  long buff_sz = 0;
  // compute buffer size needed
  for (int dir=0; dir<app->cdim; ++dir) {
    long vol = GKYL_MAX2(app->lower_skin[dir].volume, app->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  f->bc_buffer = mkarr(app->use_gpu, 8*app->confBasis.num_basis, buff_sz);
  
  for (int d=0; d<app->cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (f->lower_bc[d] == GKYL_FIELD_COPY)
      bctype = GKYL_BC_COPY;
    else if (f->lower_bc[d] == GKYL_FIELD_PEC_WALL)
      bctype = GKYL_BC_MAXWELL_PEC;
    else if (f->lower_bc[d] == GKYL_FIELD_SYM_WALL)
      bctype = GKYL_BC_MAXWELL_SYM;
    else if (f->lower_bc[d] == GKYL_FIELD_RESERVOIR)
      bctype = GKYL_BC_MAXWELL_RESERVOIR;

    f->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev.confBasis,
      &app->lower_skin[d], &app->lower_ghost[d], f->em->ncomp, app->cdim, app->use_gpu);

    // Upper BC updater. Copy BCs by default.
    if (f->upper_bc[d] == GKYL_FIELD_COPY)
      bctype = GKYL_BC_COPY;
    else if (f->upper_bc[d] == GKYL_FIELD_PEC_WALL)
      bctype = GKYL_BC_MAXWELL_PEC;
    else if (f->upper_bc[d] == GKYL_FIELD_SYM_WALL)
      bctype = GKYL_BC_MAXWELL_SYM;
    else if (f->upper_bc[d] == GKYL_FIELD_RESERVOIR)
      bctype = GKYL_BC_MAXWELL_RESERVOIR;

    f->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.confBasis,
      &app->upper_skin[d], &app->upper_ghost[d], f->em->ncomp, app->cdim, app->use_gpu);
  }

  gkyl_dg_eqn_release(eqn);

  return f;
}

void
vm_field_apply_ic(gkyl_vlasov_app *app, struct vm_field *field, double t0)
{
  if (!app->has_field) return;
  
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    poly_order+1, 8, field->info.init, field->info.ctx);

  // run updater; need to project onto extended range for ease of handling
  // subsequent operations over extended range such as magnetic field unit vector computation
  // This is needed to fill the corner cells as the corner cells may not be filled by
  // boundary conditions and we cannot divide by 0 anywhere or the weak divisions will fail
  gkyl_proj_on_basis_advance(proj, t0, &app->local_ext, field->em_host);
  gkyl_proj_on_basis_release(proj);

  if (app->use_gpu)
    gkyl_array_copy(field->em, field->em_host);

  // pre-compute external EM field and applied current if present
  // pre-computation necessary in case external EM field or applied current
  // are time-independent and not computed in the time-stepping loop
  vm_field_calc_ext_em(app, field, t0);
  vm_field_calc_app_current(app, field, t0);
}

void
vm_field_calc_ext_em(gkyl_vlasov_app *app, struct vm_field *field, double tm)
{
  if (field->has_ext_em) {
    gkyl_proj_on_basis_advance(field->ext_em_proj, tm, &app->local_ext, field->ext_em_host);
    if (app->use_gpu) // note: ext_em_host is same as ext_em when not on GPUs
      gkyl_array_copy(field->ext_em, field->ext_em_host);
  }
}

void
vm_field_calc_app_current(gkyl_vlasov_app *app, struct vm_field *field, double tm)
{
  if (field->has_app_current) {
    gkyl_proj_on_basis_advance(field->app_current_proj, tm, &app->local_ext, field->app_current_host);
    if (app->use_gpu) // note: app_current_host is same as app_current when not on GPUs
      gkyl_array_copy(field->app_current, field->app_current_host);
  }
}

void
vm_field_accumulate_current(gkyl_vlasov_app *app, 
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[], 
  struct gkyl_array *emout)
{
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *s = &app->species[i];
    double qbyeps = s->info.charge/app->field->info.epsilon0; 

    vm_species_moment_calc(&s->m1i, s->local, app->local, fin[i]);
    gkyl_array_accumulate_range(emout, -qbyeps, s->m1i.marr, &app->local);
  } 
  // Accumulate applied current to electric field terms
  if (app->field->has_app_current)
    gkyl_array_accumulate_range(emout, -1.0/app->field->info.epsilon0, app->field->app_current, &app->local);
}

// Compute the RHS for field update, returning maximum stable
// time-step.
double
vm_field_rhs(gkyl_vlasov_app *app, struct vm_field *field,
  const struct gkyl_array *em, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  
  double omegaCfl = 1/DBL_MAX;
  
  gkyl_array_clear(field->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  if (!field->info.is_static) {
    gkyl_hyper_dg_advance(field->slvr, &app->local, em, field->cflrate, rhs);
    
    gkyl_array_reduce_range(field->omegaCfl_ptr, field->cflrate, GKYL_MAX, &app->local);

    app->stat.nfield_omega_cfl += 1;
    struct timespec tm = gkyl_wall_clock();
    
    double omegaCfl_ho[1];
    if (app->use_gpu)
      gkyl_cu_memcpy(omegaCfl_ho, field->omegaCfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
    else
      omegaCfl_ho[0] = field->omegaCfl_ptr[0];
    omegaCfl = omegaCfl_ho[0];

    app->stat.field_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  }

  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);
  
  return app->cfl/omegaCfl;
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for EM fields
void
vm_field_apply_bc(gkyl_vlasov_app *app, const struct vm_field *field, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();  
  
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext,
    num_periodic_dir, app->periodic_dirs, f);
  
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d)
    is_np_bc[app->periodic_dirs[d]] = 0;

  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {

      switch (field->lower_bc[d]) {
        case GKYL_FIELD_COPY:
        case GKYL_FIELD_PEC_WALL:
        case GKYL_FIELD_SYM_WALL:
        case GKYL_FIELD_RESERVOIR:
          gkyl_bc_basic_advance(field->bc_lo[d], field->bc_buffer, f);
          break;

        default:
          break;
      }

      switch (field->upper_bc[d]) {
        case GKYL_FIELD_COPY:
        case GKYL_FIELD_PEC_WALL:
        case GKYL_FIELD_SYM_WALL:
        case GKYL_FIELD_RESERVOIR:
          gkyl_bc_basic_advance(field->bc_up[d], field->bc_buffer, f);
          break;
          
        default:
          break;
      }   
    }
  }

  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, f);

  app->stat.field_bc_tm += gkyl_time_diff_now_sec(wst);
}

void
vm_field_calc_energy(gkyl_vlasov_app *app, double tm, const struct vm_field *field)
{
  for (int i=0; i<6; ++i)
    gkyl_dg_calc_l2_range(app->confBasis, i, field->em_energy, i, field->em, app->local);
  gkyl_array_scale_range(field->em_energy, app->grid.cellVolume, &app->local);
  
  double energy[6] = { 0.0 };
  if (app->use_gpu) {
    gkyl_array_reduce_range(field->em_energy_red, field->em_energy, GKYL_SUM, &app->local);
    gkyl_cu_memcpy(energy, field->em_energy_red, sizeof(double[6]), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    gkyl_array_reduce_range(energy, field->em_energy, GKYL_SUM, &app->local);
  }

  double energy_global[6] = { 0.0 };
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 6, energy, energy_global);
  
  gkyl_dynvec_append(field->integ_energy, tm, energy_global);
}

// release resources for field
void
vm_field_release(const gkyl_vlasov_app* app, struct vm_field *f)
{
  gkyl_array_release(f->em);
  gkyl_array_release(f->em1);
  gkyl_array_release(f->emnew);
  gkyl_array_release(f->tot_em);
  
  gkyl_array_release(f->bc_buffer);
  gkyl_array_release(f->cflrate);
  gkyl_array_release(f->em_energy);
  gkyl_dynvec_release(f->integ_energy);

  if (f->has_ext_em) {
    gkyl_array_release(f->ext_em);
    if (app->use_gpu)
      gkyl_array_release(f->ext_em_host);

    gkyl_proj_on_basis_release(f->ext_em_proj);
  }

  if (f->has_app_current) {
    gkyl_array_release(f->app_current);
    if (app->use_gpu)
      gkyl_array_release(f->app_current_host);

    gkyl_proj_on_basis_release(f->app_current_proj);
  }

  gkyl_hyper_dg_release(f->slvr);

  if (app->use_gpu) {
    gkyl_array_release(f->em_host);
    gkyl_cu_free(f->omegaCfl_ptr);
    gkyl_cu_free(f->em_energy_red);
  }
  else {
    gkyl_free(f->omegaCfl_ptr);
  }
  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    gkyl_bc_basic_release(f->bc_lo[d]);
    gkyl_bc_basic_release(f->bc_up[d]);
  }

  gkyl_free(f);
}

