#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov_priv.h>

// initialize field object
struct vm_field* 
vm_field_new(struct gkyl_vm *vm, struct gkyl_vlasov_app *app)
{
  struct vm_field *f = gkyl_malloc(sizeof(struct vm_field));

  f->info = vm->field;

  if (vm->field.job_pool)
    // use specified job pool if it exists ...
    f->job_pool = gkyl_job_pool_acquire(vm->field.job_pool);
  else
    // .. or app job pool if it does not
    f->job_pool = gkyl_job_pool_acquire(app->job_pool);
  
  // allocate EM arrays
  f->em = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  f->em1 = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  f->emnew = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  f->qmem = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);

  f->em_host = f->em;  
  if (app->use_gpu)
    f->em_host = mkarr(false, 8*app->confBasis.num_basis, app->local_ext.volume);

  // allocate buffer for applying BCs (used for both periodic and copy BCs)
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->cdim; ++d) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  f->bc_buffer = mkarr(app->use_gpu, 8*app->confBasis.num_basis, buff_sz);

  // allocate cflrate (scalar array)
  f->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu)
    f->omegaCfl_ptr = gkyl_cu_malloc_host(sizeof(double));
  else
    f->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  // equation object
  double c = 1/sqrt(f->info.epsilon0*f->info.mu0);
  double ef = f->info.elcErrorSpeedFactor, mf = f->info.mgnErrorSpeedFactor;

  struct gkyl_dg_eqn *eqn;
  eqn = gkyl_dg_maxwell_new(&app->confBasis, c, ef, mf, app->use_gpu);

  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_flags[GKYL_MAX_DIM] = {0, 0, 0};

  // Maxwell solver
  f->slvr = gkyl_hyper_dg_new(&app->grid, &app->confBasis, eqn,
    app->cdim, up_dirs, zero_flux_flags, 1, app->use_gpu);

  gkyl_dg_eqn_release(eqn);

  return f;
}

void
vm_field_apply_ic(gkyl_vlasov_app *app, struct vm_field *field, double t0)
{
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    poly_order+1, 8, field->info.init, field->info.ctx);

  gkyl_proj_on_basis_advance(proj, t0, &app->local, field->em_host);
  gkyl_proj_on_basis_release(proj);

  if (app->use_gpu)
    gkyl_array_copy(field->em, field->em_host);
}

// Compute the RHS for field update, returning maximum stable
// time-step.
double
vm_field_rhs(gkyl_vlasov_app *app, struct vm_field *field,
  const struct gkyl_array *em, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  
  gkyl_array_clear(field->cflrate, 0.0);

  gkyl_array_clear(rhs, 0.0);
  if (app->use_gpu)
    gkyl_hyper_dg_advance_cu(field->slvr, &app->local, em, field->cflrate, rhs);
  else
    gkyl_hyper_dg_advance(field->slvr, &app->local, em, field->cflrate, rhs);

  gkyl_array_reduce_range(field->omegaCfl_ptr, field->cflrate, GKYL_MAX, app->local);
  double omegaCfl = field->omegaCfl_ptr[0];

  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);
  
  return app->cfl/omegaCfl;
}

// Apply periodic BCs on EM fields
void
vm_field_apply_periodic_bc(gkyl_vlasov_app *app, const struct vm_field *field,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(field->bc_buffer->data, f, app->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, field->bc_buffer->data, app->skin_ghost.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(field->bc_buffer->data, f, app->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, field->bc_buffer->data, app->skin_ghost.lower_ghost[dir]);
}

// Apply copy BCs on EM fields
void
vm_field_apply_copy_bc(gkyl_vlasov_app *app, const struct vm_field *field,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(field->bc_buffer->data, f, app->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, field->bc_buffer->data, app->skin_ghost.lower_ghost[dir]);

  gkyl_array_copy_to_buffer(field->bc_buffer->data, f, app->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, field->bc_buffer->data, app->skin_ghost.upper_ghost[dir]);
}

// Determine which directions are periodic and which directions are copy,
// and then apply boundary conditions for EM fields
void
vm_field_apply_bc(gkyl_vlasov_app *app, const struct vm_field *field, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim, is_copy[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    vm_field_apply_periodic_bc(app, field, app->periodic_dirs[d], f);
    is_copy[app->periodic_dirs[d]] = 0;
  }
  for (int d=0; d<cdim; ++d)
    if (is_copy[d])
      vm_field_apply_copy_bc(app, field, d, f);
}

// release resources for field
void
vm_field_release(const gkyl_vlasov_app* app, struct vm_field *f)
{
  gkyl_job_pool_release(f->job_pool);
  
  gkyl_array_release(f->em);
  gkyl_array_release(f->em1);
  gkyl_array_release(f->emnew);
  gkyl_array_release(f->qmem);
  gkyl_array_release(f->bc_buffer);
  gkyl_array_release(f->cflrate);

  gkyl_hyper_dg_release(f->slvr);

  if (app->use_gpu) {
    gkyl_array_release(f->em_host);
    gkyl_cu_free_host(f->omegaCfl_ptr);
  }
  else {
    gkyl_free(f->omegaCfl_ptr);
  }

  gkyl_free(f);
}

