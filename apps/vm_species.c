#include <assert.h>
#include <gkyl_vlasov_priv.h>

// initialize species object
void
vm_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_species *s)
{
  if (s->info.job_pool)
    // use specified job pool if it exists ...
    s->job_pool = gkyl_job_pool_acquire(s->info.job_pool);
  else
    // .. or app job pool if it does not
    s->job_pool = gkyl_job_pool_acquire(app->job_pool);
  
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = vm->cells[d];
    lower[d] = vm->lower[d];
    upper[d] = vm->upper[d];
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    cells[cdim+d] = s->info.cells[d];
    lower[cdim+d] = s->info.lower[d];
    upper[cdim+d] = s->info.upper[d];
    ghost[cdim+d] = 0; // no ghost-cells in velocity space
  }
  gkyl_rect_grid_init(&s->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->local_ext, &s->local);
  
  skin_ghost_ranges_init(&s->skin_ghost, &s->local_ext, ghost);
  
  // allocate distribution function arrays
  s->f = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  s->f1 = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  s->fnew = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);

  s->f_host = s->f;
  if (app->use_gpu)
    s->f_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);

  // allocate cflrate (scalar array)
  s->cflrate = mkarr(app->use_gpu, 1, s->local_ext.volume);
  if (app->use_gpu)
    s->omegaCfl_ptr = gkyl_cu_malloc_host(sizeof(double));
  else
    s->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  // allocate buffer for applying periodic BCs
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<cdim; ++d) {
    long vol = s->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  s->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  
  // allocate data for momentum (for use in current accumulation)
  vm_species_moment_init(app, s, &s->m1i, "M1i");
  
  int ndm = s->info.num_diag_moments;
  // allocate data for diagnostic moments
  s->moms = gkyl_malloc(sizeof(struct vm_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    vm_species_moment_init(app, s, &s->moms[m], s->info.diag_moments[m]);

  // determine field-type to use in dg_vlasov
  enum gkyl_field_id field_id = app->has_field ? app->field->info.field_id : GKYL_FIELD_NULL;

  // create equation object
  s->eqn = gkyl_dg_vlasov_new(&app->confBasis, &app->basis, &app->local, field_id, app->use_gpu);

  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    up_dirs[d] = d;
    zero_flux_flags[d] = 0;
  }
  int num_up_dirs = cdim;
  // update velocity space only when field is present
  if (app->has_field) {
    for (int d=cdim; d<pdim; ++d) {
      up_dirs[d] = d;
      zero_flux_flags[d] = 1; // zero-flux BCs in vel-space
    }
    num_up_dirs = pdim;
  }
  // create solver
  s->slvr = gkyl_hyper_dg_new(&s->grid, &app->basis, s->eqn,
    num_up_dirs, up_dirs, zero_flux_flags, 1, app->use_gpu);

  // determine collision type to use in vlasov update
  s->collision_id = s->info.collision_id;
  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    vm_species_lbo_init(app, s, &s->lbo);
  }
}

void
vm_species_apply_ic(gkyl_vlasov_app *app, struct vm_species *species, double t0)
{
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&species->grid, &app->basis,
    poly_order+1, 1, species->info.init, species->info.ctx);

  // run updater
  gkyl_proj_on_basis_advance(proj, t0, &species->local, species->f_host);
  gkyl_proj_on_basis_release(proj);    

  if (app->use_gpu) // note: f_host is same as f when not on GPUs
    gkyl_array_copy(species->f, species->f_host);
}

// Compute the RHS for species update, returning maximum stable
// time-step.
double
vm_species_rhs(gkyl_vlasov_app *app, struct vm_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *qmem, struct gkyl_array *rhs)
{
  gkyl_array_clear(species->cflrate, 0.0);
  if (qmem)
    gkyl_vlasov_set_auxfields(species->eqn, 
      (struct gkyl_dg_vlasov_auxfields) { .qmem = qmem }); // must set EM fields to use
  
  gkyl_array_clear(rhs, 0.0);

  struct timespec wst = gkyl_wall_clock();  
  if (app->use_gpu)
    gkyl_hyper_dg_advance_cu(species->slvr, &species->local, fin, species->cflrate, rhs);
  else
    gkyl_hyper_dg_advance(species->slvr, &species->local, fin, species->cflrate, rhs);
  app->stat.species_rhs_tm += gkyl_time_diff_now_sec(wst);

  wst = gkyl_wall_clock();
  if (species->collision_id == GKYL_LBO_COLLISIONS)
    vm_species_lbo_rhs(app, species, &species->lbo, fin, rhs);
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);      

  gkyl_array_reduce_range(species->omegaCfl_ptr, species->cflrate, GKYL_MAX, species->local);
  double omegaCfl = species->omegaCfl_ptr[0];
  
  return app->cfl/omegaCfl;
}

// Apply periodic BCs on distribution function
void
vm_species_apply_periodic_bc(gkyl_vlasov_app *app, const struct vm_species *species,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.lower_ghost[dir]);
}

// Apply copy BCs on distribution function
void
vm_species_apply_copy_bc(gkyl_vlasov_app *app, const struct vm_species *species,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.lower_ghost[dir]);

  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.upper_ghost[dir]);
}

// Determine which directions are periodic and which directions are copy,
// and then apply boundary conditions for distribution function
void
vm_species_apply_bc(gkyl_vlasov_app *app, const struct vm_species *species, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim, is_copy[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    vm_species_apply_periodic_bc(app, species, app->periodic_dirs[d], f);
    is_copy[app->periodic_dirs[d]] = 0;
  }
  for (int d=0; d<cdim; ++d)
    if (is_copy[d])
      vm_species_apply_copy_bc(app, species, d, f);
}

// release resources for species
void
vm_species_release(const gkyl_vlasov_app* app, const struct vm_species *s)
{
  gkyl_job_pool_release(s->job_pool);
  
  // release various arrays
  gkyl_array_release(s->f);
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->cflrate);
  gkyl_array_release(s->bc_buffer);
  
  if (app->use_gpu)
    gkyl_array_release(s->f_host);
  
  // release moment data
  vm_species_moment_release(app, &s->m1i);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    vm_species_moment_release(app, &s->moms[i]);
  gkyl_free(s->moms);

  gkyl_dg_eqn_release(s->eqn);
  gkyl_hyper_dg_release(s->slvr);  

  if (s->collision_id == GKYL_LBO_COLLISIONS)
    vm_species_lbo_release(app, &s->lbo);
  
  if (app->use_gpu)
    gkyl_cu_free_host(s->omegaCfl_ptr);
  else
    gkyl_free(s->omegaCfl_ptr);
}
