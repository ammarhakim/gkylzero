#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_priv.h>

// context for use in Wall BCs
struct maxwell_wall_bc_ctx {
  int dir; // direction for BCs
  const struct gkyl_basis *basis; // basis function
};

enum { M_EX, M_EY, M_EZ, M_BX, M_BY, M_BZ }; // components of EM field
static const int m_flip_even[3][3] = { // zero tangent E and zero normal B
  {M_BX, M_EY, M_EZ},
  {M_BY, M_EX, M_EZ},
  {M_BZ, M_EX, M_EY},
};
static const int m_flip_odd[3][3] = { // zero gradient
  { M_EX, M_BY, M_BZ },
  { M_EY, M_BX, M_BZ },
  { M_EZ, M_BX, M_BY },
};

static void
maxwell_wall_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct maxwell_wall_bc_ctx *mc = ctx;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  const int *feven = m_flip_even[dir];
  const int *fodd = m_flip_odd[dir];

  for (int i=0; i<3; ++i) {
    int eloc = nbasis*feven[i], oloc = nbasis*fodd[i];
    mc->basis->flip_even_sign(dir, &inp[eloc], &out[eloc]);
    mc->basis->flip_odd_sign(dir, &inp[oloc], &out[oloc]);
  }
  // correction potentials
  int eloc = nbasis*6, oloc = nbasis*7;
  mc->basis->flip_even_sign(dir, &inp[eloc], &out[eloc]);
  mc->basis->flip_odd_sign(dir, &inp[oloc], &out[oloc]);
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

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int dir=0; dir<app->cdim; ++dir) {
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
  
  double omegaCfl = 1/DBL_MAX;
  
  gkyl_array_clear(field->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  if (!field->info.is_static) {
    if (app->use_gpu)
      gkyl_hyper_dg_advance_cu(field->slvr, &app->local, em, field->cflrate, rhs);
    else
      gkyl_hyper_dg_advance(field->slvr, &app->local, em, field->cflrate, rhs);
    
    gkyl_array_reduce_range(field->omegaCfl_ptr, field->cflrate, GKYL_MAX, app->local);
    omegaCfl = field->omegaCfl_ptr[0];
  }

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
  int dir, enum vm_domain_edge edge, struct gkyl_array *f)
{
  if (edge == VM_EDGE_LOWER) {
    gkyl_array_copy_to_buffer(field->bc_buffer->data, f, app->skin_ghost.lower_skin[dir]);
    gkyl_array_copy_from_buffer(f, field->bc_buffer->data, app->skin_ghost.lower_ghost[dir]);
  }

  if (edge == VM_EDGE_UPPER) {
    gkyl_array_copy_to_buffer(field->bc_buffer->data, f, app->skin_ghost.upper_skin[dir]);
    gkyl_array_copy_from_buffer(f, field->bc_buffer->data, app->skin_ghost.upper_ghost[dir]);
  }
}

void
vm_field_apply_pec_bc(gkyl_vlasov_app *app, const struct vm_field *field,
  int dir, enum vm_domain_edge edge, struct gkyl_array *f)
{
  if (edge == VM_EDGE_LOWER) {
    gkyl_array_copy_to_buffer_fn(field->bc_buffer->data, f, app->skin_ghost.lower_skin[dir],
      maxwell_wall_bc, &(struct maxwell_wall_bc_ctx) { .dir = dir, .basis = &app->confBasis }
    );
    gkyl_array_copy_from_buffer(f, field->bc_buffer->data, app->skin_ghost.lower_ghost[dir]);
  }

  if (edge == VM_EDGE_UPPER) {
    gkyl_array_copy_to_buffer_fn(field->bc_buffer->data, f, app->skin_ghost.upper_skin[dir],
      maxwell_wall_bc, &(struct maxwell_wall_bc_ctx) { .dir = dir, .basis = &app->confBasis }
    );
    gkyl_array_copy_from_buffer(f, field->bc_buffer->data, app->skin_ghost.upper_ghost[dir]);
  }  
}

// Determine which directions are periodic and which directions are copy,
// and then apply boundary conditions for EM fields
void
vm_field_apply_bc(gkyl_vlasov_app *app, const struct vm_field *field, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d) {
    vm_field_apply_periodic_bc(app, field, app->periodic_dirs[d], f);
    is_np_bc[app->periodic_dirs[d]] = 0;
  }
  for (int d=0; d<cdim; ++d)
    if (is_np_bc[d]) {

      switch (field->lower_bc[d]) {
        case GKYL_FIELD_COPY:
          vm_field_apply_copy_bc(app, field, d, VM_EDGE_LOWER, f);
          break;
        case GKYL_FIELD_PEC_WALL:
          vm_field_apply_pec_bc(app, field, d, VM_EDGE_LOWER, f);
          break;
      }

      switch (field->upper_bc[d]) {
        case GKYL_FIELD_COPY:
          vm_field_apply_copy_bc(app, field, d, VM_EDGE_UPPER, f);
          break;
        case GKYL_FIELD_PEC_WALL:
          vm_field_apply_pec_bc(app, field, d, VM_EDGE_UPPER, f);
          break;
      }      
    }
}

// release resources for field
void
vm_field_release(const gkyl_vlasov_app* app, struct vm_field *f)
{
  gkyl_array_release(f->em);
  gkyl_array_release(f->em1);
  gkyl_array_release(f->emnew);
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

