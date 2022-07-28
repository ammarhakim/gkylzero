#include "gkyl_dg_bin_ops.h"
#include "gkyl_util.h"
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_vlasov_priv.h>
#include <gkyl_bc_basic.h>
#include <time.h>

// Projection functions for p/(m*gamma) = v in special relativistic systems
// Simplifies to p/sqrt(m^2 + p^2) where c = 1
static void 
ev_p_over_gamma_1p(double t, const double *xn, double *out, void *ctx)
{
  struct gamma_ctx *gtx = (struct gamma_ctx*) ctx;
  double mass = gtx->mass;
  out[0] = xn[0]/sqrt(1.0 + xn[0]*xn[0]);
}
static void 
ev_p_over_gamma_2p(double t, const double *xn, double *out, void *ctx)
{
  struct gamma_ctx *gtx = (struct gamma_ctx*) ctx;
  double mass = gtx->mass;
  out[0] = xn[0]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1]);
  out[1] = xn[1]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1]);
}
static void 
ev_p_over_gamma_3p(double t, const double *xn, double *out, void *ctx)
{
  struct gamma_ctx *gtx = (struct gamma_ctx*) ctx;
  double mass = gtx->mass;
  out[0] = xn[0]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
  out[1] = xn[1]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
  out[2] = xn[2]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
}

static const evalf_t p_over_gamma_func[3] = {ev_p_over_gamma_1p, ev_p_over_gamma_2p, ev_p_over_gamma_3p};

// function to evaluate acceleration (this is needed as accel function
// provided by the user returns 3 components, while the Vlasov solver
// expects 8 components to match the EM field)
static void
eval_accel(double t, const double *xn, double *aout, void *ctx)
{
  struct vm_eval_accel_ctx *a_ctx = ctx;
  double a[3]; // output acceleration
  a_ctx->accel_func(t, xn, a, a_ctx->accel_ctx);
  
  for (int i=0; i<3; ++i) aout[i] = a[i];
  for (int i=3; i<8; ++i) aout[i] = 0.0;
}

// initialize species object
void
vm_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_species *s)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = vm->cells[d];
    lower[d] = vm->lower[d];
    upper[d] = vm->upper[d];
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    // full phase space grid
    cells[cdim+d] = s->info.cells[d];
    lower[cdim+d] = s->info.lower[d];
    upper[cdim+d] = s->info.upper[d];
    ghost[cdim+d] = 0; // no ghost-cells in velocity space

    // only velocity space
    cells_vel[d] = s->info.cells[d];
    lower_vel[d] = s->info.lower[d];
    upper_vel[d] = s->info.upper[d];
    ghost_vel[d] = 0; // no ghost-cells in velocity space
  }
  // full phase space grid
  gkyl_rect_grid_init(&s->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->local_ext, &s->local);
  
  skin_ghost_ranges_init(&s->skin_ghost, &s->local_ext, ghost);
  
  // velocity space grid
  gkyl_rect_grid_init(&s->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&s->grid_vel, ghost_vel, &s->local_ext_vel, &s->local_vel);

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
    s->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
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
  
  // determine field-type 
  s->field_id = app->has_field ? app->field->info.field_id : GKYL_FIELD_NULL;

  // allocate array to store q/m*(E,B) or potential depending on equation system
  s->qmem = 0;
  s->fac_phi = 0;
  if (s->field_id  == GKYL_FIELD_E_B || s->field_id  == GKYL_FIELD_SR_E_B)
    s->qmem = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  else if (s->field_id == GKYL_FIELD_PHI || s->field_id == GKYL_FIELD_PHI_A)
    s->fac_phi = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  // allocate array to store vector potential if present
  s->vecA = 0;
  if (s->field_id == GKYL_FIELD_PHI_A)
    s->vecA = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);

  // allocate array to store p/gamma (velocity) if present
  // Since p/gamma is a geometric quantity, can pre-compute it here
  s->p_over_gamma = 0;
  if (s->field_id  == GKYL_FIELD_SR_E_B) {
    s->p_over_gamma = mkarr(app->use_gpu, vdim*app->velBasis.num_basis, s->local_vel.volume);
    s->p_over_gamma_host = s->p_over_gamma;
    if (app->use_gpu)
      s->p_over_gamma_host = mkarr(false, vdim*app->velBasis.num_basis, s->local_vel.volume);
    struct gamma_ctx ctx = { .mass = s->info.mass };
    gkyl_proj_on_basis *p_over_gamma_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &s->grid_vel,
        .basis = &app->velBasis,
        .qtype = GKYL_GAUSS_LOBATTO_QUAD,
        .num_quad = 8,
        .num_ret_vals = vdim,
        .eval = p_over_gamma_func[vdim-1],
        .ctx = &ctx
      }
    );  
    // run updater
    gkyl_proj_on_basis_advance(p_over_gamma_proj, 0.0, &s->local_vel, s->p_over_gamma_host);
    gkyl_proj_on_basis_release(p_over_gamma_proj);    
    if (app->use_gpu) // note: p_over_gamma_host is same as p_over_gamma when not on GPUs
      gkyl_array_copy(s->p_over_gamma, s->p_over_gamma_host);
  }

  // create solver
  s->slvr = gkyl_dg_updater_vlasov_new(&s->grid, &app->confBasis, &app->basis, 
    &app->local, &s->local_vel, s->field_id, app->use_gpu);

  // acquire equation object
  s->eqn_vlasov = gkyl_dg_updater_vlasov_acquire_eqn(s->slvr);
  
  // allocate data for momentum (for use in current accumulation)
  vm_species_moment_init(app, s, &s->m1i, "M1i");
  
  int ndm = s->info.num_diag_moments;
  // allocate data for diagnostic moments
  s->moms = gkyl_malloc(sizeof(struct vm_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    vm_species_moment_init(app, s, &s->moms[m], s->info.diag_moments[m]);

  // allocate data for integrated moments
  vm_species_moment_init(app, s, &s->integ_moms, "Integrated");

  if (app->use_gpu)
    s->red_integ_diag = gkyl_cu_malloc(sizeof(double[2+GKYL_MAX_DIM]));
  
  // allocate dynamic-vector to store all-reduced integrated moments
  s->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  s->is_first_integ_write_call = true;

  s->has_accel = false;
  // setup applied acceleration
  if (s->info.accel) {
    s->has_accel = true;
    // we need to ensure applied acceleration has same shape as EM
    // field as it will get added to qmem
    s->accel = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);

    s->accel_host = s->accel;
    if (app->use_gpu)
      s->accel_host = mkarr(false, 8*app->confBasis.num_basis, app->local_ext.volume);

    s->accel_ctx = (struct vm_eval_accel_ctx) {
      .accel_func = s->info.accel, .accel_ctx = s->info.accel_ctx
    };
    s->accel_proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis, app->confBasis.poly_order+1,
      8, eval_accel, &s->accel_ctx);
  }

  s->source_id = s->info.source.source_id;
  // setup constant source
  if (s->source_id) {
    vm_species_source_init(app, s, &s->src);
  }

  // determine collision type to use in vlasov update
  s->collision_id = s->info.collisions.collision_id;
  s->collides_with_fluid = false;
  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    if (vm_find_fluid_species(app, s->info.collisions.collide_with_fluid)) {
      s->collides_with_fluid = true;
      // index in fluid_species struct of fluid species kinetic species is colliding with
      s->fluid_index = vm_find_fluid_species_idx(app, s->info.collisions.collide_with_fluid);
    }
    vm_species_lbo_init(app, s, &s->lbo, s->collides_with_fluid);
  }

  // initialize boundary flux object
  if (s->calc_bflux)
    vm_species_bflux_init(app, s, &s->bflux);

  // setup mirror force from fluid species if present
  s->has_mirror_force = false;
  if (vm_find_fluid_species(app, s->info.mirror_force.fluid_mirror_force)) {
    s->has_mirror_force = true;
    // index in fluid_species struct of fluid species kinetic species 
    s->fluid_index = vm_find_fluid_species_idx(app, s->info.mirror_force.fluid_mirror_force);

    // gradB and Bmag are both scalar fields in configuration space
    s->magB = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    s->gradB = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    // Additional arrays for computing mirror force and current density for EM field coupling
    s->n = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    
    s->Tperp_mem = 0;
    if (app->use_gpu)
      s->Tperp_mem = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
    else
      s->Tperp_mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);
    
    s->Tperp = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    s->mirror_force = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    s->m1i_no_J = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

    s->gradB_host = s->gradB;
    s->magB_host = s->magB;
    s->n_host = s->n;
    s->Tperp_host = s->Tperp;
    s->mirror_force_host = s->mirror_force;
    s->m1i_no_J_host = s->m1i_no_J;
    if (app->use_gpu) {
      s->gradB_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
      s->magB_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
      s->n_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
      s->Tperp_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
      s->mirror_force_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
      s->m1i_no_J_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    }

    gkyl_proj_on_basis *gradB_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid,
        .basis = &app->confBasis,
        .qtype = GKYL_GAUSS_LOBATTO_QUAD,
        .num_quad = 8,
        .num_ret_vals = 1,
        .eval = s->info.mirror_force.gradB,
        .ctx = s->info.mirror_force.gradB_ctx
      }
    );

    // run updater
    gkyl_proj_on_basis_advance(gradB_proj, 0.0, &app->local, s->gradB_host);
    gkyl_proj_on_basis_release(gradB_proj);    
    if (app->use_gpu) // note: p_over_gamma_host is same as p_over_gamma when not on GPUs
      gkyl_array_copy(s->gradB, s->gradB_host);

    gkyl_proj_on_basis *magB_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid,
        .basis = &app->confBasis,
        .qtype = GKYL_GAUSS_LOBATTO_QUAD,
        .num_quad = 8,
        .num_ret_vals = 1,
        .eval = s->info.mirror_force.magB,
        .ctx = s->info.mirror_force.magB_ctx
      }
    );
    gkyl_proj_on_basis_advance(magB_proj, 0.0, &app->local, s->magB_host);
    gkyl_proj_on_basis_release(magB_proj);    
    if (app->use_gpu) // note: p_over_gamma_host is same as p_over_gamma when not on GPUs
      gkyl_array_copy(s->magB, s->magB_host);
  }

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int dir=0; dir<app->cdim; ++dir) {
    if (is_np[dir]) {
      const enum gkyl_species_bc_type *bc;
      if (dir == 0)
        bc = s->info.bcx;
      else if (dir == 1)
        bc = s->info.bcy;
      else
        bc = s->info.bcz;

      s->lower_bc[dir] = bc[0];
      s->upper_bc[dir] = bc[1];
    }
  }
  for (int d=0; d<app->cdim; ++d) {
    // Lower BC updater.
    enum gkyl_bc_basic_type bctype;
    if (s->lower_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype = GKYL_BC_REFLECT;
    } else if (s->lower_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype = GKYL_BC_ABSORB;
    }
    s->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, &s->local_ext, ghost, bctype,
                                    app->basis_on_dev.basis, app->cdim, app->use_gpu);
    // Upper BC updater.
    if (s->upper_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype = GKYL_BC_REFLECT;
    } else if (s->upper_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype = GKYL_BC_ABSORB;
    }
    s->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, &s->local_ext, ghost, bctype,
                                    app->basis_on_dev.basis, app->cdim, app->use_gpu);
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

  // we are pre-computing acceleration for now as it is time-independent
  vm_species_calc_accel(app, species, t0);

  // we are pre-computing source for now as it is time-independent
  vm_species_calc_source(app, species, t0);
}

void
vm_species_calc_accel(gkyl_vlasov_app *app, struct vm_species *species, double tm)
{
  if (species->has_accel) {
    gkyl_proj_on_basis_advance(species->accel_proj, tm, &app->local_ext, species->accel_host);
    if (app->use_gpu) // note: accel_host is same as accel when not on GPUs
      gkyl_array_copy(species->accel, species->accel_host);
  }
}

void
vm_species_calc_source(gkyl_vlasov_app *app, struct vm_species *species, double tm)
{
  if (species->source_id) {
    gkyl_proj_on_basis_advance(species->src.source_proj, tm, &species->local_ext, species->src.source_host);
    if (app->use_gpu) // note: source_host is same as source when not on GPUs
      gkyl_array_copy(species->src.source, species->src.source_host);
  }
}

// Compute the RHS for species update, returning maximum stable
// time-step.
double
vm_species_rhs(gkyl_vlasov_app *app, struct vm_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs,
  const struct gkyl_array *fluidin[])
{
  gkyl_array_clear(species->cflrate, 0.0);
  if (species->field_id  == GKYL_FIELD_E_B || species->field_id  == GKYL_FIELD_SR_E_B) {
    double qbym = species->info.charge/species->info.mass;
    gkyl_array_set(species->qmem, qbym, em);

    if (species->has_accel)
      gkyl_array_accumulate(species->qmem, 1.0, species->accel);

    if (species->has_mirror_force) {
      // get n = J*M0*magB (since J = 1/B)
      gkyl_dg_mul_op_range(app->confBasis, 0, species->n, 0,
        species->magB, 0, species->lbo.m0, app->local);
      // get J*T_perp = J*p_perp/n
      gkyl_dg_div_op_range(species->Tperp_mem, app->confBasis, 0, species->Tperp, 0,
        fluidin[species->fluid_index], 0, species->n, app->local);
      // get mirror force = J*T_perp*grad(B)
      gkyl_dg_mul_op_range(app->confBasis, 0, species->mirror_force, 0,
        species->gradB, 0, species->Tperp, app->local);

      gkyl_array_accumulate_range(species->qmem, -1.0, species->mirror_force, app->local);
    }
  }
  
  gkyl_array_clear(rhs, 0.0);
 
  if (app->use_gpu)
    gkyl_dg_updater_vlasov_advance_cu(species->slvr, species->field_id, &species->local, 
      species->qmem, species->p_over_gamma, 
      fin, species->cflrate, rhs);
  else
    gkyl_dg_updater_vlasov_advance(species->slvr, species->field_id, &species->local, 
      species->qmem, species->p_over_gamma, 
      fin, species->cflrate, rhs);

  if (species->collision_id == GKYL_LBO_COLLISIONS)
    vm_species_lbo_rhs(app, species, &species->lbo, fin, rhs);
  
  if (species->source_id)
    vm_species_source_rhs(app, species, &species->src, fin, rhs);

  // bflux calculation needs to be after source. The source uses bflux from the previous stage.
  // There's also an order of operations issue since the source may use bflux from a different
  // species, using the previous step insures all species bflux are updated before the source.
  if (species->calc_bflux)
    vm_species_bflux_rhs(app, species, &species->bflux, fin, rhs);
  
  app->stat.nspecies_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omegaCfl_ptr, species->cflrate, GKYL_MAX, species->local);

  double omegaCfl_ho[1];
  if (app->use_gpu)
    gkyl_cu_memcpy(omegaCfl_ho, species->omegaCfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    omegaCfl_ho[0] = species->omegaCfl_ptr[0];
  double omegaCfl = omegaCfl_ho[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  
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
  int dir, enum vm_domain_edge edge, struct gkyl_array *f)
{
  if (edge == VM_EDGE_LOWER) {
    gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.lower_skin[dir]);
    gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.lower_ghost[dir]);
  }

  if (edge == VM_EDGE_UPPER) {
    gkyl_array_copy_to_buffer(species->bc_buffer->data, f, species->skin_ghost.upper_skin[dir]);
    gkyl_array_copy_from_buffer(f, species->bc_buffer->data, species->skin_ghost.upper_ghost[dir]);
  }
}

// Determine which directions are periodic and which directions are copy,
// and then apply boundary conditions for distribution function
void
vm_species_apply_bc(gkyl_vlasov_app *app, const struct vm_species *species, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d) {
    vm_species_apply_periodic_bc(app, species, app->periodic_dirs[d], f);
    is_np_bc[app->periodic_dirs[d]] = 0;
  }
  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {

      switch (species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
          vm_species_apply_copy_bc(app, species, d, VM_EDGE_LOWER, f);
          break;
        case GKYL_SPECIES_REFLECT:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer, f);
          break;
        default:
          break;
      }

      switch (species->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
          vm_species_apply_copy_bc(app, species, d, VM_EDGE_UPPER, f);
          break;
        case GKYL_SPECIES_REFLECT:
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer, f);
          break;
        default:
          break;
      }      
    }
  }
}

void
vm_species_coll_tm(gkyl_vlasov_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      struct gkyl_dg_updater_lbo_vlasov_tm tm =
        gkyl_dg_updater_lbo_vlasov_get_tm(app->species[i].lbo.coll_slvr);
      app->stat.species_lbo_coll_diff_tm[i] = tm.diff_tm;
      app->stat.species_lbo_coll_drag_tm[i] = tm.drag_tm;
    }
  }
}

void
vm_species_tm(gkyl_vlasov_app *app)
{
  app->stat.species_rhs_tm = 0.0;
  for (int i=0; i<app->num_species; ++i) {
    struct gkyl_dg_updater_vlasov_tm tm =
      gkyl_dg_updater_vlasov_get_tm(app->species[i].slvr);
    app->stat.species_rhs_tm += tm.vlasov_tm;
  }
}

// release resources for species
void
vm_species_release(const gkyl_vlasov_app* app, const struct vm_species *s)
{
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
  vm_species_moment_release(app, &s->integ_moms);

  gkyl_dynvec_release(s->integ_diag);

  gkyl_dg_eqn_release(s->eqn_vlasov);
  gkyl_dg_updater_vlasov_release(s->slvr);

  // Release arrays for different types of Vlasov equations
  if (s->field_id  == GKYL_FIELD_E_B || s->field_id  == GKYL_FIELD_SR_E_B)
    gkyl_array_release(s->qmem);
  else if (s->field_id == GKYL_FIELD_PHI || s->field_id == GKYL_FIELD_PHI_A)
    gkyl_array_release(s->fac_phi);

  if (s->field_id == GKYL_FIELD_PHI_A)
    gkyl_array_release(s->vecA);

  if (s->field_id  == GKYL_FIELD_SR_E_B) {
    gkyl_array_release(s->p_over_gamma);
    if (app->use_gpu)
      gkyl_array_release(s->p_over_gamma_host);
  }
  
  if (s->has_accel) {
    gkyl_array_release(s->accel);
    if (app->use_gpu)
      gkyl_array_release(s->accel_host);

    gkyl_proj_on_basis_release(s->accel_proj);
  }

  if (s->source_id) {
    vm_species_source_release(app, &s->src);
  }

  if (s->has_mirror_force) {
    gkyl_array_release(s->gradB);
    gkyl_array_release(s->magB);
    gkyl_array_release(s->n);
    gkyl_dg_bin_op_mem_release(s->Tperp_mem);
    gkyl_array_release(s->Tperp);
    gkyl_array_release(s->mirror_force);
    gkyl_array_release(s->m1i_no_J);
    if (app->use_gpu) {
      gkyl_array_release(s->gradB_host);
      gkyl_array_release(s->magB_host);
      gkyl_array_release(s->n_host);
      gkyl_array_release(s->Tperp_host);
      gkyl_array_release(s->mirror_force_host);
      gkyl_array_release(s->m1i_no_J_host);
    }
  }

  if (s->collision_id == GKYL_LBO_COLLISIONS)
    vm_species_lbo_release(app, &s->lbo);

  if (s->calc_bflux)
    vm_species_bflux_release(app, &s->bflux);

  for (int d=0; d<app->cdim; ++d) {
    if (s->lower_bc[d]==GKYL_SPECIES_REFLECT || s->lower_bc[d]==GKYL_SPECIES_ABSORB)
      gkyl_bc_basic_release(s->bc_lo[d]);
    if (s->upper_bc[d]==GKYL_SPECIES_REFLECT || s->upper_bc[d]==GKYL_SPECIES_ABSORB)
      gkyl_bc_basic_release(s->bc_up[d]);
  }
  
  if (app->use_gpu) {
    gkyl_cu_free(s->omegaCfl_ptr);
    gkyl_cu_free(s->red_integ_diag);
  }
  else
    gkyl_free(s->omegaCfl_ptr);
}
