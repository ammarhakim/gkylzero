#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_priv.h>

// initialize fluid species object
void 
vm_fluid_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_fluid_species *f)
{
  int cdim = app->cdim;
  int vdim = app->vdim;
  // allocate fluid arrays
  f->fluid = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  f->fluid1 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  f->fluidnew = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  f->fluid_host = f->fluid;  
  if (app->use_gpu)
    f->fluid_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);

  // allocate buffer for applying BCs (used for both periodic and copy BCs)
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->cdim; ++d) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  f->bc_buffer = mkarr(app->use_gpu, app->confBasis.num_basis, buff_sz);

  // allocate cflrate (scalar array)
  f->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu)
    f->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
  else
    f->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  // allocate array to store advection velocity
  f->u = mkarr(app->use_gpu, cdim*app->confBasis.num_basis, app->local_ext.volume);
  
  // allocate array to store diffusion tensor
  f->D = mkarr(app->use_gpu, cdim*app->confBasis.num_basis, app->local_ext.volume);
  f->D_host = f->D;
  if (app->use_gpu)
    f->D_host = mkarr(false, cdim*app->confBasis.num_basis, app->local_ext.volume);
  
  // equation objects
  f->advect_eqn = gkyl_dg_advection_new(&app->confBasis, &app->local, app->use_gpu);
  f->diff_eqn = gkyl_dg_diffusion_new(&app->confBasis, &app->local, app->use_gpu);

  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_flags[GKYL_MAX_DIM] = {0, 0, 0};

  // fluid solvers
  f->advect_slvr = gkyl_hyper_dg_new(&app->grid, &app->confBasis, f->advect_eqn,
    app->cdim, up_dirs, zero_flux_flags, 1, app->use_gpu);
  f->diff_slvr = gkyl_hyper_dg_new(&app->grid, &app->confBasis, f->diff_eqn,
    app->cdim, up_dirs, zero_flux_flags, 1, app->use_gpu);

  f->has_advect = false;
  f->advects_with_species = false;
  // setup applied advection or advection with other species
  if (f->info.advection.velocity) {
    f->has_advect = true;
    // we need to ensure applied advection has same shape as current                                             
    f->advect = mkarr(app->use_gpu, cdim*app->confBasis.num_basis, app->local_ext.volume);

    f->advect_host = f->advect;
    if (app->use_gpu)
      f->advect_host = mkarr(false, cdim*app->confBasis.num_basis, app->local_ext.volume);

    f->advect_ctx = (struct vm_eval_advect_ctx) {
      .advect_func = f->info.advection.velocity, .advect_ctx = f->info.advection.velocity_ctx
    };

    f->advect_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid,
        .basis = &app->confBasis,
        .qtype = f->info.advection.qtype,
        .num_quad = app->confBasis.poly_order+1,
        .num_ret_vals = cdim,
        .eval = f->info.advection.velocity,
        .ctx = f->info.advection.velocity_ctx
      }
    );
  }
  else {
    f->advects_with_species = true;
    f->advection_species = vm_find_species(app, f->info.advection.advect_with);
    f->other_advect = f->advection_species->lbo.u_drift;   
    // determine collision type to use in vlasov update
    f->collision_id = f->info.advection.collision_id;
    if (f->collision_id == GKYL_LBO_COLLISIONS) {
      f->other_nu = f->advection_species->lbo.nu_sum;
      f->other_m0 = f->advection_species->lbo.m0;
      f->other_nu_vthsq = f->advection_species->lbo.nu_vthsq;     
      // allocate arrays to store collisional relaxation terms (nu*n*vthsq and nu*n*T_perp or nu*n*T_z)
      f->nu_fluid = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
      f->nu_n_vthsq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    }  
  }

  f->has_diffusion = false;
  if (f->info.diffusion.D) {
    f->has_diffusion = true;
    f->diff_ctx = (struct vm_eval_diffusion_ctx) {
      .diff_func = f->info.diffusion.D, .diff_ctx = f->info.diffusion.D_ctx
    };
    f->diff_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid,
        .basis = &app->confBasis,
        .qtype = GKYL_GAUSS_LOBATTO_QUAD,
        .num_quad = app->confBasis.poly_order+1,
        .num_ret_vals = cdim,
        .eval = f->info.diffusion.D,
        .ctx = f->info.diffusion.D_ctx
      }
    );
  }

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int dir=0; dir<app->cdim; ++dir) {
    if (is_np[dir]) {
      const enum gkyl_fluid_species_bc_type *bc;
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
  for (int d=0; d<3; ++d)
    f->absorb_bc_func[d] = gkyl_advection_absorb_bc_create(f->advect_eqn, d,
      app->basis_on_dev.confBasis);
}

void
vm_fluid_species_apply_ic(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double t0)
{
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    poly_order+1, 1, fluid_species->info.init, fluid_species->info.ctx);

  // run updater
  gkyl_proj_on_basis_advance(proj, t0, &app->local, fluid_species->fluid_host);
  gkyl_proj_on_basis_release(proj);

  if (app->use_gpu)
    gkyl_array_copy(fluid_species->fluid, fluid_species->fluid_host);

  // we are computing acceleration for now as it is time-independent
  vm_fluid_species_calc_advect(app, fluid_species, t0);
  // project diffusion tensor
  vm_fluid_species_calc_diff(app, fluid_species, t0);
}

void
vm_fluid_species_calc_advect(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm)
{
  if (fluid_species->has_advect) {
    gkyl_proj_on_basis_advance(fluid_species->advect_proj, tm, &app->local_ext, fluid_species->advect_host);
    if (app->use_gpu) // note: advect_host is same as advect when not on GPUs
      gkyl_array_copy(fluid_species->advect, fluid_species->advect_host);
  }
}

void
vm_fluid_species_calc_diff(gkyl_vlasov_app* app, struct vm_fluid_species* fluid_species, double tm)
{
  if (fluid_species->has_diffusion) {
    gkyl_proj_on_basis_advance(fluid_species->diff_proj, tm, &app->local_ext, fluid_species->D_host);
    if (app->use_gpu) // note: D_host is same as D when not on GPUs
      gkyl_array_copy(fluid_species->D, fluid_species->D_host);
  }
}

// Compute the RHS for fluid species update, returning maximum stable
// time-step.
double
vm_fluid_species_rhs(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  const struct gkyl_array *fluid, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  
  double omegaCfl = 1/DBL_MAX;
  
  gkyl_array_clear(fluid_species->cflrate, 0.0);
  gkyl_array_clear(fluid_species->u, 0.0);
  if (fluid_species->has_advect)
    gkyl_array_accumulate(fluid_species->u, 1.0, fluid_species->advect);

  if (fluid_species->advects_with_species) {
    // Need to apply boundary conditions to the drift velocity
    vm_fluid_species_apply_bc(app, fluid_species, fluid_species->other_advect);    
    gkyl_array_accumulate(fluid_species->u, 1.0, fluid_species->other_advect);
  }
  
  gkyl_advection_set_auxfields(fluid_species->advect_eqn,
    (struct gkyl_dg_advection_auxfields) { .u = fluid_species->u  }); // set total advection
  gkyl_diffusion_set_auxfields(fluid_species->diff_eqn,
    (struct gkyl_dg_diffusion_auxfields) { .D = fluid_species->D  }); // set total advection
  
  gkyl_array_clear(rhs, 0.0);

  if (app->use_gpu) {
    gkyl_hyper_dg_advance_cu(fluid_species->advect_slvr, &app->local, fluid, fluid_species->cflrate, rhs);
    if (fluid_species->has_diffusion)
      gkyl_hyper_dg_advance_cu(fluid_species->diff_slvr, &app->local, fluid, fluid_species->cflrate, rhs);
  } else {
    gkyl_hyper_dg_advance(fluid_species->advect_slvr, &app->local, fluid, fluid_species->cflrate, rhs);
    if (fluid_species->has_diffusion)
      gkyl_hyper_dg_advance(fluid_species->diff_slvr, &app->local, fluid, fluid_species->cflrate, rhs);
  }

  // accumulate nu*n*T - nu*fluid_species
  // where fluid_species = nT_perp or nT_z
  if (fluid_species->collision_id == GKYL_LBO_COLLISIONS) {
    gkyl_dg_mul_op_range(app->confBasis, 0, fluid_species->nu_fluid, 0,
      fluid_species->other_nu, 0, fluid, app->local);
    gkyl_dg_mul_op_range(app->confBasis, 0, fluid_species->nu_n_vthsq, 0,
      fluid_species->other_m0, 0, fluid_species->other_nu_vthsq, app->local);
    gkyl_array_accumulate(rhs, 1.0, fluid_species->nu_n_vthsq);
    gkyl_array_accumulate(rhs, -1.0, fluid_species->nu_fluid);
  }

  gkyl_array_reduce_range(fluid_species->omegaCfl_ptr, fluid_species->cflrate, GKYL_MAX, app->local);

  double omegaCfl_ho[1];
  if (app->use_gpu)
    gkyl_cu_memcpy(omegaCfl_ho, fluid_species->omegaCfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    omegaCfl_ho[0] = fluid_species->omegaCfl_ptr[0];
  omegaCfl = omegaCfl_ho[0];

  app->stat.fluid_species_rhs_tm += gkyl_time_diff_now_sec(wst);
  
  return app->cfl/omegaCfl;
}

// Apply periodic BCs on fluid species
void
vm_fluid_species_apply_periodic_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(fluid_species->bc_buffer->data, f, app->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, fluid_species->bc_buffer->data, app->skin_ghost.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(fluid_species->bc_buffer->data, f, app->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, fluid_species->bc_buffer->data, app->skin_ghost.lower_ghost[dir]);
}


// Apply copy BCs on fluid species
void
vm_fluid_species_apply_copy_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species,
  int dir, enum vm_domain_edge edge, struct gkyl_array *f)
{
  if (edge == VM_EDGE_LOWER) {
    gkyl_array_copy_to_buffer(fluid_species->bc_buffer->data, f, app->skin_ghost.lower_skin[dir]);
    gkyl_array_copy_from_buffer(f, fluid_species->bc_buffer->data, app->skin_ghost.lower_ghost[dir]);
  }

  if (edge == VM_EDGE_UPPER) {
    gkyl_array_copy_to_buffer(fluid_species->bc_buffer->data, f, app->skin_ghost.upper_skin[dir]);
    gkyl_array_copy_from_buffer(f, fluid_species->bc_buffer->data, app->skin_ghost.upper_ghost[dir]);
  }
}

// Apply absorbing BCs on fluid species
void
vm_fluid_species_apply_absorb_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species,
  int dir, enum vm_domain_edge edge, struct gkyl_array *f)
{
  
  if (edge == VM_EDGE_LOWER) {
    gkyl_array_copy_to_buffer_fn(fluid_species->bc_buffer->data, f, app->skin_ghost.lower_skin[dir],
      fluid_species->absorb_bc_func[dir]->on_dev
    );
    gkyl_array_copy_from_buffer(f, fluid_species->bc_buffer->data, app->skin_ghost.lower_ghost[dir]);
  }

  if (edge == VM_EDGE_UPPER) {
    gkyl_array_copy_to_buffer_fn(fluid_species->bc_buffer->data, f, app->skin_ghost.upper_skin[dir],
      fluid_species->absorb_bc_func[dir]->on_dev
    );
    gkyl_array_copy_from_buffer(f, fluid_species->bc_buffer->data, app->skin_ghost.upper_ghost[dir]);
  }  
}

// Determine which directions are periodic and which directions are copy,
// and then apply boundary conditions for fluid species
void
vm_fluid_species_apply_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d) {
    vm_fluid_species_apply_periodic_bc(app, fluid_species, app->periodic_dirs[d], f);
    is_np_bc[app->periodic_dirs[d]] = 0;
  }
  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {

      switch (fluid_species->lower_bc[d]) {
        case GKYL_FLUID_SPECIES_COPY:
          vm_fluid_species_apply_copy_bc(app, fluid_species, d, VM_EDGE_LOWER, f);
          break;
        case GKYL_FLUID_SPECIES_ABSORB:
          vm_fluid_species_apply_absorb_bc(app, fluid_species, d, VM_EDGE_LOWER, f);
          break;
      }

      switch (fluid_species->upper_bc[d]) {
        case GKYL_FLUID_SPECIES_COPY:
          vm_fluid_species_apply_copy_bc(app, fluid_species, d, VM_EDGE_UPPER, f);
          break;
        case GKYL_FLUID_SPECIES_ABSORB:
          vm_fluid_species_apply_absorb_bc(app, fluid_species, d, VM_EDGE_UPPER, f);
          break;
      }      
    }
  }
}

// release resources for fluid species
void
vm_fluid_species_release(const gkyl_vlasov_app* app, struct vm_fluid_species *f)
{
  gkyl_array_release(f->fluid);
  gkyl_array_release(f->fluid1);
  gkyl_array_release(f->fluidnew);
  gkyl_array_release(f->bc_buffer);
  gkyl_array_release(f->cflrate);

  gkyl_dg_eqn_release(f->advect_eqn);
  gkyl_hyper_dg_release(f->advect_slvr);
  gkyl_dg_eqn_release(f->diff_eqn);
  gkyl_hyper_dg_release(f->diff_slvr);

  gkyl_array_release(f->u);
  if (f->has_advect) {
    gkyl_array_release(f->advect);
    if (app->use_gpu)
      gkyl_array_release(f->advect_host);

    gkyl_proj_on_basis_release(f->advect_proj);
  }

  gkyl_array_release(f->D);
  if (app->use_gpu)
    gkyl_array_release(f->D_host);
  if (f->has_diffusion)
    gkyl_proj_on_basis_release(f->diff_proj);

  if (f->collision_id == GKYL_LBO_COLLISIONS) {
    gkyl_array_release(f->nu_fluid);
    gkyl_array_release(f->nu_n_vthsq);
  }

  if (app->use_gpu) {
    gkyl_array_release(f->fluid_host);
    gkyl_cu_free(f->omegaCfl_ptr);
  }
  else {
    gkyl_free(f->omegaCfl_ptr);
  }
  for (int d=0; d<3; ++d)
    gkyl_advection_bc_release(f->absorb_bc_func[d]);
}

