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
  int num_eqn = f->info.num_eqn ? f->info.num_eqn : 1;
  // allocate fluid arrays
  f->fluid = mkarr(app->use_gpu, num_eqn*app->confBasis.num_basis, app->local_ext.volume);
  f->fluid1 = mkarr(app->use_gpu, num_eqn*app->confBasis.num_basis, app->local_ext.volume);
  f->fluidnew = mkarr(app->use_gpu, num_eqn*app->confBasis.num_basis, app->local_ext.volume);

  f->fluid_host = f->fluid;
  if (app->use_gpu)
    f->fluid_host = mkarr(false, num_eqn*app->confBasis.num_basis, app->local_ext.volume);

  // allocate buffer for applying BCs
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->cdim; ++d) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  f->bc_buffer = mkarr(app->use_gpu, num_eqn*app->confBasis.num_basis, buff_sz);

  // allocate cflrate (scalar array)
  f->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu)
    f->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
  else
    f->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  f->source_id = f->info.source.source_id;

  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_flags[GKYL_MAX_DIM] = {0, 0, 0};

  // pre-allocated memory for weak division
  f->u_mem = 0;
  if (app->use_gpu)
    f->u_mem = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
  else
    f->u_mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);

  // initialize pointers to flow velocity and pressure
  f->u = 0;
  f->u_host = 0;
  f->u_bc_buffer = 0;
  f->p = 0;
  f->p_host = 0;
  f->p_bc_buffer = 0;

  // initialize pointers to pkpm variables: 
  // div_b (divergence of magnetic field unit vector)
  // bb_grad_u (bb : grad(u))
  // div_p (divergence of the pressure tensor), 
  // p_force (total pressure forces in kinetic equation 1/rho (b . div(p) - p_perp div(b)))
  // u_perp_i (perpendicular bulk velocity u - u : bb)
  // rhou_perp_i (perpendicular momentum density rhou - rhou : bb)
  // p_perp (perpendicular pressure)
  f->div_b = 0;
  f->bb_grad_u = 0;
  f->div_p = 0;
  f->p_force = 0;
  f->u_perp_i = 0;
  f->rhou_perp_i = 0;
  f->p_perp = 0;

  f->param = 0.0;
  f->has_advect = false;
  f->advects_with_species = false;
  f->collision_id = GKYL_NO_COLLISIONS;
  // fluid solvers
  if (f->info.vt) {
    f->param = f->info.vt; // parameter for isothermal Euler is vt, thermal velocity
    f->eqn_id = GKYL_EQN_ISO_EULER;
    // allocate array to store fluid velocity
    f->u = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
    // allocate buffer for applying boundary conditions to primitive variables (u)
    f->u_bc_buffer = mkarr(app->use_gpu, 3*app->confBasis.num_basis, buff_sz);
  }
  else if (f->info.gas_gamma) {
    f->param = f->info.gas_gamma; // parameter for Euler is gas_gamma, adiabatic index
    f->eqn_id = GKYL_EQN_EULER;
    // allocate array to store fluid velocity and pressure
    f->u = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
    f->p = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    // allocate buffer for applying boundary conditions to primitive variables (u and p)
    f->u_bc_buffer = mkarr(app->use_gpu, 3*app->confBasis.num_basis, buff_sz);
    f->p_bc_buffer = mkarr(app->use_gpu, app->confBasis.num_basis, buff_sz);
  }
  else if (f->info.advection.velocity) {
    f->eqn_id = GKYL_EQN_ADVECTION;

    // allocate array to store advection velocity
    f->u = mkarr(app->use_gpu, cdim*app->confBasis.num_basis, app->local_ext.volume);
    f->u_bc_buffer = mkarr(app->use_gpu, cdim*app->confBasis.num_basis, buff_sz);
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
  }
  else {
    f->eqn_id = GKYL_EQN_EULER_PKPM;
    // allocate array to store fluid velocity and pressure tensor
    f->u = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
    f->p = mkarr(app->use_gpu, 6*app->confBasis.num_basis, app->local_ext.volume);
    // allocate buffer for applying boundary conditions to primitive variables (u and p)
    f->u_bc_buffer = mkarr(app->use_gpu, 3*app->confBasis.num_basis, buff_sz);
    f->p_bc_buffer = mkarr(app->use_gpu, 6*app->confBasis.num_basis, buff_sz);

    f->u_host = f->u;
    f->p_host = f->p;
    if (app->use_gpu) {
      f->u_host = mkarr(false, 3*app->confBasis.num_basis, app->local_ext.volume);
      f->p_host = mkarr(false, 6*app->confBasis.num_basis, app->local_ext.volume);
    }

    // allocate array to pkpm variables
    f->div_b = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    f->bb_grad_u = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    f->div_p = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
    f->p_force = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    f->u_perp_i = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
    f->rhou_perp_i = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
    f->p_perp = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    // Need additional boundary buffer for p_perp
    f->p_perp_bc_buffer = mkarr(app->use_gpu, app->confBasis.num_basis, buff_sz);

    f->pkpm_species = vm_find_species(app, f->info.pkpm_species);
    // index in fluid_species struct of fluid species kinetic species is colliding with
    f->species_index = vm_find_species_idx(app, f->info.pkpm_species);
  }

  f->advect_slvr = gkyl_dg_updater_fluid_new(&app->grid, &app->confBasis,
    &app->local, f->eqn_id, f->param, app->use_gpu);

  f->has_diffusion = false;
  f->diffusion_id = GKYL_NO_DIFFUSION;
  if (f->info.diffusion.D) {
    f->has_diffusion = true;
    // allocate array to store diffusion tensor
    int szD;
    if (f->info.diffusion.anisotropic) { // allocate space for mix terms {
      f->diffusion_id = GKYL_ANISO_DIFFUSION;
      szD = cdim*(cdim+1)/2;
    }
    else if (f->eqn_id == GKYL_EQN_ISO_EULER) {
      f->diffusion_id = GKYL_EULER_ISO_DIFFUSION;
      szD = 3;
    }
    else if (f->eqn_id == GKYL_EQN_EULER) {
      f->diffusion_id = GKYL_EULER_DIFFUSION;
      szD = 4;
    }
    else {
      f->diffusion_id = GKYL_ISO_DIFFUSION;
      szD = cdim;
    }

    f->D = mkarr(app->use_gpu, szD*app->confBasis.num_basis, app->local_ext.volume);
    f->D_host = f->D;
    if (app->use_gpu)
      f->D_host = mkarr(false, cdim*app->confBasis.num_basis, app->local_ext.volume);

    f->diff_ctx = (struct vm_eval_diffusion_ctx) {
      .diff_func = f->info.diffusion.D, .diff_ctx = f->info.diffusion.D_ctx
    };
    f->diff_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid,
        .basis = &app->confBasis,
        .qtype = GKYL_GAUSS_LOBATTO_QUAD,
        .num_quad = app->confBasis.poly_order+1,
        .num_ret_vals = szD,
        .eval = f->info.diffusion.D,
        .ctx = f->info.diffusion.D_ctx
      }
    );

    f->diff_slvr = gkyl_dg_updater_diffusion_new(&app->grid, &app->confBasis,
      &app->local, f->diffusion_id, app->use_gpu);
  }

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int dir=0; dir<app->cdim; ++dir) {
    f->lower_bc[dir] = f->upper_bc[dir] = GKYL_SPECIES_COPY;
    if (is_np[dir]) {
      const enum gkyl_species_bc_type *bc;
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

  int ghost[GKYL_MAX_DIM] = {0.0};
  for (int d=0; d<app->cdim; ++d)
    ghost[d] = 1;

  for (int d=0; d<app->cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (f->lower_bc[d] == GKYL_SPECIES_COPY)
      bctype = GKYL_BC_COPY;
    else if (f->lower_bc[d] == GKYL_SPECIES_ABSORB)
      bctype = GKYL_BC_ABSORB;

    f->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, &app->local_ext, ghost, bctype,
                                    app->basis_on_dev.confBasis, f->fluid->ncomp, app->cdim, app->use_gpu);
    f->bc_u_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, &app->local_ext, ghost, bctype,
                                    app->basis_on_dev.confBasis, f->u->ncomp, app->cdim, app->use_gpu);
    if (f->eqn_id == GKYL_EQN_EULER || f->eqn_id == GKYL_EQN_EULER_PKPM)
      f->bc_p_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, &app->local_ext, ghost, bctype,
                                      app->basis_on_dev.confBasis, f->p->ncomp, app->cdim, app->use_gpu);
    if (f->eqn_id == GKYL_EQN_EULER_PKPM)
      f->bc_p_perp_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, &app->local_ext, ghost, bctype,
                                      app->basis_on_dev.confBasis, f->p_perp->ncomp, app->cdim, app->use_gpu);

    // Upper BC updater. Copy BCs by default.
    if (f->upper_bc[d] == GKYL_SPECIES_COPY)
      bctype = GKYL_BC_COPY;
    else if (f->upper_bc[d] == GKYL_SPECIES_ABSORB)
      bctype = GKYL_BC_ABSORB;

    f->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, &app->local_ext, ghost, bctype,
                                    app->basis_on_dev.confBasis, f->fluid->ncomp, app->cdim, app->use_gpu);
    f->bc_u_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, &app->local_ext, ghost, bctype,
                                    app->basis_on_dev.confBasis, f->u->ncomp, app->cdim, app->use_gpu);
    if (f->eqn_id == GKYL_EQN_EULER || f->eqn_id == GKYL_EQN_EULER_PKPM)
      f->bc_p_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, &app->local_ext, ghost, bctype,
                                      app->basis_on_dev.confBasis, f->p->ncomp, app->cdim, app->use_gpu);
    if (f->eqn_id == GKYL_EQN_EULER_PKPM)
      f->bc_p_perp_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, &app->local_ext, ghost, bctype,
                                      app->basis_on_dev.confBasis, f->p_perp->ncomp, app->cdim, app->use_gpu);
  }
}

void
vm_fluid_species_apply_ic(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double t0)
{
  int poly_order = app->poly_order;
  int num_eqn = fluid_species->info.num_eqn ? fluid_species->info.num_eqn : 1;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    poly_order+1, num_eqn, fluid_species->info.init, fluid_species->info.ctx);

  // run updater
  gkyl_proj_on_basis_advance(proj, t0, &app->local, fluid_species->fluid_host);
  gkyl_proj_on_basis_release(proj);

  if (app->use_gpu)
    gkyl_array_copy(fluid_species->fluid, fluid_species->fluid_host);

  // we are computing acceleration for now as it is time-independent
  vm_fluid_species_calc_advect(app, fluid_species, t0);
  // project diffusion tensor
  vm_fluid_species_calc_diff(app, fluid_species, t0);

  // compute primitive variables at t = 0
  vm_fluid_species_prim_vars(app, fluid_species, fluid_species->fluid);

  // we are pre-computing source for now as it is time-independent
  vm_fluid_species_source_calc(app, fluid_species, t0);

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

void
vm_fluid_species_prim_vars(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  const struct gkyl_array *fluid)
{
  gkyl_array_clear(fluid_species->u, 0.0);

  // Compute bulk flow velocity, either from external advection or from state variables (rho*u = rhou)
  // Also compute pressure if present in equation system (Euler or Euler PKPM)
  if (fluid_species->eqn_id == GKYL_EQN_ISO_EULER) {
    gkyl_calc_prim_vars_u_from_statevec(fluid_species->u_mem, app->confBasis, &app->local,
      fluid, fluid_species->u);
  }
  else if (fluid_species->eqn_id == GKYL_EQN_EULER) {
    gkyl_calc_prim_vars_u_from_statevec(fluid_species->u_mem, app->confBasis, &app->local,
      fluid, fluid_species->u);
    // param is gas_gamma needed to compute pressure from energy
    gkyl_calc_prim_vars_p_from_statevec(app->confBasis, &app->local, fluid_species->param,
      fluid_species->u, fluid, fluid_species->p);
  }
  else if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM) {
    gkyl_calc_prim_vars_pkpm(app->confBasis, &app->local, app->field->bvar,
      fluid_species->pkpm_species->pkpm_moms.marr, fluid, 
      fluid_species->u, fluid_species->u_perp_i, fluid_species->rhou_perp_i,
      fluid_species->p_perp, fluid_species->p);
  }
  else {
    if (fluid_species->has_advect)
      gkyl_array_accumulate(fluid_species->u, 1.0, fluid_species->advect);
    if (fluid_species->advects_with_species)
      gkyl_array_accumulate(fluid_species->u, 1.0, fluid_species->other_advect);
  }
  // Note: Also applies boundary conditions to PKPM moments for use in fluid solve
  vm_fluid_species_prim_vars_apply_bc(app, fluid_species);

  if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM) {
    // calculate gradient quantities using recovery
    gkyl_array_clear(fluid_species->div_b, 0.0);
    gkyl_array_clear(fluid_species->bb_grad_u, 0.0);
    gkyl_array_clear(fluid_species->div_p, 0.0);
    gkyl_calc_prim_vars_pkpm_recovery(app->grid.dx, app->confBasis, &app->local, 
      app->field->bvar, fluid_species->u, fluid_species->p, 
      fluid_species->div_b, fluid_species->bb_grad_u, fluid_species->div_p);

    // calculate total pressure forces for kinetic equation coupling
    gkyl_array_clear(fluid_species->p_force, 0.0);
    gkyl_calc_prim_vars_pkpm_p_force(app->confBasis, &app->local, 
      app->field->bvar, fluid_species->div_p, fluid_species->pkpm_species->pkpm_moms.marr, 
      fluid_species->p_perp, fluid_species->div_b, 
      fluid_species->p_force);
  }
}

// Compute the RHS for fluid species update, returning maximum stable
// time-step.
double
vm_fluid_species_rhs(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  const struct gkyl_array *fluid, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  double omegaCfl = 1/DBL_MAX;

  gkyl_array_clear(fluid_species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  if (app->use_gpu) {
    if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM)
      gkyl_dg_updater_fluid_advance_cu(fluid_species->advect_slvr, fluid_species->eqn_id,
        &app->local, fluid_species->u, fluid_species->div_p, 
        fluid_species->u_perp_i, fluid_species->p_perp, 
        fluid, fluid_species->cflrate, rhs);
    else
      gkyl_dg_updater_fluid_advance_cu(fluid_species->advect_slvr, fluid_species->eqn_id,
        &app->local, fluid_species->u, fluid_species->p, 
        0, 0,  
        fluid, fluid_species->cflrate, rhs);

    if (fluid_species->has_diffusion)
      gkyl_dg_updater_diffusion_advance_cu(fluid_species->diff_slvr, fluid_species->diffusion_id,
        &app->local, fluid_species->D, fluid_species->u, fluid, fluid_species->cflrate, rhs);
  }
  else {
    if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM)
      gkyl_dg_updater_fluid_advance(fluid_species->advect_slvr, fluid_species->eqn_id,
        &app->local, fluid_species->u, fluid_species->div_p, 
        fluid_species->u_perp_i, fluid_species->p_perp,  
        fluid, fluid_species->cflrate, rhs);
    else 
      gkyl_dg_updater_fluid_advance(fluid_species->advect_slvr, fluid_species->eqn_id,
        &app->local, fluid_species->u, fluid_species->p, 
        0, 0, 
        fluid, fluid_species->cflrate, rhs);

    if (fluid_species->has_diffusion)
      gkyl_dg_updater_diffusion_advance(fluid_species->diff_slvr, fluid_species->diffusion_id,
        &app->local, fluid_species->D, fluid_species->u, fluid, fluid_species->cflrate, rhs);
  }

  // Relax T_perp, d T_perp/dt = nu*n*(T - T_perp)
  if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM) {
    double qbym = fluid_species->pkpm_species->info.charge/fluid_species->pkpm_species->info.mass;
    gkyl_array_set(fluid_species->pkpm_species->qmem, qbym, em);

    // Accumulate applied acceleration and/or q/m*(external electromagnetic)
    // fields onto qmem to get the total acceleration
    if (fluid_species->pkpm_species->has_accel)
      gkyl_array_accumulate(fluid_species->pkpm_species->qmem, 1.0, fluid_species->pkpm_species->accel);
    if (app->field->has_ext_em)
      gkyl_array_accumulate(fluid_species->pkpm_species->qmem, qbym, app->field->ext_em);

    gkyl_calc_prim_vars_pkpm_source(app->confBasis, &app->local, fluid_species->pkpm_species->qmem,
      fluid_species->pkpm_species->lbo.nu_sum, fluid_species->pkpm_species->lbo.nu_prim_moms,
      fluid_species->pkpm_species->pkpm_moms.marr, fluid,
      fluid_species->rhou_perp_i, fluid_species->p_perp,
      rhs);
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
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(fluid_species->bc_lo[d], fluid_species->bc_buffer, f);
          break;
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
        case GKYL_SPECIES_FIXED_FUNC:
          assert(false);
          break;
        default:
          break;
      }

      switch (fluid_species->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(fluid_species->bc_up[d], fluid_species->bc_buffer, f);
          break;
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
        case GKYL_SPECIES_FIXED_FUNC:
          assert(false);
          break;
        default:
          break;
      }
    }
  }
}

// Apply boundary conditions to primitive variables (bulk velocity u and pressure p)
void
vm_fluid_species_prim_vars_apply_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species)
{
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d) {
    gkyl_array_copy_to_buffer(fluid_species->u_bc_buffer->data, fluid_species->u, app->skin_ghost.lower_skin[d]);
    gkyl_array_copy_from_buffer(fluid_species->u, fluid_species->u_bc_buffer->data, app->skin_ghost.upper_ghost[d]);

    gkyl_array_copy_to_buffer(fluid_species->u_bc_buffer->data, fluid_species->u, app->skin_ghost.upper_skin[d]);
    gkyl_array_copy_from_buffer(fluid_species->u, fluid_species->u_bc_buffer->data, app->skin_ghost.lower_ghost[d]);
    // Apply boundary conditions to pressure if pressure present
    if (fluid_species->eqn_id == GKYL_EQN_EULER || fluid_species->eqn_id == GKYL_EQN_EULER_PKPM) {
      gkyl_array_copy_to_buffer(fluid_species->p_bc_buffer->data, fluid_species->p, app->skin_ghost.lower_skin[d]);
      gkyl_array_copy_from_buffer(fluid_species->p, fluid_species->p_bc_buffer->data, app->skin_ghost.upper_ghost[d]);

      gkyl_array_copy_to_buffer(fluid_species->p_bc_buffer->data, fluid_species->p, app->skin_ghost.upper_skin[d]);
      gkyl_array_copy_from_buffer(fluid_species->p, fluid_species->p_bc_buffer->data, app->skin_ghost.lower_ghost[d]);
    }
    if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM) {
      // Can reuse u_bc_buffer because u_bc_buffer is size 3*num_basis, same u_perp
      gkyl_array_copy_to_buffer(fluid_species->u_bc_buffer->data, fluid_species->u_perp_i, app->skin_ghost.lower_skin[d]);
      gkyl_array_copy_from_buffer(fluid_species->u_perp_i, fluid_species->u_bc_buffer->data, app->skin_ghost.upper_ghost[d]);

      gkyl_array_copy_to_buffer(fluid_species->u_bc_buffer->data, fluid_species->u_perp_i, app->skin_ghost.upper_skin[d]);
      gkyl_array_copy_from_buffer(fluid_species->u_perp_i, fluid_species->u_bc_buffer->data, app->skin_ghost.lower_ghost[d]);  

      gkyl_array_copy_to_buffer(fluid_species->p_perp_bc_buffer->data, fluid_species->p_perp, app->skin_ghost.lower_skin[d]);
      gkyl_array_copy_from_buffer(fluid_species->p_perp, fluid_species->p_perp_bc_buffer->data, app->skin_ghost.upper_ghost[d]);

      gkyl_array_copy_to_buffer(fluid_species->p_perp_bc_buffer->data, fluid_species->p_perp, app->skin_ghost.upper_skin[d]);
      gkyl_array_copy_from_buffer(fluid_species->p_perp, fluid_species->p_perp_bc_buffer->data, app->skin_ghost.lower_ghost[d]);  
    }
    is_np_bc[app->periodic_dirs[d]] = 0;
  }
  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {

      switch (fluid_species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(fluid_species->bc_u_lo[d], fluid_species->u_bc_buffer, fluid_species->u);
          if (fluid_species->eqn_id == GKYL_EQN_EULER || fluid_species->eqn_id == GKYL_EQN_EULER_PKPM)
            gkyl_bc_basic_advance(fluid_species->bc_p_lo[d], fluid_species->p_bc_buffer, fluid_species->p);
          if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM) {
            gkyl_bc_basic_advance(fluid_species->bc_u_lo[d], fluid_species->u_bc_buffer, fluid_species->u_perp_i);
            gkyl_bc_basic_advance(fluid_species->bc_p_perp_lo[d], fluid_species->p_perp_bc_buffer, fluid_species->p_perp);
          }
          break;
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }

      switch (fluid_species->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(fluid_species->bc_u_up[d], fluid_species->u_bc_buffer, fluid_species->u);
          if (fluid_species->eqn_id == GKYL_EQN_EULER || fluid_species->eqn_id == GKYL_EQN_EULER_PKPM)
            gkyl_bc_basic_advance(fluid_species->bc_p_up[d], fluid_species->p_bc_buffer, fluid_species->p);
          if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM) {
            gkyl_bc_basic_advance(fluid_species->bc_u_up[d], fluid_species->u_bc_buffer, fluid_species->u_perp_i);
            gkyl_bc_basic_advance(fluid_species->bc_p_perp_up[d], fluid_species->p_perp_bc_buffer, fluid_species->p_perp);
          }
          break;
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
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

  gkyl_dg_bin_op_mem_release(f->u_mem);

  gkyl_dg_updater_fluid_release(f->advect_slvr);
  if (f->has_diffusion) {
    gkyl_array_release(f->D);
    if (app->use_gpu)
      gkyl_array_release(f->D_host);

    gkyl_proj_on_basis_release(f->diff_proj);
    gkyl_dg_updater_diffusion_release(f->diff_slvr);
  }

  gkyl_array_release(f->u);
  gkyl_array_release(f->u_bc_buffer);
  if (f->eqn_id == GKYL_EQN_EULER || f->eqn_id == GKYL_EQN_EULER_PKPM) {
    gkyl_array_release(f->p);
    gkyl_array_release(f->p_bc_buffer);
  }
  if (f->eqn_id == GKYL_EQN_EULER_PKPM) {
    gkyl_array_release(f->div_b);
    gkyl_array_release(f->bb_grad_u);
    gkyl_array_release(f->div_p);
    gkyl_array_release(f->p_force);
    gkyl_array_release(f->u_perp_i);
    gkyl_array_release(f->rhou_perp_i);
    gkyl_array_release(f->p_perp);
  }
  if (f->has_advect) {
    gkyl_array_release(f->advect);
    if (app->use_gpu)
      gkyl_array_release(f->advect_host);

    gkyl_proj_on_basis_release(f->advect_proj);
  }

  if (f->source_id) {
    vm_fluid_species_source_release(app, &f->src);
  }

  if (app->use_gpu) {
    gkyl_array_release(f->fluid_host);
    gkyl_array_release(f->u_host);
    gkyl_cu_free(f->omegaCfl_ptr);
    if (f->eqn_id == GKYL_EQN_EULER || f->eqn_id == GKYL_EQN_EULER_PKPM) 
      gkyl_array_release(f->p_host);
  }
  else {
    gkyl_free(f->omegaCfl_ptr);
  }
  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    gkyl_bc_basic_release(f->bc_lo[d]);
    gkyl_bc_basic_release(f->bc_up[d]);
    gkyl_bc_basic_release(f->bc_u_lo[d]);
    gkyl_bc_basic_release(f->bc_u_up[d]);
    if (f->eqn_id == GKYL_EQN_EULER || f->eqn_id == GKYL_EQN_EULER_PKPM) {
      gkyl_bc_basic_release(f->bc_p_lo[d]);
      gkyl_bc_basic_release(f->bc_p_up[d]);
    }
    if (f->eqn_id == GKYL_EQN_EULER_PKPM) {
      gkyl_bc_basic_release(f->bc_p_perp_lo[d]);
      gkyl_bc_basic_release(f->bc_p_perp_up[d]);
    }
  }
}
