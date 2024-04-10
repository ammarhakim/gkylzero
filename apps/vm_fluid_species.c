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
  f->eqn_type = f->info.equation->type;
  f->num_equations = f->info.equation->num_equations;
  f->equation = gkyl_wv_eqn_acquire(f->info.equation);

  // allocate fluid arrays
  f->fluid = mkarr(app->use_gpu, f->num_equations*app->confBasis.num_basis, app->local_ext.volume);
  f->fluid1 = mkarr(app->use_gpu, f->num_equations*app->confBasis.num_basis, app->local_ext.volume);
  f->fluidnew = mkarr(app->use_gpu, f->num_equations*app->confBasis.num_basis, app->local_ext.volume);

  f->fluid_host = f->fluid;
  if (app->use_gpu)
    f->fluid_host = mkarr(false, f->num_equations*app->confBasis.num_basis, app->local_ext.volume);

  // Duplicate copy of fluid data in case time step fails.
  // Needed because of implicit source split which modifies solution and 
  // is always successful, so if a time step fails due to the SSP RK3 
  // we must restore the old solution before restarting the time step
  f->fluid_dup = mkarr(app->use_gpu, f->num_equations*app->confBasis.num_basis, app->local_ext.volume);

  // allocate cflrate (scalar array)
  f->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu)
    f->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
  else
    f->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_flags[GKYL_MAX_DIM] = {0, 0, 0};

  // fluid solvers
  if (f->info.advection.velocity) {
    f->eqn_type = GKYL_EQN_ADVECTION;

    // setup applied advection 
    f->app_advect = mkarr(app->use_gpu, cdim*app->confBasis.num_basis, app->local_ext.volume);

    f->app_advect_host = f->app_advect;
    if (app->use_gpu)
      f->app_advect_host = mkarr(false, cdim*app->confBasis.num_basis, app->local_ext.volume);

    gkyl_proj_on_basis *app_advect_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid,
        .basis = &app->confBasis,
        .qtype = GKYL_GAUSS_LOBATTO_QUAD,
        .num_quad = 8,
        .num_ret_vals = cdim,
        .eval = f->info.advection.velocity,
        .ctx = f->info.advection.velocity_ctx
      }
    );

    gkyl_proj_on_basis_advance(app_advect_proj, 0.0, &app->local_ext, f->app_advect_host);
    if (app->use_gpu) // note: app_advect_host is same as advect when not on GPUs
      gkyl_array_copy(f->app_advect, f->app_advect_host);
    // Free projection object
    gkyl_proj_on_basis_release(app_advect_proj);

    struct gkyl_dg_advection_auxfields aux_inp = {.u_i = f->app_advect};
    f->advect_slvr = gkyl_dg_updater_fluid_new(&app->grid, &app->confBasis,
      &app->local, f->equation, app->geom, &aux_inp, app->use_gpu);
  }
  else {
    // allocate array to store fluid velocity (ux, uy, uz) and pressure
    // For isothermal Euler, p : (vth*rho)
    // For Euler, p : (gamma - 1)*(E - 1/2 rho u^2)
    f->u = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
    f->p = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    // boolean array for if we are only using the cell average for primitive variables
    f->cell_avg_prim = mk_int_arr(app->use_gpu, 1, app->local_ext.volume);

    // Allocate arrays for kinetic energy at old and new time steps.
    // These are used because implicit source solve updates momentum but does not affect 
    // the pressure, so we can construct the updated energy from the updated momentum.
    f->ke_old = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);    
    f->ke_new = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);    

    int Nbasis_surf = app->confBasis.num_basis/(app->confBasis.poly_order + 1); // *only valid for tensor bases for cdim > 1*
    // Surface primitive variables (2*cdim*3 components). Ordered as:
    // [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
    //  ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 
    //  ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr] 
    f->u_surf = mkarr(app->use_gpu, 2*cdim*3*Nbasis_surf, app->local_ext.volume);
    // Surface pressure tensor (2*cdim components). Ordered as:
    // [p_xl, p_xr, p_yl, p_yr, p_zl, p_zr]
    f->p_surf = mkarr(app->use_gpu, 2*cdim*Nbasis_surf, app->local_ext.volume);

    // Check if limiter_fac is specified for adjusting how much diffusion is applied through slope limiter
    // If not specified, set to 0.0 and updater sets default behavior (1/sqrt(3); see gkyl_dg_calc_fluid_vars.h)
    double limiter_fac = f->info.limiter_fac == 0 ? 0.0 : f->info.limiter_fac;

    // updater for computing fluid variables: flow velocity and pressure
    // also stores kernels for computing source terms, integrated variables
    // Two instances, one over extended range and one over local range for ease of handling boundary conditions
    f->calc_fluid_vars_ext = gkyl_dg_calc_fluid_vars_new(f->equation, &app->confBasis, &app->local_ext, limiter_fac, app->use_gpu);
    f->calc_fluid_vars = gkyl_dg_calc_fluid_vars_new(f->equation, &app->confBasis, &app->local, limiter_fac, app->use_gpu); 

    struct gkyl_dg_euler_auxfields aux_inp = {.u = f->u, .p = f->p, 
      .u_surf = f->u_surf, .p_surf = f->p_surf};
    f->advect_slvr = gkyl_dg_updater_fluid_new(&app->grid, &app->confBasis,
      &app->local, f->equation, app->geom, &aux_inp, app->use_gpu);    
  }

  f->has_diffusion = false;
  f->diffD = NULL;
  if (f->info.diffusion.Dij) {
    f->has_diffusion = true;
    // allocate space for full diffusion tensor 
    int szD = cdim*(cdim+1)/2;

    f->diffD = mkarr(app->use_gpu, szD*app->confBasis.num_basis, app->local_ext.volume);
    struct gkyl_array *diffD_host = f->diffD;
    if (app->use_gpu)
      diffD_host = mkarr(false, szD*app->confBasis.num_basis, app->local_ext.volume);

    gkyl_proj_on_basis *diff_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
        .grid = &app->grid,
        .basis = &app->confBasis,
        .qtype = GKYL_GAUSS_LOBATTO_QUAD,
        .num_quad = 8,
        .num_ret_vals = szD,
        .eval = f->info.diffusion.Dij,
        .ctx = f->info.diffusion.Dij_ctx
      }
    );
    gkyl_proj_on_basis_advance(diff_proj, 0.0, &app->local_ext, diffD_host);
    if (app->use_gpu) {// note: diffD_host is same as diffD when not on GPUs
      gkyl_array_copy(f->diffD, diffD_host);
      gkyl_array_release(diffD_host);
    }
    // Free projection object
    gkyl_proj_on_basis_release(diff_proj);

    f->diff_slvr_gen = gkyl_dg_updater_diffusion_gen_new(&app->grid, &app->confBasis,
      &app->local, app->use_gpu);
  }
  else if (f->info.diffusion.D) {
    f->has_diffusion = true;
    f->info.diffusion.order = f->info.diffusion.order<2? 2 : f->info.diffusion.order;

    int szD = cdim;
    f->diffD = mkarr(app->use_gpu, szD, 1);
    struct gkyl_array *diffD_host = f->diffD;
    if (app->use_gpu)
      diffD_host = mkarr(false, szD, 1);
    // Set diffusion coefficient in each direction to input value.
    gkyl_array_clear(diffD_host, 0.);
    for (int d=0; d<cdim; d++) gkyl_array_shiftc(diffD_host, f->info.diffusion.D, d);

    if (app->use_gpu) {// note: diffD_host is same as diffD when not on GPUs
      gkyl_array_copy(f->diffD, diffD_host);
      gkyl_array_release(diffD_host);
    }

    const bool is_zero_flux[GKYL_MAX_CDIM] = {false};

    f->diff_slvr = gkyl_dg_updater_diffusion_fluid_new(&app->grid, &app->confBasis,
      true, f->num_equations, NULL, f->info.diffusion.order, &app->local, is_zero_flux, app->use_gpu);
  }

  // array for storing integrated moments in each cell
  f->integ_mom = mkarr(app->use_gpu, 6, app->local_ext.volume);
  if (app->use_gpu) {
    f->red_integ_diag = gkyl_cu_malloc(sizeof(double[6]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments 
  f->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, 6);
  f->is_first_integ_write_call = true;

  // Initialize applied acceleration for use in force update. 
  // Always used by fluid implicit sources, so always initialize.
  f->app_accel = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(f->app_accel, 0.0);
  f->has_app_accel = false;
  f->app_accel_evolve = false;
  // setup applied acceleration
  if (f->info.app_accel) {
    f->has_app_accel = true;
    if (f->info.app_accel_evolve) {
      f->app_accel_evolve = f->info.app_accel_evolve;
    }

    f->app_accel_host = f->app_accel;
    if (app->use_gpu) {
      f->app_accel_host = mkarr(false, 3*app->confBasis.num_basis, app->local_ext.volume);
    }
    f->app_accel_proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis, app->confBasis.poly_order+1,
      3, f->info.app_accel, f->info.app_accel_ctx);
  }

  // set species source id
  f->source_id = f->info.source.source_id;

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

  // allocate buffer for applying BCs 
  long buff_sz = 0;
  // compute buffer size needed
  for (int dir=0; dir<app->cdim; ++dir) {
    long vol = GKYL_MAX2(app->lower_skin[dir].volume, app->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  f->bc_buffer = mkarr(app->use_gpu, f->num_equations*app->confBasis.num_basis, buff_sz);

  // Certain operations fail if absorbing BCs used because absorbing BCs 
  // means the mass density is 0 in the ghost cells (divide by zero)
  f->bc_is_absorb = false;
  for (int d=0; d<app->cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (f->lower_bc[d] == GKYL_SPECIES_COPY) {
      bctype = GKYL_BC_COPY;
    }
    else if (f->lower_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype = GKYL_BC_ABSORB;
      f->bc_is_absorb = true;
    }
    else if (f->lower_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype = GKYL_BC_EULER_REFLECT;
    }
    else if (f->lower_bc[d] == GKYL_SPECIES_NO_SLIP) {
      bctype = GKYL_BC_EULER_NO_SLIP;
    }

    f->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev.confBasis,
      &app->lower_skin[d], &app->lower_ghost[d], f->fluid->ncomp, app->cdim, app->use_gpu);

    // Upper BC updater. Copy BCs by default.
    if (f->upper_bc[d] == GKYL_SPECIES_COPY) {
      bctype = GKYL_BC_COPY;
    }
    else if (f->upper_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype = GKYL_BC_ABSORB;
      f->bc_is_absorb = true;
    }
    else if (f->upper_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype = GKYL_BC_EULER_REFLECT;
    }
    else if (f->upper_bc[d] == GKYL_SPECIES_NO_SLIP) {
      bctype = GKYL_BC_EULER_NO_SLIP;
    }

    f->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.confBasis,
      &app->upper_skin[d], &app->upper_ghost[d], f->fluid->ncomp, app->cdim, app->use_gpu);
  }
}

void
vm_fluid_species_apply_ic(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double t0)
{
  int poly_order = app->poly_order;

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    poly_order+1, fluid_species->num_equations, fluid_species->info.init, fluid_species->info.ctx);

  // run updater
  gkyl_proj_on_basis_advance(proj, t0, &app->local_ext, fluid_species->fluid_host);
  gkyl_proj_on_basis_release(proj);

  if (app->use_gpu) {
    gkyl_array_copy(fluid_species->fluid, fluid_species->fluid_host);
  }
  // Apply limiter at t=0 to insure slopes are well-behaved at beginning of simulation
  vm_fluid_species_limiter(app, fluid_species, fluid_species->fluid);

  // Pre-compute applied acceleration in case it's time-independent
  vm_fluid_species_calc_app_accel(app, fluid_species, t0);

  // we are pre-computing source for now as it is time-independent
  vm_fluid_species_source_calc(app, fluid_species, t0);
}

void
vm_fluid_species_calc_app_accel(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm)
{
  if (fluid_species->has_app_accel) {
    gkyl_proj_on_basis_advance(fluid_species->app_accel_proj, tm, &app->local_ext, fluid_species->app_accel_host);
    // note: app_accel_host is same as app_accel when not on GPUs
    if (app->use_gpu) {
      gkyl_array_copy(fluid_species->app_accel, fluid_species->app_accel_host);
    }
  }
}

void
vm_fluid_species_prim_vars(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  const struct gkyl_array *fluid)
{
  struct timespec tm = gkyl_wall_clock();

  // Compute flow velocity in both the volume and on surfaces
  if (fluid_species->bc_is_absorb) {
    gkyl_dg_calc_fluid_vars_advance(fluid_species->calc_fluid_vars,
      fluid, fluid_species->cell_avg_prim, fluid_species->u, fluid_species->u_surf); 
  }
  else {
    gkyl_dg_calc_fluid_vars_advance(fluid_species->calc_fluid_vars_ext,
      fluid, fluid_species->cell_avg_prim, fluid_species->u, fluid_species->u_surf); 
  }

  // Compute scalar pressure in the volume and at needed surfaces
  gkyl_dg_calc_fluid_vars_pressure(fluid_species->calc_fluid_vars, 
    &app->local_ext, fluid, fluid_species->u, fluid_species->p, fluid_species->p_surf);

  app->stat.fluid_species_vars_tm += gkyl_time_diff_now_sec(tm);
}

void
vm_fluid_species_limiter(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  struct gkyl_array *fluid)
{
  struct timespec tm = gkyl_wall_clock();

  // Limit the slopes of the solution
  gkyl_dg_calc_fluid_vars_limiter(fluid_species->calc_fluid_vars, &app->local, fluid);

  app->stat.fluid_species_vars_tm += gkyl_time_diff_now_sec(tm);

  // Apply boundary conditions after limiting solution
  vm_fluid_species_apply_bc(app, fluid_species, fluid);
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

  gkyl_dg_updater_fluid_advance(fluid_species->advect_slvr, 
    &app->local, fluid, fluid_species->cflrate, rhs);

  // Accumulate explicit source contribution, e.g., external forces
  // Only done if there are external forces and no EM fields, as fluid-EM coupling
  // is handled by implicit source solve, see vm_fluid_em_coupling.c
  if (fluid_species->has_app_accel && !app->has_field) {
    gkyl_dg_calc_fluid_vars_source(fluid_species->calc_fluid_vars, &app->local, 
      fluid_species->app_accel, fluid, rhs); 
  }

  if (fluid_species->has_diffusion) {
    if (fluid_species->info.diffusion.Dij) {
      gkyl_dg_updater_diffusion_gen_advance(fluid_species->diff_slvr_gen,
        &app->local, fluid_species->diffD, fluid, fluid_species->cflrate, rhs);
    }
    else if (fluid_species->info.diffusion.D) {
      gkyl_dg_updater_diffusion_fluid_advance(fluid_species->diff_slvr,
        &app->local, fluid_species->diffD, fluid, fluid_species->cflrate, rhs);
    }
  }

  gkyl_array_reduce_range(fluid_species->omegaCfl_ptr, fluid_species->cflrate, GKYL_MAX, &app->local);

  double omegaCfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omegaCfl_ho, fluid_species->omegaCfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else {
    omegaCfl_ho[0] = fluid_species->omegaCfl_ptr[0];
  }
  omegaCfl = omegaCfl_ho[0];

  app->stat.fluid_species_rhs_tm += gkyl_time_diff_now_sec(wst);

  return app->cfl/omegaCfl;
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for fluid species
void
vm_fluid_species_apply_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species, struct gkyl_array *f)
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

      switch (fluid_species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_ABSORB:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_NO_SLIP:
          gkyl_bc_basic_advance(fluid_species->bc_lo[d], fluid_species->bc_buffer, f);
          break;
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
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_NO_SLIP:
          gkyl_bc_basic_advance(fluid_species->bc_up[d], fluid_species->bc_buffer, f);
          break;
        case GKYL_SPECIES_WEDGE:
        case GKYL_SPECIES_FIXED_FUNC:
          assert(false);
          break;
        default:
          break;
      }
    }
  }

  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, f);

  app->stat.fluid_species_bc_tm += gkyl_time_diff_now_sec(wst);
}

// release resources for fluid species
void
vm_fluid_species_release(const gkyl_vlasov_app* app, struct vm_fluid_species *f)
{
  gkyl_wv_eqn_release(f->equation);

  gkyl_array_release(f->fluid);
  gkyl_array_release(f->fluid1);
  gkyl_array_release(f->fluidnew);
  gkyl_array_release(f->bc_buffer);
  gkyl_array_release(f->cflrate);
  gkyl_array_release(f->fluid_dup);

  gkyl_dg_updater_fluid_release(f->advect_slvr);
  if (f->has_diffusion) {
    gkyl_array_release(f->diffD);
    if (f->info.diffusion.Dij)
      gkyl_dg_updater_diffusion_gen_release(f->diff_slvr_gen);
    else if (f->info.diffusion.D)
      gkyl_dg_updater_diffusion_fluid_release(f->diff_slvr);
  }

  if (f->eqn_type == GKYL_EQN_ADVECTION) {
    gkyl_array_release(f->app_advect);
    if (app->use_gpu)
      gkyl_array_release(f->app_advect_host);
  }
  else if (f->eqn_type == GKYL_EQN_EULER || f->eqn_type == GKYL_EQN_ISO_EULER) {
    gkyl_array_release(f->u);
    gkyl_array_release(f->p);
    gkyl_array_release(f->cell_avg_prim);
    gkyl_array_release(f->u_surf);
    gkyl_array_release(f->p_surf);
    gkyl_array_release(f->ke_old);
    gkyl_array_release(f->ke_new);
    gkyl_dg_calc_fluid_vars_release(f->calc_fluid_vars);
    gkyl_dg_calc_fluid_vars_release(f->calc_fluid_vars_ext);
  }

  gkyl_array_release(f->integ_mom);
  gkyl_dynvec_release(f->integ_diag);

  gkyl_array_release(f->app_accel);
  if (f->has_app_accel) {
    if (app->use_gpu) {
      gkyl_array_release(f->app_accel_host);
    }
    gkyl_proj_on_basis_release(f->app_accel_proj);
  }

  if (f->source_id) {
    vm_fluid_species_source_release(app, &f->src);
  }

  if (app->use_gpu) {
    gkyl_array_release(f->fluid_host);
    gkyl_cu_free(f->omegaCfl_ptr);
    gkyl_cu_free(f->red_integ_diag);
  }
  else {
    gkyl_free(f->omegaCfl_ptr);
  }
  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    gkyl_bc_basic_release(f->bc_lo[d]);
    gkyl_bc_basic_release(f->bc_up[d]);
  }
}
