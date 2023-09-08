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
    long vol = GKYL_MAX(app->skin_ghost.lower_skin[d].volume, app->skin_ghost.upper_skin[d].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  f->bc_buffer = mkarr(app->use_gpu, num_eqn*app->confBasis.num_basis, buff_sz);

  // allocate cflrate (scalar array)
  f->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu)
    f->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
  else
    f->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_flags[GKYL_MAX_DIM] = {0, 0, 0};

  // initial pointers to primitive variables, pressure, and boolean array for if we are only using the cell average for primitive variables
  // For isothermal Euler, prim : (ux, uy, uz), p : (vth*rho)
  // For Euler, prim : (ux, uy, uz, T/m), p : (gamma - 1)*(E - 1/2 rho u^2)
  f->prim = 0;
  f->p = 0;
  f->cell_avg_prim = 0;

  f->param = 0.0;
  // fluid solvers
  if (f->info.vt || f->info.gas_gamma) {
    if (f->info.vt) {
      f->param = f->info.vt; // parameter for isothermal Euler is vt, thermal velocity
      f->eqn_id = GKYL_EQN_ISO_EULER;
    }
    else if (f->info.gas_gamma) {
      f->param = f->info.gas_gamma; // parameter for Euler is gas_gamma, adiabatic index
      f->eqn_id = GKYL_EQN_EULER;
    }
    // allocate array to store fluid velocity (ux, uy, uz) and pressure
    // For isothermal Euler, p : (vth*rho)
    // For Euler, p : (gamma - 1)*(E - 1/2 rho u^2)
    f->prim = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
    f->p = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    // boolean array for if we are only using the cell average for primitive variables
    f->cell_avg_prim = mk_int_arr(app->use_gpu, 1, app->local_ext.volume);
    struct gkyl_dg_euler_auxfields aux_inp = {.u_i = f->prim, .p_ij = f->p};
    f->advect_slvr = gkyl_dg_updater_fluid_new(&app->grid, &app->confBasis,
      &app->local, f->eqn_id, f->param, &aux_inp, app->use_gpu);
  }
  else if (f->info.advection.velocity) {
    f->eqn_id = GKYL_EQN_ADVECTION;

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
      &app->local, f->eqn_id, f->param, &aux_inp, app->use_gpu);
  }
  else {
    f->eqn_id = GKYL_EQN_EULER_PKPM;

    f->pkpm_species = vm_find_species(app, f->info.pkpm_species);
    // index in fluid_species struct of fluid species kinetic species is colliding with
    f->species_index = vm_find_species_idx(app, f->info.pkpm_species);

    struct gkyl_dg_euler_pkpm_auxfields aux_inp = {.vlasov_pkpm_moms = f->pkpm_species->pkpm_moms.marr, 
      .pkpm_prim = f->pkpm_species->pkpm_prim, .pkpm_prim_surf = f->pkpm_species->pkpm_prim_surf, 
      .pkpm_p_ij = f->pkpm_species->pkpm_p_ij, .pkpm_p_ij_surf = f->pkpm_species->pkpm_p_ij_surf, 
      .pkpm_lax = f->pkpm_species->pkpm_lax};
    f->advect_slvr = gkyl_dg_updater_fluid_new(&app->grid, &app->confBasis,
      &app->local, f->eqn_id, f->param, &aux_inp, app->use_gpu);
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
    int num_eqn = 1;
    if (f->eqn_id == GKYL_EQN_ISO_EULER) 
      num_eqn = 4;
    else if (f->eqn_id == GKYL_EQN_EULER) 
      num_eqn = 5;
    else if (f->eqn_id == GKYL_EQN_EULER_PKPM) 
      num_eqn = 3;

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

    f->diff_slvr = gkyl_dg_updater_diffusion_fluid_new(&app->grid, &app->confBasis,
      true, num_eqn, NULL, f->info.diffusion.order, &app->local, app->use_gpu);
  }

  // array for storing integrated moments in each cell
  f->integ_mom = mkarr(app->use_gpu, 6, app->local_ext.volume);
  if (app->use_gpu) {
    f->red_integ_diag = gkyl_cu_malloc(sizeof(double[6]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments 
  f->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, 6);
  f->is_first_integ_write_call = true;

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

  int ghost[GKYL_MAX_DIM] = {0.0};
  for (int d=0; d<app->cdim; ++d)
    ghost[d] = 1;

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
    else if (f->lower_bc[d] == GKYL_SPECIES_REFLECT && f->eqn_id == GKYL_EQN_EULER_PKPM) {
      bctype = GKYL_BC_PKPM_MOM_REFLECT;
    }
    else if (f->lower_bc[d] == GKYL_SPECIES_NO_SLIP && f->eqn_id == GKYL_EQN_EULER_PKPM) {
      bctype = GKYL_BC_PKPM_MOM_NO_SLIP;
    }

    // Create local lower skin and ghost ranges
    gkyl_skin_ghost_ranges(&f->lower_skin[d], &f->lower_ghost[d], d, GKYL_LOWER_EDGE, &app->local_ext, ghost);
    f->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev.confBasis,
      &f->lower_skin[d], &f->lower_ghost[d], f->fluid->ncomp, app->cdim, app->use_gpu);

    // Upper BC updater. Copy BCs by default.
    if (f->upper_bc[d] == GKYL_SPECIES_COPY) {
      bctype = GKYL_BC_COPY;
    }
    else if (f->upper_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype = GKYL_BC_ABSORB;
      f->bc_is_absorb = true;
    }
    else if (f->upper_bc[d] == GKYL_SPECIES_REFLECT && f->eqn_id == GKYL_EQN_EULER_PKPM) {
      bctype = GKYL_BC_PKPM_MOM_REFLECT;
    }
    else if (f->upper_bc[d] == GKYL_SPECIES_NO_SLIP && f->eqn_id == GKYL_EQN_EULER_PKPM) {
      bctype = GKYL_BC_PKPM_MOM_NO_SLIP;
    }

    // Create local upper skin and ghost ranges
    gkyl_skin_ghost_ranges(&f->upper_skin[d], &f->upper_ghost[d], d, GKYL_UPPER_EDGE, &app->local_ext, ghost);
    f->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.confBasis,
      &f->upper_skin[d], &f->upper_ghost[d], f->fluid->ncomp, app->cdim, app->use_gpu);
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

  // we are pre-computing source for now as it is time-independent
  vm_fluid_species_source_calc(app, fluid_species, t0);
}

void
vm_fluid_species_prim_vars(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  const struct gkyl_array *fluid)
{
  // Once merged with main, can utilize dg_prim_vars_type infrastructure written for
  // GPU hackathon (JJ : 08/03/23)
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

  if (fluid_species->has_diffusion) {
    if (fluid_species->info.diffusion.Dij) 
      gkyl_dg_updater_diffusion_gen_advance(fluid_species->diff_slvr_gen,
        &app->local, fluid_species->diffD, fluid, fluid_species->cflrate, rhs);
    else if (fluid_species->info.diffusion.D)
      gkyl_dg_updater_diffusion_fluid_advance(fluid_species->diff_slvr,
        &app->local, fluid_species->diffD, fluid, fluid_species->cflrate, rhs);
  }

  // Accumulate source contribution if PKPM -> adds forces (E + u x B) to momentum equation RHS
  if (fluid_species->eqn_id == GKYL_EQN_EULER_PKPM) {
    double qbym = fluid_species->pkpm_species->info.charge/fluid_species->pkpm_species->info.mass;
    gkyl_array_set(fluid_species->pkpm_species->qmem, qbym, em);

    // Accumulate applied acceleration and/or q/m*(external electromagnetic)
    // fields onto qmem to get the total acceleration
    if (fluid_species->pkpm_species->has_accel)
      gkyl_array_accumulate(fluid_species->pkpm_species->qmem, 1.0, fluid_species->pkpm_species->accel);
    if (app->field->has_ext_em)
      gkyl_array_accumulate(fluid_species->pkpm_species->qmem, qbym, app->field->ext_em);

    gkyl_dg_calc_pkpm_vars_source(fluid_species->pkpm_species->calc_pkpm_vars, &app->local, 
      fluid_species->pkpm_species->qmem, fluid_species->pkpm_species->pkpm_moms.marr, fluid, rhs);
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
  gkyl_array_release(f->fluid);
  gkyl_array_release(f->fluid1);
  gkyl_array_release(f->fluidnew);
  gkyl_array_release(f->bc_buffer);
  gkyl_array_release(f->cflrate);

  gkyl_dg_updater_fluid_release(f->advect_slvr);
  if (f->has_diffusion) {
    gkyl_array_release(f->diffD);
    if (f->info.diffusion.Dij)
      gkyl_dg_updater_diffusion_gen_release(f->diff_slvr_gen);
    else if (f->info.diffusion.D)
      gkyl_dg_updater_diffusion_fluid_release(f->diff_slvr);
  }

  if (f->eqn_id == GKYL_EQN_ADVECTION) {
    gkyl_array_release(f->app_advect);
    if (app->use_gpu)
      gkyl_array_release(f->app_advect_host);
  }
  else if (f->eqn_id == GKYL_EQN_EULER || f->eqn_id == GKYL_EQN_ISO_EULER) {
    gkyl_array_release(f->prim);
    gkyl_array_release(f->p);
    gkyl_array_release(f->cell_avg_prim);
  }

  gkyl_array_release(f->integ_mom);
  gkyl_dynvec_release(f->integ_diag);

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
