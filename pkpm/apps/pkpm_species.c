#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_pkpm_priv.h>

#include <assert.h>
#include <time.h>

// initialize species object
void
pkpm_species_init(struct gkyl_pkpm *pkpm, struct gkyl_pkpm_app *app, struct pkpm_species *s)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim + vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d = 0; d < cdim; ++d) {
    cells[d] = pkpm->cells[d];
    lower[d] = pkpm->lower[d];
    upper[d] = pkpm->upper[d];
    ghost[d] = 1;
  }
  for (int d = 0; d < vdim; ++d) {
    // full phase space grid
    cells[cdim + d] = s->info.cells[d];
    lower[cdim + d] = s->info.lower[d];
    upper[cdim + d] = s->info.upper[d];
    ghost[cdim + d] = 0; // no ghost-cells in velocity space

    // only velocity space
    cells_vel[d] = s->info.cells[d];
    lower_vel[d] = s->info.lower[d];
    upper_vel[d] = s->info.upper[d];
    ghost_vel[d] = 0; // no ghost-cells in velocity space
  }
  // full phase space grid
  gkyl_rect_grid_init(&s->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->global_ext, &s->global);

  // velocity space grid
  gkyl_rect_grid_init(&s->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&s->grid_vel, ghost_vel, &s->local_ext_vel, &s->local_vel);

  // phase-space communicator
  s->comm = gkyl_comm_extend_comm(app->comm, &s->local_vel);

  // create local and local_ext from app local range
  struct gkyl_range local;
  // local = conf-local X local_vel
  gkyl_range_ten_prod(&local, &app->local, &s->local_vel);
  gkyl_create_ranges(&local, ghost, &s->local_ext, &s->local);

  // allocate distribution function arrays
  // PKPM has two distribution functions F_0 and T_perp*G, first two Laguerre moments
  s->f = mkarr(app->use_gpu, 2 * app->basis.num_basis, s->local_ext.volume);
  s->f1 = mkarr(app->use_gpu, 2 * app->basis.num_basis, s->local_ext.volume);
  s->fnew = mkarr(app->use_gpu, 2 * app->basis.num_basis, s->local_ext.volume);
  // allocate momentum arrays (rho ux, rho uy, rho uz)
  s->fluid = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, app->local_ext.volume);
  s->fluid1 = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, app->local_ext.volume);
  s->fluidnew = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, app->local_ext.volume);

  // create host arrays if on GPUs so initialization occurs host-side
  s->f_host = s->f;
  s->fluid_host = s->fluid;
  if (app->use_gpu) {
    s->f_host = mkarr(false, 2 * app->basis.num_basis, s->local_ext.volume);
    s->fluid_host = mkarr(false, 3 * app->confBasis.num_basis, app->local_ext.volume);
  }

  // Duplicate copy of momentum data in case time step fails.
  // Needed because of implicit source split which modifies solution and
  // is always successful, so if a time step fails due to the SSP RK3
  // we must restore the old solution before restarting the time step
  s->fluid_dup = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, app->local_ext.volume);

  // Wave equation object for upwinding fluid equations
  // We use the 10 moment system since PKPM model generates a full pressure tensor
  // k0 = 0.0 and use_grad_closure = false because we do not need a closure
  // (the kinetic equations are the closure)
  s->equation = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, app->use_gpu);

  // allocate cflrate for distribution function and momentum updates (scalar arrays)
  s->cflrate_f = mkarr(app->use_gpu, 1, s->local_ext.volume);
  s->cflrate_fluid = mkarr(app->use_gpu, 1, app->local_ext.volume);

  if (app->use_gpu) {
    s->omegaCfl_ptr_dist = gkyl_cu_malloc(sizeof(double));
    s->omegaCfl_ptr_fluid = gkyl_cu_malloc(sizeof(double));
  } else {
    s->omegaCfl_ptr_dist = gkyl_malloc(sizeof(double));
    s->omegaCfl_ptr_fluid = gkyl_malloc(sizeof(double));
  }

  // allocate array to store q/m*(E,B) if using a fully explicit update
  s->qmem = mkarr(app->use_gpu, 8 * app->confBasis.num_basis, app->local_ext.volume);

  int nc = app->confBasis.num_basis;
  // pkpm moments for the update [rho, p_par, p_perp, M1]
  pkpm_species_moment_init(app, s, &s->pkpm_moms, false);
  // pkpm moments for diagnostics [rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp]
  // (the is_integrated flag selects the diagnostic moment set)
  pkpm_species_moment_init(app, s, &s->pkpm_moms_diag, true);

  // Refactored update variables (p=1 configuration space only).
  s->pkpm_u = mkarr(app->use_gpu, 3 * nc, app->local_ext.volume);
  s->pkpm_grad_u = mkarr(app->use_gpu, cdim * 3 * nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_grad_u, 0.0);
  s->pkpm_theta = mkarr(app->use_gpu, nc, app->local_ext.volume);
  s->pkpm_a_P = mkarr(app->use_gpu, nc, app->local_ext.volume);
  s->pkpm_rhs_moms = mkarr(app->use_gpu, 4 * nc, app->local_ext.volume);
  s->pkpm_p_ij = mkarr(app->use_gpu, 6 * nc, app->local_ext.volume);
  s->pkpm_nu_vthsq_sum = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_nu_vthsq_sum, 0.0);
  s->cell_avg_prim = mk_int_arr(app->use_gpu, 1, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_a_P, 0.0);
  gkyl_array_clear(s->pkpm_theta, 0.0);
  s->num_null_rho = 0;

  // Check if limiter_fac is specified for adjusting how much diffusion is applied through slope limiter
  // If not specified, set to 0.0 and updater sets default behavior (1/sqrt(3); see gkyl_dg_calc_pkpm_vars.h)
  double limiter_fac = s->info.limiter_fac == 0 ? 0.0 : s->info.limiter_fac;
  s->limit_fluid = s->info.limit_fluid;

  // Updaters for the primitive variables (weak division u, nodal theta, the a_P
  // closure, limiter pressure projection, integrated variables, sources, diagnostics).
  // Two instances, one over the extended range and one over the local range, for ease of handling boundary conditions.
  s->calc_pkpm_vars_ext = gkyl_dg_calc_pkpm_vars_new(
    &app->grid, &app->confBasis, &app->local_ext, s->equation, app->geom, limiter_fac, app->use_gpu
  );
  s->calc_pkpm_vars = gkyl_dg_calc_pkpm_vars_new(
    &app->grid, &app->confBasis, &app->local, s->equation, app->geom, limiter_fac, app->use_gpu
  );

  // array for storing integrated fluid variables in each cell
  // integ_pkpm_mom : integral[rho, rhoux, rhouy, rhouz, rhoux^2, rhouy^2, rhouz^2, p_par, p_perp]
  s->integ_pkpm_mom = mkarr(app->use_gpu, 9, app->local_ext.volume);
  // arrays for I/O, fluid_io and pkpm_vars_io
  // fluid_io : [rho, rhoux, rhouy, rhouz, P_xx+rhoux^2, P_xy+rhouxuy, P_xz+rhouxuz, P_yy+rhouy^2, P_yz+rhouyuz, P_zz+rhouz^2]
  // pkpm_vars_io : [ux, uy, uz, T_perp/m, a_P, div(b), bb:grad(u), T_perp/m div(b)]
  s->fluid_io = mkarr(app->use_gpu, 10 * nc, app->local_ext.volume);
  s->pkpm_vars_io = mkarr(app->use_gpu, 8 * nc, app->local_ext.volume);
  s->fluid_io_host = s->fluid_io;
  s->pkpm_vars_io_host = s->pkpm_vars_io;
  if (app->use_gpu) {
    s->fluid_io_host = mkarr(false, 10 * nc, app->local_ext.volume);
    s->pkpm_vars_io_host = mkarr(false, 8 * nc, app->local_ext.volume);
  }

  // Surface-flux updaters and stored modal fluxes (each cell owns its lower faces).
  // The mass-flux constraint leaves the domain faces of non-periodic directions alone:
  // there the ghost state is a boundary construct and the kinetic flux is the physical one.
  struct gkyl_dg_pkpm_conf_flux_surf_inp conf_flux_inp = {
    .phase_grid = &s->grid,
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .acoustic_floor = false,
    .use_gpu = app->use_gpu,
  };
  for (int d = 0; d < cdim; ++d) {
    conf_flux_inp.skip_domain_faces[d] = true;
  }
  for (int d = 0; d < app->num_periodic_dir; ++d) {
    conf_flux_inp.skip_domain_faces[app->periodic_dirs[d]] = false;
  }
  s->calc_conf_flux = gkyl_dg_pkpm_conf_flux_surf_inew(&conf_flux_inp);
  int nkc = gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(s->calc_conf_flux);
  s->conf_flux_surf = mkarr(app->use_gpu, 2 * cdim * nkc, s->local_ext.volume);
  int nmf = gkyl_dg_pkpm_conf_flux_surf_num_face_basis(s->calc_conf_flux);
  s->pkpm_kin_mass_flux = mkarr(app->use_gpu, cdim * nmf, s->local_ext.volume);
  s->pkpm_kin_mass_flux_face = mkarr(app->use_gpu, cdim * nmf, app->local_ext.volume);
  s->pkpm_roe_mass_flux = mkarr(app->use_gpu, cdim * nmf, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_roe_mass_flux, 0.0);

  struct gkyl_dg_pkpm_vel_flux_surf_inp vel_flux_inp = {
    .phase_grid = &s->grid,
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .use_gpu = app->use_gpu,
  };
  s->calc_vel_flux = gkyl_dg_pkpm_vel_flux_surf_inew(&vel_flux_inp);
  int nkv = gkyl_dg_pkpm_vel_flux_surf_num_surf_basis(s->calc_vel_flux);
  s->vel_flux_surf = mkarr(app->use_gpu, 2 * nkv, s->local_ext.volume);

  struct gkyl_dg_pkpm_fluid_flux_surf_inp fluid_flux_inp = {
    .conf_grid = &app->grid,
    .conf_basis = &app->confBasis,
    .wv_eqn = s->equation,
    .geom = app->geom,
    .use_gpu = app->use_gpu,
  };
  s->calc_fluid_flux = gkyl_dg_pkpm_fluid_flux_surf_inew(&fluid_flux_inp);
  int nkf = gkyl_dg_pkpm_fluid_flux_surf_num_surf_basis(s->calc_fluid_flux);
  assert(nkf == nmf); // the two flux updaters share the configuration face basis
  s->fluid_flux_surf = mkarr(app->use_gpu, 3 * cdim * nkf, app->local_ext.volume);

  // Collision coefficients consumed by the configuration pass (G source): copied from
  // the LBO arrays at each stage when collisions are on, zeros otherwise.
  s->pkpm_nu_sum = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_clear(s->pkpm_nu_sum, 0.0);

  // Equation objects and hyper_dg solvers of the three passes.
  struct gkyl_dg_pkpm_kinetic_conf_inp kinetic_conf_inp = {
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .conf_range = &app->local,
    .phase_range = &s->local,
    .auxfields =
      {
        .bvar = app->field->bvar,
        .grad_b = app->field->grad_b,
        .u = s->pkpm_u,
        .grad_u = s->pkpm_grad_u,
        .theta = s->pkpm_theta,
        .nu_sum = s->pkpm_nu_sum,
        .nu_vthsq_sum = s->pkpm_nu_vthsq_sum,
        .conf_flux_surf = s->conf_flux_surf,
      },
    .use_gpu = app->use_gpu,
  };
  s->eqn_kinetic_conf = gkyl_dg_pkpm_kinetic_conf_inew(&kinetic_conf_inp);
  struct gkyl_dg_pkpm_kinetic_vel_inp kinetic_vel_inp = {
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .conf_range = &app->local,
    .phase_range = &s->local,
    .auxfields = {.a_P = s->pkpm_a_P, .vel_flux_surf = s->vel_flux_surf},
    .use_gpu = app->use_gpu,
  };
  s->eqn_kinetic_vel = gkyl_dg_pkpm_kinetic_vel_inew(&kinetic_vel_inp);
  struct gkyl_dg_pkpm_fluid_inp fluid_inp = {
    .conf_basis = &app->confBasis,
    .conf_range = &app->local,
    .auxfields =
      {
        .moms = s->pkpm_moms.marr,
        .bvar = app->field->bvar,
        .u = s->pkpm_u,
        .fluid_flux_surf = s->fluid_flux_surf,
      },
    .use_gpu = app->use_gpu,
  };
  s->eqn_fluid = gkyl_dg_pkpm_fluid_inew(&fluid_inp);

  int up_dirs_conf[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_conf[2 * GKYL_MAX_DIM] = {0};
  s->slvr_kinetic_conf = gkyl_hyper_dg_new(
    &s->grid, &app->basis, s->eqn_kinetic_conf, cdim, up_dirs_conf, zero_flux_conf, 1, app->use_gpu
  );
  int up_dirs_vel[GKYL_MAX_DIM] = {cdim}, zero_flux_vel[2 * GKYL_MAX_DIM] = {0};
  zero_flux_vel[cdim] = zero_flux_vel[cdim + pdim] = 1; // zero flux at the velocity domain edges
  s->slvr_kinetic_vel = gkyl_hyper_dg_new(
    &s->grid, &app->basis, s->eqn_kinetic_vel, 1, up_dirs_vel, zero_flux_vel, 1, app->use_gpu
  );
  s->slvr_fluid = gkyl_hyper_dg_new(
    &app->grid, &app->confBasis, s->eqn_fluid, cdim, up_dirs_conf, zero_flux_conf, 1, app->use_gpu
  );
  s->kinetic_tm = 0.0;
  s->fluid_tm = 0.0;

  // array for storing F_0^2 (0th Laguerre coefficient) in each cell
  s->L2_f = mkarr(app->use_gpu, 1, s->local_ext.volume);
  if (app->use_gpu) {
    s->red_L2_f = gkyl_cu_malloc(sizeof(double));
    s->red_integ_diag = gkyl_cu_malloc(sizeof(double[9]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments and F_0^2
  s->integ_L2_f = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  s->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, 9);
  s->is_first_integ_L2_write_call = true;
  s->is_first_integ_write_call = true;

  // Initialize applied acceleration for use in force update.
  s->app_accel = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(s->app_accel, 0.0);
  s->has_app_accel = false;
  s->app_accel_evolve = false;
  // setup applied acceleration
  if (s->info.app_accel) {
    s->has_app_accel = true;
    if (s->info.app_accel_evolve) {
      s->app_accel_evolve = s->info.app_accel_evolve;
    }

    s->app_accel_host = s->app_accel;
    if (app->use_gpu) {
      s->app_accel_host = mkarr(false, 3 * app->confBasis.num_basis, app->local_ext.volume);
    }
    s->app_accel_proj = gkyl_proj_on_basis_new(
      &app->grid, &app->confBasis, app->confBasis.poly_order + 1, 3, s->info.app_accel,
      s->info.app_accel_ctx
    );
  }

  // determine collision type to use in PKPM update
  s->collision_id = s->info.collisions.collision_id;
  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    pkpm_species_lbo_init(app, s, &s->lbo);
  }

  // initialize diffusion if present
  s->has_diffusion = false;
  if (s->info.diffusion.D) {
    s->has_diffusion = true;
    s->info.diffusion.order = s->info.diffusion.order < 2 ? 2 : s->info.diffusion.order;
    int num_eqn = 3;

    int szD = cdim;
    s->diffD = mkarr(app->use_gpu, szD, 1);
    struct gkyl_array *diffD_host = s->diffD;
    if (app->use_gpu) {
      diffD_host = mkarr(false, szD, 1);
    }
    // Set diffusion coefficient in each direction to input value.
    gkyl_array_clear(diffD_host, 0.);
    for (int d = 0; d < cdim; d++) {
      gkyl_array_shiftc(diffD_host, s->info.diffusion.D, d);
    }

    if (app->use_gpu) { // note: diffD_host is same as diffD when not on GPUs
      gkyl_array_copy(s->diffD, diffD_host);
      gkyl_array_release(diffD_host);
    }

    bool is_zero_flux[GKYL_MAX_CDIM] = {false};
    s->diff_slvr = gkyl_dg_updater_diffusion_fluid_new(
      &app->grid, &app->confBasis, num_eqn, NULL, s->info.diffusion.order, is_zero_flux,
      app->use_gpu
    );
  }

  // determine which directions are not periodic and get non-periodic BC info
  // also create local skin-ghost ranges for ease of handling non-periodic BCs
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np[app->periodic_dirs[d]] = 0;
  }

  for (int dir = 0; dir < cdim; ++dir) {
    s->lower_bc[dir] = s->upper_bc[dir] = GKYL_SPECIES_COPY;
    if (is_np[dir]) {
      const enum gkyl_species_bc_type *bc;
      if (dir == 0) {
        bc = s->info.bcx;
      } else if (dir == 1) {
        bc = s->info.bcy;
      } else {
        bc = s->info.bcz;
      }

      s->lower_bc[dir] = bc[0];
      s->upper_bc[dir] = bc[1];
    }
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(
      &s->lower_skin_dist[dir], &s->lower_ghost_dist[dir], dir, GKYL_LOWER_EDGE, &s->local_ext,
      ghost
    );
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(
      &s->upper_skin_dist[dir], &s->upper_ghost_dist[dir], dir, GKYL_UPPER_EDGE, &s->local_ext,
      ghost
    );
  }

  // allocate buffer for applying BCs
  long buff_sz_dist = 0;
  long buff_sz_fluid = 0;
  // compute buffer size needed
  for (int dir = 0; dir < cdim; ++dir) {
    long vol_dist = GKYL_MAX2(s->lower_skin_dist[dir].volume, s->upper_skin_dist[dir].volume);
    buff_sz_dist = buff_sz_dist > vol_dist ? buff_sz_dist : vol_dist;
    long vol_fluid = GKYL_MAX2(app->lower_skin[dir].volume, app->upper_skin[dir].volume);
    buff_sz_fluid = buff_sz_fluid > vol_fluid ? buff_sz_fluid : vol_fluid;
  }

  s->bc_buffer_dist = mkarr(app->use_gpu, 2 * app->basis.num_basis, buff_sz_dist);
  // buffer arrays for fixed function boundary conditions on distribution function
  s->bc_buffer_lo_fixed_dist = mkarr(app->use_gpu, 2 * app->basis.num_basis, buff_sz_dist);
  s->bc_buffer_up_fixed_dist = mkarr(app->use_gpu, 2 * app->basis.num_basis, buff_sz_dist);

  s->bc_buffer_fluid = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, buff_sz_fluid);
  // buffer arrays for fixed function boundary conditions on momentum
  s->bc_buffer_lo_fixed_fluid = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, buff_sz_fluid);
  s->bc_buffer_up_fixed_fluid = mkarr(app->use_gpu, 3 * app->confBasis.num_basis, buff_sz_fluid);

  // Certain operations fail if absorbing BCs used because absorbing BCs
  // means the mass density is 0 in the ghost cells (divide by zero)
  s->bc_is_absorb = false;
  for (int d = 0; d < cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype_dist = GKYL_BC_COPY;
    enum gkyl_bc_basic_type bctype_fluid = GKYL_BC_COPY;
    if (s->lower_bc[d] == GKYL_SPECIES_COPY) {
      bctype_dist = GKYL_BC_COPY;
      bctype_fluid = GKYL_BC_COPY;
    } else if (s->lower_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype_dist = GKYL_BC_ABSORB;
      bctype_fluid = GKYL_BC_ABSORB;
      s->bc_is_absorb = true;
    } else if (s->lower_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype_dist = GKYL_BC_PKPM_SPECIES_REFLECT;
      bctype_fluid = GKYL_BC_PKPM_MOM_REFLECT;
    } else if (s->lower_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      bctype_dist = GKYL_BC_FIXED_FUNC;
      bctype_fluid = GKYL_BC_FIXED_FUNC;
    }

    // Distribution function non-periodic lower boundary conditions
    s->bc_lo_dist[d] = gkyl_bc_basic_new(
      d, GKYL_LOWER_EDGE, bctype_dist, app->basis_on_dev.basis, &s->lower_skin_dist[d],
      &s->lower_ghost_dist[d], s->f->ncomp, app->cdim, app->use_gpu
    );
    // Momentum non-periodic lower boundary conditions
    s->bc_lo_fluid[d] = gkyl_bc_basic_new(
      d, GKYL_LOWER_EDGE, bctype_fluid, app->basis_on_dev.confBasis, &app->lower_skin[d],
      &app->lower_ghost[d], s->fluid->ncomp, app->cdim, app->use_gpu
    );

    // Upper BC updater. Copy BCs by default.
    if (s->upper_bc[d] == GKYL_SPECIES_COPY) {
      bctype_dist = GKYL_BC_COPY;
      bctype_fluid = GKYL_BC_COPY;
    } else if (s->upper_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype_dist = GKYL_BC_ABSORB;
      bctype_fluid = GKYL_BC_ABSORB;
      s->bc_is_absorb = true;
    } else if (s->upper_bc[d] == GKYL_SPECIES_REFLECT) {
      bctype_dist = GKYL_BC_PKPM_SPECIES_REFLECT;
      bctype_fluid = GKYL_BC_PKPM_MOM_REFLECT;
    } else if (s->upper_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      bctype_dist = GKYL_BC_FIXED_FUNC;
      bctype_fluid = GKYL_BC_FIXED_FUNC;
    }

    // Distribution function non-periodic upper boundary conditions
    s->bc_up_dist[d] = gkyl_bc_basic_new(
      d, GKYL_UPPER_EDGE, bctype_dist, app->basis_on_dev.basis, &s->upper_skin_dist[d],
      &s->upper_ghost_dist[d], s->f->ncomp, app->cdim, app->use_gpu
    );
    // Momentum non-periodic upper boundary conditions
    s->bc_up_fluid[d] = gkyl_bc_basic_new(
      d, GKYL_UPPER_EDGE, bctype_fluid, app->basis_on_dev.confBasis, &app->upper_skin[d],
      &app->upper_ghost[d], s->fluid->ncomp, app->cdim, app->use_gpu
    );
  }
}

void
pkpm_species_apply_ic(gkyl_pkpm_app *app, struct pkpm_species *species, double t0)
{
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj_dist;
  proj_dist = gkyl_proj_on_basis_new(
    &species->grid, &app->basis, 8, 2, species->info.init_dist, species->info.ctx_dist
  );
  gkyl_proj_on_basis *proj_fluid;
  proj_fluid = gkyl_proj_on_basis_new(
    &app->grid, &app->confBasis, poly_order + 1, 3, species->info.init_fluid,
    species->info.ctx_fluid
  );

  // run updaters; need to project onto extended range for ease of handling
  // subsequent operations over extended range such as primitive variable computations
  // This is needed to fill the corner cells as the corner cells may not be filled by
  // boundary conditions and we cannot divide by 0 anywhere or the weak divisions will fail
  gkyl_proj_on_basis_advance(proj_dist, t0, &species->local_ext, species->f_host);
  gkyl_proj_on_basis_release(proj_dist);
  gkyl_proj_on_basis_advance(proj_fluid, t0, &app->local_ext, species->fluid_host);
  gkyl_proj_on_basis_release(proj_fluid);

  // note: f_host and fluid_host are the same as f and fluid respectively when not on GPUs
  if (app->use_gpu) {
    gkyl_array_copy(species->f, species->f_host);
    gkyl_array_copy(species->fluid, species->fluid_host);
  }

  // we are pre-computing acceleration for now in case is time-independent
  pkpm_species_calc_app_accel(app, species, t0);

  // copy contents of initial conditions into buffer if specific BCs require them
  // *only works in x dimension for now*
  gkyl_bc_basic_buffer_fixed_func(
    species->bc_lo_dist[0], species->bc_buffer_lo_fixed_dist, species->f
  );
  gkyl_bc_basic_buffer_fixed_func(
    species->bc_up_dist[0], species->bc_buffer_up_fixed_dist, species->f
  );
  gkyl_bc_basic_buffer_fixed_func(
    species->bc_lo_fluid[0], species->bc_buffer_lo_fixed_fluid, species->fluid
  );
  gkyl_bc_basic_buffer_fixed_func(
    species->bc_up_fluid[0], species->bc_buffer_up_fixed_fluid, species->fluid
  );
}

void
pkpm_species_calc_app_accel(gkyl_pkpm_app *app, struct pkpm_species *species, double tm)
{
  if (species->has_app_accel) {
    gkyl_proj_on_basis_advance(
      species->app_accel_proj, tm, &app->local_ext, species->app_accel_host
    );
    if (app->use_gpu) {
      // note: app_accel_host is same as app_accel when not on GPUs
      gkyl_array_copy(species->app_accel, species->app_accel_host);
    }
  }
}

void
pkpm_species_calc_pkpm_vars(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  struct timespec tm = gkyl_wall_clock();

  // Compute rho, p_par, p_perp, M1 over the extended range (ghost cells hold valid
  // data after the boundary conditions), then u by weak division and theta = p_perp/rho
  // by the nodal quotient. With absorbing BCs rho vanishes in the ghost cells, so the
  // quotients are restricted to the local range.
  pkpm_species_moment_calc(&species->pkpm_moms, species->local_ext, app->local_ext, fin);
  if (species->bc_is_absorb) {
    gkyl_dg_calc_pkpm_vars_u(
      species->calc_pkpm_vars, species->pkpm_moms.marr, fluidin, species->cell_avg_prim,
      species->pkpm_u
    );
    gkyl_dg_calc_pkpm_vars_theta(
      species->calc_pkpm_vars, &app->local, species->pkpm_moms.marr, species->pkpm_theta
    );
  } else {
    gkyl_dg_calc_pkpm_vars_u(
      species->calc_pkpm_vars_ext, species->pkpm_moms.marr, fluidin, species->cell_avg_prim,
      species->pkpm_u
    );
    gkyl_dg_calc_pkpm_vars_theta(
      species->calc_pkpm_vars_ext, &app->local_ext, species->pkpm_moms.marr, species->pkpm_theta
    );
  }
  // Weak gradient of u for the kinetic sources (bb:grad u, div u); needs u in the ghost
  // cells (with absorbing boundaries the ghost u is the last value computed there).
  gkyl_dg_calc_pkpm_vars_grad_u(species->calc_pkpm_vars, &app->local, species->pkpm_u, species->pkpm_grad_u);

  app->stat.species_pkpm_vars_tm += gkyl_time_diff_now_sec(tm);
}

// Face fluxes of the momentum equation (Roe on the nodal states), including the mass
// flux of the same Riemann solve to which the kinetic configuration-space flux is
// constrained. Needs the stage moments, b and u.
static void
pkpm_species_fluid_flux(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fluidin
)
{
  gkyl_dg_pkpm_fluid_flux_surf_advance(
    species->calc_fluid_flux, &app->local, &app->local_ext, species->pkpm_moms.marr,
    app->field->bvar, species->pkpm_u, fluidin, species->cflrate_fluid, species->fluid_flux_surf,
    species->pkpm_roe_mass_flux
  );
}

// Configuration-space pass of the kinetic update into rhs (cleared here), followed by
// the closure of the pressure force from the first moment of that partial right-hand side.
// The kinetic face fluxes are constrained so that the density flux equals the mass flux
// of the momentum solve (pkpm_species_fluid_flux must have run for this stage).
static void
pkpm_species_kinetic_conf_pass(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs
)
{
  gkyl_array_clear(rhs, 0.0);
  gkyl_dg_pkpm_conf_flux_surf_advance(
    species->calc_conf_flux, &app->local, &species->local, &species->local_ext, app->field->bvar,
    species->pkpm_u, species->pkpm_theta, fin, species->cflrate_f, species->conf_flux_surf,
    species->pkpm_kin_mass_flux
  );
  gkyl_dg_pkpm_conf_flux_surf_constrain(
    species->calc_conf_flux, &app->local, &species->local, &species->local_ext,
    species->info.mass, species->pkpm_moms.marr, species->pkpm_roe_mass_flux,
    species->pkpm_kin_mass_flux, species->pkpm_kin_mass_flux_face, fin, species->conf_flux_surf
  );
  gkyl_hyper_dg_advance(species->slvr_kinetic_conf, &species->local, fin, species->cflrate_f, rhs);

  // Closure: rho a_P = -M1(rhs_A) weakly, from the moments of the partial right-hand side.
  gkyl_array_clear(species->pkpm_rhs_moms, 0.0);
  pkpm_species_moment_calc_into(
    &species->pkpm_moms, species->local, app->local, rhs, species->pkpm_rhs_moms
  );
  species->num_null_rho += gkyl_dg_calc_pkpm_vars_p_force(
    species->calc_pkpm_vars, &app->local, species->pkpm_moms.marr, species->pkpm_rhs_moms,
    species->pkpm_a_P
  );
}

void
pkpm_species_calc_diag_vars(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  pkpm_species_calc_pkpm_vars(app, species, fin, fluidin);
  // The configuration pass needs the collision coefficients of this stage.
  if (species->collision_id == GKYL_LBO_COLLISIONS) {
    gkyl_array_set(species->pkpm_nu_sum, 1.0, species->lbo.nu_sum);
    gkyl_array_set_offset(
      species->pkpm_nu_vthsq_sum, 1.0, species->lbo.nu_prim_moms, app->vdim * app->confBasis.num_basis
    );
  }
  gkyl_array_clear(species->cflrate_f, 0.0);
  gkyl_array_clear(species->cflrate_fluid, 0.0);
  pkpm_species_fluid_flux(app, species, fluidin);
  // fnew is free scratch between time steps.
  pkpm_species_kinetic_conf_pass(app, species, fin, species->fnew);
}

void
pkpm_fluid_species_limiter(
  gkyl_pkpm_app *app, struct pkpm_species *species, struct gkyl_array *fin, struct gkyl_array *fluid
)
{
  if (species->limit_fluid) {
    struct timespec tm = gkyl_wall_clock();

    // Compute the PKPM moments from the kinetic equation and the flow velocity
    pkpm_species_moment_calc(&species->pkpm_moms, species->local, app->local, fin);
    gkyl_dg_calc_pkpm_vars_u(
      species->calc_pkpm_vars, species->pkpm_moms.marr, fluid, species->cell_avg_prim,
      species->pkpm_u
    );
    // The characteristic limiter consumes a p=1 projection of the pressure tensor,
    // built from the same stage moments and field being limited (limiter scratch only).
    gkyl_dg_calc_pkpm_vars_pressure_proj(
      species->calc_pkpm_vars, &app->local_ext, species->pkpm_moms.marr, app->field->bvar,
      species->pkpm_p_ij
    );

    // Limit the slopes of the solution of the fluid system
    gkyl_dg_calc_pkpm_vars_limiter(
      species->calc_pkpm_vars, &app->local, species->pkpm_u, species->pkpm_moms.marr,
      species->pkpm_p_ij, fluid
    );

    app->stat.species_pkpm_vars_tm += gkyl_time_diff_now_sec(tm);

    // Apply boundary conditions after limiting solution
    pkpm_fluid_species_apply_bc(app, species, fluid);
  }
}

// Compute the RHS for species update, returning maximum stable
// time-step.
double
pkpm_species_rhs(
  gkyl_pkpm_app *app, struct pkpm_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, const struct gkyl_array *em, struct gkyl_array *rhs_f,
  struct gkyl_array *rhs_fluid
)
{
  gkyl_array_clear(species->cflrate_f, 0.0);
  gkyl_array_clear(species->cflrate_fluid, 0.0);
  gkyl_array_clear(rhs_fluid, 0.0);

  // Collision coefficients of this stage for the G source (zeros without collisions).
  if (species->collision_id == GKYL_LBO_COLLISIONS) {
    gkyl_array_set(species->pkpm_nu_sum, 1.0, species->lbo.nu_sum);
    gkyl_array_set_offset(
      species->pkpm_nu_vthsq_sum, 1.0, species->lbo.nu_prim_moms, app->vdim * app->confBasis.num_basis
    );
  }

  // Momentum face fluxes first: the kinetic configuration pass is constrained to
  // their mass flux.
  struct timespec tm = gkyl_wall_clock();
  pkpm_species_fluid_flux(app, species, fluidin);
  species->fluid_tm += gkyl_time_diff_now_sec(tm);

  tm = gkyl_wall_clock();
  // Kinetic equations: configuration pass (clears rhs_f), closure of a_P from its
  // first moment, then the velocity pass accumulating into the same rhs_f.
  pkpm_species_kinetic_conf_pass(app, species, fin, rhs_f);
  gkyl_dg_pkpm_vel_flux_surf_advance(
    species->calc_vel_flux, &app->local, &species->local, app->field->bvar, app->field->grad_b,
    species->pkpm_u, species->pkpm_grad_u, species->pkpm_theta, species->pkpm_a_P, fin,
    species->cflrate_f, species->vel_flux_surf
  );
  gkyl_hyper_dg_advance(species->slvr_kinetic_vel, &species->local, fin, species->cflrate_f, rhs_f);
  species->kinetic_tm += gkyl_time_diff_now_sec(tm);

  if (species->collision_id == GKYL_LBO_COLLISIONS) {
    pkpm_species_lbo_rhs(app, species, &species->lbo, fin, rhs_f);
  }

  // Momentum equation: volume and surface terms from the stored face fluxes.
  tm = gkyl_wall_clock();
  gkyl_hyper_dg_advance(species->slvr_fluid, &app->local, fluidin, species->cflrate_fluid, rhs_fluid);
  species->fluid_tm += gkyl_time_diff_now_sec(tm);

  if (species->has_diffusion) {
    gkyl_dg_updater_diffusion_fluid_advance(
      species->diff_slvr, &app->local, species->diffD, fluidin, species->cflrate_fluid, rhs_fluid
    );
  }

  // If PKPM update is fully explicit, include the fluid-EM coupling in fluid update
  if (app->use_explicit_source) {
    double qbym = species->info.charge / species->info.mass;
    gkyl_array_set(species->qmem, qbym, em);

    // Accumulate applied acceleration and/or q/m*(external electromagnetic)
    // fields onto qmem to get the total acceleration
    if (species->has_app_accel) {
      gkyl_array_accumulate_range(species->qmem, 1.0, species->app_accel, &app->local);
    }
    if (app->field->has_ext_em) {
      gkyl_array_accumulate_range(species->qmem, qbym, app->field->ext_em, &app->local);
    }

    gkyl_dg_calc_pkpm_vars_source(
      species->calc_pkpm_vars, &app->local, species->qmem, species->pkpm_moms.marr, fluidin,
      rhs_fluid
    );
  }

  app->stat.n_species_omega_cfl += 1;
  tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omegaCfl_ptr_dist, species->cflrate_f, GKYL_MAX, &species->local);
  gkyl_array_reduce_range(
    species->omegaCfl_ptr_fluid, species->cflrate_fluid, GKYL_MAX, &app->local
  );

  double omegaCfl_ho_dist[1];
  double omegaCfl_ho_fluid[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omegaCfl_ho_dist, species->omegaCfl_ptr_dist, sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(
      omegaCfl_ho_fluid, species->omegaCfl_ptr_fluid, sizeof(double), GKYL_CU_MEMCPY_D2H
    );
  } else {
    omegaCfl_ho_dist[0] = species->omegaCfl_ptr_dist[0];
    omegaCfl_ho_fluid[0] = species->omegaCfl_ptr_fluid[0];
  }
  double omegaCfl = fmax(omegaCfl_ho_dist[0], omegaCfl_ho_fluid[0]);

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);

  return app->cfl / omegaCfl;
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for distribution function
void
pkpm_species_apply_bc(gkyl_pkpm_app *app, const struct pkpm_species *species, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();

  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(
    species->comm, &species->local, &species->local_ext, num_periodic_dir, app->periodic_dirs, f
  );

  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np_bc[app->periodic_dirs[d]] = 0;
  }

  for (int d = 0; d < cdim; ++d) {
    if (is_np_bc[d]) {
      switch (species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo_dist[d], species->bc_buffer_dist, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_lo_dist[d], species->bc_buffer_lo_fixed_dist, f);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }

      switch (species->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_up_dist[d], species->bc_buffer_dist, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_up_dist[d], species->bc_buffer_up_fixed_dist, f);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }
    }
  }

  gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext, f);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for distribution function
void
pkpm_fluid_species_apply_bc(
  gkyl_pkpm_app *app, const struct pkpm_species *species, struct gkyl_array *fluid
)
{
  struct timespec wst = gkyl_wall_clock();

  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(
    app->comm, &app->local, &app->local_ext, num_periodic_dir, app->periodic_dirs, fluid
  );

  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np_bc[app->periodic_dirs[d]] = 0;
  }

  for (int d = 0; d < cdim; ++d) {
    if (is_np_bc[d]) {
      switch (species->lower_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo_fluid[d], species->bc_buffer_fluid, fluid);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_lo_fluid[d], species->bc_buffer_lo_fixed_fluid, fluid);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }

      switch (species->upper_bc[d]) {
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_up_fluid[d], species->bc_buffer_fluid, fluid);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_up_fluid[d], species->bc_buffer_up_fixed_fluid, fluid);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }
    }
  }

  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, fluid);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

void
pkpm_species_calc_L2(gkyl_pkpm_app *app, double tm, const struct pkpm_species *species)
{
  gkyl_dg_calc_l2_range(&app->basis, 0, species->L2_f, 0, species->f, species->local);
  gkyl_array_scale_range(species->L2_f, species->grid.cellVolume, &species->local);

  double L2[1] = {0.0};
  if (app->use_gpu) {
    gkyl_array_reduce_range(species->red_L2_f, species->L2_f, GKYL_SUM, &species->local);
    gkyl_cu_memcpy(L2, species->red_L2_f, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    gkyl_array_reduce_range(L2, species->L2_f, GKYL_SUM, &species->local);
  }
  double L2_global[1] = {0.0};
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, L2, L2_global);

  gkyl_dynvec_append(species->integ_L2_f, tm, L2_global);
}

void
pkpm_species_coll_tm(gkyl_pkpm_app *app)
{
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      struct gkyl_dg_updater_lbo_pkpm_tm tm =
        gkyl_dg_updater_lbo_pkpm_get_tm(app->species[i].lbo.coll_slvr);
      app->stat.species_lbo_coll_diff_tm[i] = tm.diff_tm;
      app->stat.species_lbo_coll_drag_tm[i] = tm.drag_tm;
    }
  }
}

void
pkpm_species_tm(gkyl_pkpm_app *app)
{
  app->stat.species_rhs_tm = 0.0;
  app->stat.fluid_species_rhs_tm = 0.0;
  for (int i = 0; i < app->num_species; ++i) {
    app->stat.species_rhs_tm += app->species[i].kinetic_tm;
    app->stat.fluid_species_rhs_tm += app->species[i].fluid_tm;
  }
}

// release resources for PKPM species
void
pkpm_species_release(const gkyl_pkpm_app *app, const struct pkpm_species *s)
{
  // release various arrays
  gkyl_array_release(s->f);
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->fluid);
  gkyl_array_release(s->fluid1);
  gkyl_array_release(s->fluidnew);
  gkyl_array_release(s->fluid_dup);

  gkyl_array_release(s->cflrate_f);
  gkyl_array_release(s->cflrate_fluid);

  gkyl_wv_eqn_release(s->equation);

  gkyl_array_release(s->bc_buffer_dist);
  gkyl_array_release(s->bc_buffer_lo_fixed_dist);
  gkyl_array_release(s->bc_buffer_up_fixed_dist);
  gkyl_array_release(s->bc_buffer_fluid);
  gkyl_array_release(s->bc_buffer_lo_fixed_fluid);
  gkyl_array_release(s->bc_buffer_up_fixed_fluid);

  gkyl_comm_release(s->comm);

  if (app->use_gpu) {
    gkyl_array_release(s->f_host);
    gkyl_array_release(s->fluid_host);
  }

  gkyl_array_release(s->qmem);

  // release moment data
  pkpm_species_moment_release(app, &s->pkpm_moms);
  pkpm_species_moment_release(app, &s->pkpm_moms_diag);

  // release pkpm arrays data and updaters
  gkyl_array_release(s->pkpm_u);
  gkyl_array_release(s->pkpm_grad_u);
  gkyl_array_release(s->pkpm_theta);
  gkyl_array_release(s->pkpm_a_P);
  gkyl_array_release(s->pkpm_rhs_moms);
  gkyl_array_release(s->pkpm_p_ij);
  gkyl_array_release(s->pkpm_nu_sum);
  gkyl_array_release(s->pkpm_nu_vthsq_sum);
  gkyl_array_release(s->cell_avg_prim);
  gkyl_array_release(s->conf_flux_surf);
  gkyl_array_release(s->vel_flux_surf);
  gkyl_array_release(s->fluid_flux_surf);
  gkyl_array_release(s->pkpm_roe_mass_flux);
  gkyl_array_release(s->pkpm_kin_mass_flux);
  gkyl_array_release(s->pkpm_kin_mass_flux_face);
  gkyl_array_release(s->integ_pkpm_mom);
  gkyl_array_release(s->fluid_io);
  gkyl_array_release(s->pkpm_vars_io);
  gkyl_dg_calc_pkpm_vars_release(s->calc_pkpm_vars);
  gkyl_dg_calc_pkpm_vars_release(s->calc_pkpm_vars_ext);
  gkyl_dg_pkpm_conf_flux_surf_release(s->calc_conf_flux);
  gkyl_dg_pkpm_vel_flux_surf_release(s->calc_vel_flux);
  gkyl_dg_pkpm_fluid_flux_surf_release(s->calc_fluid_flux);
  gkyl_hyper_dg_release(s->slvr_kinetic_conf);
  gkyl_hyper_dg_release(s->slvr_kinetic_vel);
  gkyl_hyper_dg_release(s->slvr_fluid);
  gkyl_dg_eqn_release(s->eqn_kinetic_conf);
  gkyl_dg_eqn_release(s->eqn_kinetic_vel);
  gkyl_dg_eqn_release(s->eqn_fluid);

  gkyl_array_release(s->L2_f);
  gkyl_dynvec_release(s->integ_L2_f);
  gkyl_dynvec_release(s->integ_diag);

  gkyl_array_release(s->app_accel);
  if (s->has_app_accel) {
    if (app->use_gpu) {
      gkyl_array_release(s->app_accel_host);
    }
    gkyl_proj_on_basis_release(s->app_accel_proj);
  }

  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    pkpm_species_lbo_release(app, &s->lbo);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d = 0; d < app->cdim; ++d) {
    gkyl_bc_basic_release(s->bc_lo_dist[d]);
    gkyl_bc_basic_release(s->bc_up_dist[d]);
    gkyl_bc_basic_release(s->bc_lo_fluid[d]);
    gkyl_bc_basic_release(s->bc_up_fluid[d]);
  }

  if (app->use_gpu) {
    gkyl_cu_free(s->omegaCfl_ptr_dist);
    gkyl_cu_free(s->omegaCfl_ptr_fluid);
    gkyl_cu_free(s->red_L2_f);
    gkyl_cu_free(s->red_integ_diag);
  } else {
    gkyl_free(s->omegaCfl_ptr_dist);
    gkyl_free(s->omegaCfl_ptr_fluid);
  }
}
