#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_vlasov_priv.h>

#include <assert.h>
#include <time.h>

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

  // determine field-type 
  s->model_id = s->info.model_id; 
  s->field_id = app->has_field ? app->field->info.field_id : GKYL_FIELD_NULL;

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

  // allocate array to store q/m*(E,B) or potentials (phi, A) depending on equation system
  if (s->field_id == GKYL_FIELD_E_B)
    s->qmem = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);
  else
    s->qmem = mkarr(app->use_gpu, 4*app->confBasis.num_basis, app->local_ext.volume);

  if (s->model_id  == GKYL_MODEL_SR) {
    // Allocate special relativistic variables, p/gamma, gamma, & 1/gamma
    s->p_over_gamma = mkarr(app->use_gpu, vdim*app->velBasis.num_basis, s->local_vel.volume);
    s->p_over_gamma_host = s->p_over_gamma;
    s->gamma = mkarr(app->use_gpu, app->velBasis.num_basis, s->local_vel.volume);
    s->gamma_host = s->gamma;
    s->gamma_inv = mkarr(app->use_gpu, app->velBasis.num_basis, s->local_vel.volume);
    s->gamma_inv_host = s->gamma_inv;
    // Projection routines are done on CPU, need to allocate host arrays if simulation on GPU
    if (app->use_gpu) {
      s->p_over_gamma_host = mkarr(false, vdim*app->velBasis.num_basis, s->local_vel.volume);
      s->gamma_host = mkarr(false, app->velBasis.num_basis, s->local_vel.volume);
      s->gamma_inv_host = mkarr(false, app->velBasis.num_basis, s->local_vel.volume);
    }
    // Project p/gamma, gamma, & 1/gamma
    gkyl_calc_sr_vars_init_p_vars(&s->grid_vel, &app->velBasis, &s->local_vel, 
      s->p_over_gamma_host, s->gamma_host, s->gamma_inv_host);
    // Copy CPU arrays to GPU 
    if (app->use_gpu) {
      gkyl_array_copy(s->p_over_gamma, s->p_over_gamma_host);
      gkyl_array_copy(s->gamma, s->gamma_host);
      gkyl_array_copy(s->gamma_inv, s->gamma_inv_host);
    }
    // Derived quantities array
    s->V_drift = mkarr(app->use_gpu, vdim*app->confBasis.num_basis, app->local_ext.volume);
    s->GammaV2 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    s->GammaV_inv = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    if (app->use_gpu)
      s->V_drift_mem = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
    else
      s->V_drift_mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);

    // by default, we do not have zero-flux boundary conditions in any direction
    bool is_zero_flux[GKYL_MAX_DIM] = {false};

    struct gkyl_dg_vlasov_sr_auxfields aux_inp = {.qmem = s->qmem, .p_over_gamma = s->p_over_gamma};
    // create solver
    s->slvr = gkyl_dg_updater_vlasov_new(&s->grid, &app->confBasis, &app->basis, 
      &app->local, &s->local_vel, &s->local, is_zero_flux, s->model_id, s->field_id, &aux_inp, app->use_gpu);
  }
  else {
    // by default, we do not have zero-flux boundary conditions in any direction
    bool is_zero_flux[GKYL_MAX_DIM] = {false};

    struct gkyl_dg_vlasov_auxfields aux_inp = {.field = s->qmem, .cot_vec = 0, 
      .alpha_surf = 0, .sgn_alpha_surf = 0, .const_sgn_alpha = 0 };
    // create solver
    s->slvr = gkyl_dg_updater_vlasov_new(&s->grid, &app->confBasis, &app->basis, 
      &app->local, &s->local_vel, &s->local, is_zero_flux, s->model_id, s->field_id, &aux_inp, app->use_gpu);
  }

  // acquire equation object
  s->eqn_vlasov = gkyl_dg_updater_vlasov_acquire_eqn(s->slvr);

  // allocate data for momentum (for use in current accumulation)
  vm_species_moment_init(app, s, &s->m1i, "M1i");
  // allocate date for density (for use in charge density accumulation and weak division for V_drift)
  vm_species_moment_init(app, s, &s->m0, "M0");
  // allocate data for integrated moments
  vm_species_moment_init(app, s, &s->integ_moms, "Integrated");

  // allocate data for diagnostic moments
  int ndm = s->info.num_diag_moments;
  s->moms = gkyl_malloc(sizeof(struct vm_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    vm_species_moment_init(app, s, &s->moms[m], s->info.diag_moments[m]);

  // array for storing f^2 in each cell
  s->L2_f = mkarr(app->use_gpu, 1, s->local_ext.volume);
  if (app->use_gpu) {
    s->red_L2_f = gkyl_cu_malloc(sizeof(double));
    s->red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments and f^2
  s->integ_L2_f = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  s->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  s->is_first_integ_L2_write_call = true;
  s->is_first_integ_write_call = true;

  s->has_accel = false;
  // setup applied acceleration
  if (s->info.accel) {
    s->has_accel = true;
    if (s->info.accel_evolve)
      s->accel_evolve = s->info.accel_evolve;
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

  // initialize projection routine for initial conditions
  vm_species_projection_init(app, s, s->info.projection, &s->proj_init);

  // set species source id
  s->source_id = s->info.source.source_id;
  
  // determine collision type to use in vlasov update
  s->collision_id = s->info.collisions.collision_id;
  s->lbo = (struct vm_lbo_collisions) { };
  s->bgk = (struct vm_bgk_collisions) { };
  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    vm_species_lbo_init(app, s, &s->lbo);
  }
  else if (s->collision_id == GKYL_BGK_COLLISIONS) {
    vm_species_bgk_init(app, s, &s->bgk);
  }

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int dir=0; dir<app->cdim; ++dir) {
    s->lower_bc[dir].type = s->upper_bc[dir].type = GKYL_SPECIES_COPY;
    if (is_np[dir]) {
      const struct gkyl_vlasov_bcs *bc;
      if (dir == 0)
        bc = &s->info.bcx;
      else if (dir == 1)
        bc = &s->info.bcy;
      else
        bc = &s->info.bcz;

      s->lower_bc[dir] = bc->lower;
      s->upper_bc[dir] = bc->upper;
    }
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->lower_skin[dir], &s->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &s->local_ext, ghost);
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->upper_skin[dir], &s->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &s->local_ext, ghost);
  }

  // intitalize boundary flux updater, needs to be done after skin and ghost ranges
  vm_species_bflux_init(app, s, &s->bflux);

  // allocate buffer for applying periodic BCs
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<cdim; ++d) {
    long vol = GKYL_MAX2(s->lower_skin[d].volume, s->upper_skin[d].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  s->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  // buffer arrays for fixed function boundary conditions on distribution function
  s->bc_buffer_lo_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  s->bc_buffer_up_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

  for (int d=0; d<cdim; ++d) {

    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (s->lower_bc[d].type == GKYL_SPECIES_EMISSION) {
      s->emit_lo = true;
      vm_species_emission_init(app, &s->bc_emission_lo, d, GKYL_LOWER_EDGE, s->lower_bc[d].aux_ctx,
        app->use_gpu);
    }
    else {
      if (s->lower_bc[d].type == GKYL_SPECIES_COPY)
        bctype = GKYL_BC_COPY;
      else if (s->lower_bc[d].type == GKYL_SPECIES_ABSORB)
        bctype = GKYL_BC_ABSORB;
      else if (s->lower_bc[d].type == GKYL_SPECIES_REFLECT)
        bctype = GKYL_BC_REFLECT;
      else if (s->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC)
        bctype = GKYL_BC_FIXED_FUNC;

      s->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev.basis,
        &s->lower_skin[d], &s->lower_ghost[d], s->f->ncomp, app->cdim, app->use_gpu);
    }

    // Upper BC updater. Copy BCs by default.
    if (s->upper_bc[d].type == GKYL_SPECIES_EMISSION) {
      s->emit_up = true;
      vm_species_emission_init(app, &s->bc_emission_up, d, GKYL_UPPER_EDGE, s->upper_bc[d].aux_ctx,
        app->use_gpu);
    }
    else {
      if (s->upper_bc[d].type == GKYL_SPECIES_COPY)
        bctype = GKYL_BC_COPY;
      else if (s->upper_bc[d].type == GKYL_SPECIES_ABSORB)
        bctype = GKYL_BC_ABSORB;
      else if (s->upper_bc[d].type == GKYL_SPECIES_REFLECT)
        bctype = GKYL_BC_REFLECT;
      else if (s->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC)
        bctype = GKYL_BC_FIXED_FUNC;

      s->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.basis,
        &s->upper_skin[d], &s->upper_ghost[d], s->f->ncomp, app->cdim, app->use_gpu);
    }
  }
}

void
vm_species_apply_ic(gkyl_vlasov_app *app, struct vm_species *species, double t0)
{
  vm_species_projection_calc(app, species, &species->proj_init, species->f, t0);

  // Pre-compute applied acceleration in case it's time-independent
  vm_species_calc_accel(app, species, t0);

  // we are pre-computing source for now as it is time-independent
  vm_species_source_calc(app, species, &species->src, t0);

  vm_species_bflux_rhs(app, species, &species->bflux, species->f, species->f1);
  
  // copy contents of initial conditions into buffer if specific BCs require them
  // *only works in x dimension for now*
  /* gkyl_bc_basic_buffer_fixed_func(species->bc_lo[0], species->bc_buffer_lo_fixed, species->f); */
  /* printf("4\n"); */
  /* gkyl_bc_basic_buffer_fixed_func(species->bc_up[0], species->bc_buffer_up_fixed, species->f); */
  /* printf("5\n"); */
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

// Compute the RHS for species update, returning maximum stable
// time-step.
double
vm_species_rhs(gkyl_vlasov_app *app, struct vm_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *em, struct gkyl_array *rhs)
{
  if (species->field_id  == GKYL_FIELD_E_B) {
    double qbym = species->info.charge/species->info.mass;
    gkyl_array_set(species->qmem, qbym, em);

    // Accumulate applied acceleration and/or q/m*(external electromagnetic)
    // fields onto qmem to get the total acceleration
    if (species->has_accel)
      gkyl_array_accumulate(species->qmem, 1.0, species->accel);
    if (app->field->has_ext_em)
      gkyl_array_accumulate(species->qmem, qbym, app->field->ext_em);
  }

  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  gkyl_dg_updater_vlasov_advance(species->slvr, &species->local, 
    fin, species->cflrate, rhs);

  if (species->collision_id == GKYL_LBO_COLLISIONS) {
    vm_species_lbo_rhs(app, species, &species->lbo, fin, rhs);
  }
  else if (species->collision_id == GKYL_BGK_COLLISIONS) {
    vm_species_bgk_rhs(app, species, &species->bgk, fin, rhs);
  }

  vm_species_bflux_rhs(app, species, &species->bflux, fin, rhs);
  
  app->stat.nspecies_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omegaCfl_ptr, species->cflrate, GKYL_MAX, &species->local);

  double omegaCfl_ho[1];
  if (app->use_gpu)
    gkyl_cu_memcpy(omegaCfl_ho, species->omegaCfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    omegaCfl_ho[0] = species->omegaCfl_ptr[0];
  double omegaCfl = omegaCfl_ho[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  
  return app->cfl/omegaCfl;
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for distribution function
void
vm_species_apply_bc(gkyl_vlasov_app *app, const struct vm_species *species, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(species->comm, &species->local, &species->local_ext,
    num_periodic_dir, app->periodic_dirs, f); 
  
  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is periodic
  for (int d=0; d<num_periodic_dir; ++d)
    is_np_bc[app->periodic_dirs[d]] = 0;

  for (int d=0; d<cdim; ++d) {
    if (is_np_bc[d]) {

      switch (species->lower_bc[d].type) {
        case GKYL_SPECIES_EMISSION:
          vm_species_emission_apply_bc(&species->bc_emission_lo, f);
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer_lo_fixed, f);
          break;
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        default:
          break;
      }

      switch (species->upper_bc[d].type) {
        case GKYL_SPECIES_EMISSION:
          vm_species_emission_apply_bc(&species->bc_emission_up, f);
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer_up_fixed, f);
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


void
vm_species_calc_L2(gkyl_vlasov_app *app, double tm, const struct vm_species *species)
{
  gkyl_dg_calc_l2_range(app->basis, 0, species->L2_f, 0, species->f, species->local);
  gkyl_array_scale_range(species->L2_f, species->grid.cellVolume, &species->local);
  
  double L2[1] = { 0.0 };
  if (app->use_gpu) {
    gkyl_array_reduce_range(species->red_L2_f, species->L2_f, GKYL_SUM, &species->local);
    gkyl_cu_memcpy(L2, species->red_L2_f, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    gkyl_array_reduce_range(L2, species->L2_f, GKYL_SUM, &species->local);
  }
  double L2_global[1] = { 0.0 };
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, L2, L2_global);
  
  gkyl_dynvec_append(species->integ_L2_f, tm, L2_global);  
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
  gkyl_array_release(s->bc_buffer_lo_fixed);
  gkyl_array_release(s->bc_buffer_up_fixed);

  vm_species_projection_release(app, &s->proj_init);

  vm_species_bflux_release(app, &s->bflux);

  gkyl_comm_release(s->comm);

  if (app->use_gpu)
    gkyl_array_release(s->f_host);

  gkyl_array_release(s->qmem);

  // Release arrays for different types of Vlasov equations
  if (s->model_id  == GKYL_MODEL_SR) {
    // release relativistic arrays data
    gkyl_array_release(s->p_over_gamma);
    gkyl_array_release(s->gamma);
    gkyl_array_release(s->gamma_inv);
    gkyl_array_release(s->V_drift);
    gkyl_array_release(s->GammaV2);
    gkyl_array_release(s->GammaV_inv);
    if (app->use_gpu) {
      gkyl_array_release(s->p_over_gamma_host);
      gkyl_array_release(s->gamma_host);
      gkyl_array_release(s->gamma_inv_host);
    }
    gkyl_dg_bin_op_mem_release(s->V_drift_mem);
  }

  // release equation object and solver
  gkyl_dg_eqn_release(s->eqn_vlasov);
  gkyl_dg_updater_vlasov_release(s->slvr);

  // release moment data
  vm_species_moment_release(app, &s->m1i);
  vm_species_moment_release(app, &s->m0);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    vm_species_moment_release(app, &s->moms[i]);
  gkyl_free(s->moms);
  vm_species_moment_release(app, &s->integ_moms); 

  gkyl_array_release(s->L2_f);
  gkyl_dynvec_release(s->integ_L2_f);
  gkyl_dynvec_release(s->integ_diag);
  
  if (s->has_accel) {
    gkyl_array_release(s->accel);
    if (app->use_gpu)
      gkyl_array_release(s->accel_host);

    gkyl_proj_on_basis_release(s->accel_proj);
  }

  if (s->source_id) {
    vm_species_source_release(app, &s->src);
  }

  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    vm_species_lbo_release(app, &s->lbo);
  }
  else if (s->collision_id == GKYL_BGK_COLLISIONS) {
    vm_species_bgk_release(app, &s->bgk);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    gkyl_bc_basic_release(s->bc_lo[d]);
    gkyl_bc_basic_release(s->bc_up[d]);
  }
  
  if (app->use_gpu) {
    gkyl_cu_free(s->omegaCfl_ptr);
    gkyl_cu_free(s->red_L2_f);
    gkyl_cu_free(s->red_integ_diag);
  }
  else {
    gkyl_free(s->omegaCfl_ptr);
  }
}
