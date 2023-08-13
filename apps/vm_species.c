#include "gkyl_dg_bin_ops.h"
#include "gkyl_util.h"
#include <assert.h>

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
  gkyl_create_grid_ranges(&s->grid, ghost, &s->local_ext, &s->local);
  
  skin_ghost_ranges_init(&s->skin_ghost, &s->local_ext, ghost);
  
  // velocity space grid
  gkyl_rect_grid_init(&s->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&s->grid_vel, ghost_vel, &s->local_ext_vel, &s->local_vel);

  // allocate buffer for applying periodic BCs
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<cdim; ++d) {
    long vol = s->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }

  // determine field-type 
  s->model_id = s->info.model_id; 
  s->field_id = app->has_field ? app->field->info.field_id : GKYL_FIELD_NULL;

  // allocate distribution function arrays
  if (s->model_id == GKYL_MODEL_PKPM) {
    // PKPM has two distribution functions F_0 and T_perp*G, first two Laguerre moments
    s->f = mkarr(app->use_gpu, 2*app->basis.num_basis, s->local_ext.volume);
    s->f1 = mkarr(app->use_gpu, 2*app->basis.num_basis, s->local_ext.volume);
    s->fnew = mkarr(app->use_gpu, 2*app->basis.num_basis, s->local_ext.volume);

    s->f_host = s->f;
    if (app->use_gpu)
      s->f_host = mkarr(false, 2*app->basis.num_basis, s->local_ext.volume);

    // Distribution function arrays for coupling different Laguerre moments
    // g_dist_source has two components: 
    // [2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)), 
    // (-vpar div(b) + bb:grad(u) - div(u) - 2 nu) T_perp/m G + 2 nu vth^2 F_0]
    // First component is mirror force source *distribution*, second component is *total* vperp characteristics source.
    s->g_dist_source = mkarr(app->use_gpu, 2*app->basis.num_basis, s->local_ext.volume); 
    s->F_k_p_1 = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume); 
    s->F_k_m_1 = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume); 

    gkyl_array_clear(s->g_dist_source, 0.0);
    gkyl_array_clear(s->F_k_p_1, 0.0);
    gkyl_array_clear(s->F_k_m_1, 0.0);

    s->bc_buffer = mkarr(app->use_gpu, 2*app->basis.num_basis, buff_sz);
    // buffer arrays for fixed boundary conditions
    // JJ 10/25/22 somewhat hacky, there is probably a better way to do this
    s->bc_buffer_lo_fixed = mkarr(app->use_gpu, 2*app->basis.num_basis, buff_sz);
    s->bc_buffer_up_fixed = mkarr(app->use_gpu, 2*app->basis.num_basis, buff_sz);
  }
  else {
    s->f = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
    s->f1 = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
    s->fnew = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);

    s->f_host = s->f;
    if (app->use_gpu)
      s->f_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);

    s->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
    // buffer arrays for fixed boundary conditions
    // JJ 10/25/22 somewhat hacky, there is probably a better way to do this
    s->bc_buffer_lo_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
    s->bc_buffer_up_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  }

  // allocate cflrate (scalar array)
  s->cflrate = mkarr(app->use_gpu, 1, s->local_ext.volume);

  if (app->use_gpu)
    s->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
  else 
    s->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  // allocate array to store q/m*(E,B) or potential depending on equation system
  s->qmem = 0;
  if (s->field_id  == GKYL_FIELD_E_B)
    s->qmem = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_ext.volume);

  // allocate array to store relativistic variables if present
  // Since p/gamma, gamma, gamma_inv are a geometric quantities, can pre-compute it here
  s->p_over_gamma = 0;
  s->gamma = 0;
  s->gamma_inv = 0;
  // V_drift, GammaV2, GammaV_inv are derived quantities from the bulk fluid velocity
  s->V_drift = 0;
  s->GammaV2 = 0;
  s->GammaV_inv = 0;
  // PKPM pointers
  s->m1i_pkpm = 0;
  s->pkpm_div_ppar = 0;
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

    struct gkyl_dg_vlasov_sr_auxfields aux_inp = {.qmem = s->qmem, 
      .p_over_gamma = s->p_over_gamma};
    // create solver
    s->slvr = gkyl_dg_updater_vlasov_new(&s->grid, &app->confBasis, &app->basis, 
      &app->local, &s->local_vel, &s->local, s->model_id, s->field_id, &aux_inp, app->use_gpu);
  }
  else if (s->model_id  == GKYL_MODEL_PKPM) {
    // Get pointer to fluid species object for coupling
    s->pkpm_fluid_species = vm_find_fluid_species(app, s->info.pkpm_fluid_species);
    s->pkpm_fluid_index = vm_find_fluid_species_idx(app, s->info.pkpm_fluid_species);

    // pkpm moments for update (rho, p_par, p_perp)
    vm_species_moment_init(app, s, &s->pkpm_moms, "PKPM");
    // pkpm moments for diagnostics (rho, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp)
    // For simplicity, is_integrated flag also used by PKPM to turn on diagnostics
    vm_species_moment_init(app, s, &s->pkpm_moms_diag, "Integrated");

    // Current density for accumulating onto electric field
    s->m1i_pkpm = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
    // div(p_par b_hat), for self-consistent total pressure force
    s->pkpm_div_ppar = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    // allocate array to store primitive moments : [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp]
    // pressure p_ij : (p_par - p_perp) b_i b_j + p_perp g_ij
    s->pkpm_prim = mkarr(app->use_gpu, 6*app->confBasis.num_basis, app->local_ext.volume);
    s->pkpm_p_ij = mkarr(app->use_gpu, 6*app->confBasis.num_basis, app->local_ext.volume);
    // boolean array for if we are only using the cell average for primitive variables
    s->cell_avg_prim = mk_int_arr(app->use_gpu, 1, app->local_ext.volume);

    // Surface primitive variables. Ordered as:
    // [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 3.0*Txx_xl/m, 3.0*Txx_xr/m, 
    //  ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 3.0*Tyy_yl/m, 3.0*Tyy_yr/m, 
    //  ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr, 3.0*Tzz_zl/m, 3.0*Tzz_zr/m] 
    int Ncomp_surf = 2*cdim*3+2*cdim;
    int Nbasis_surf = app->confBasis.num_basis/(app->confBasis.poly_order + 1); // *only valid for tensor bases for cdim > 1*
    s->pkpm_prim_surf = mkarr(app->use_gpu, Ncomp_surf*Nbasis_surf, app->local_ext.volume);
    // Surface expansion of Lax penalization lambda_i = |u_i| + sqrt(3*P_ii/rho)
    s->pkpm_lax = mkarr(app->use_gpu, 2*cdim*Nbasis_surf, app->local_ext.volume);

    // allocate array for pkpm acceleration variables, stored in pkpm_accel: 
    // 0: div_b (divergence of magnetic field unit vector)
    // 1: bb_grad_u (bb : grad(u))
    // 2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
    // 3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - nu + nu rho vth^2/p_perp)
    // 4: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
    s->pkpm_accel = mkarr(app->use_gpu, 5*app->confBasis.num_basis, app->local_ext.volume); 

    // updater for computing pkpm variables 
    // pressure, primitive variables, and acceleration variables
    // also stores kernels for computing source terms, integrated variables
    // Two instances, one over extended range and one over local range for ease of handling boundary conditions
    s->calc_pkpm_vars_ext = gkyl_dg_calc_pkpm_vars_new(&app->grid, &app->confBasis, &app->local_ext, app->use_gpu);
    s->calc_pkpm_vars = gkyl_dg_calc_pkpm_vars_new(&app->grid, &app->confBasis, &app->local, app->use_gpu); 
    // updater for computing pkpm distribution function variables
    // div(p_par b_hat) for self-consistent total pressure force and distribution function sources for
    // Laguerre couplings and vperp characteristics 
    s->calc_pkpm_dist_vars = gkyl_dg_calc_pkpm_dist_vars_new(&s->grid, &app->confBasis, app->use_gpu);

    // allocate arrays for integrated quantities and I/O of fluid variables
    // since these require kinetic species information, kinetic species handles integrations and I/O
    // of fluid variables 
    // array for storing integrated fluid variables in each cell
    s->integ_pkpm_mom = mkarr(app->use_gpu, 9, app->local_ext.volume);
    // arrays for I/O, fluid_io 
    s->fluid_io = mkarr(app->use_gpu, 10*app->confBasis.num_basis, app->local_ext.volume);
    s->pkpm_vars_io = mkarr(app->use_gpu, 10*app->confBasis.num_basis, app->local_ext.volume);
    s->fluid_io_host = s->fluid_io;
    s->pkpm_vars_io_host = s->pkpm_vars_io;
    if (app->use_gpu) {
      s->fluid_io_host = mkarr(false, 10*app->confBasis.num_basis, app->local_ext.volume);
      s->pkpm_vars_io_host = mkarr(false, 10*app->confBasis.num_basis, app->local_ext.volume);
    }
    struct gkyl_dg_vlasov_pkpm_auxfields aux_inp = {.bvar = app->field->bvar, 
      .pkpm_prim = s->pkpm_prim, .pkpm_prim_surf = s->pkpm_prim_surf, 
      .pkpm_accel_vars = s->pkpm_accel, .pkpm_lax = s->pkpm_lax, 
      .g_dist_source = s->g_dist_source};
    // create solver
    s->slvr = gkyl_dg_updater_vlasov_new(&s->grid, &app->confBasis, &app->basis, 
      &app->local, &s->local_vel, &s->local, s->model_id, s->field_id, &aux_inp, app->use_gpu);
  }
  else {
    struct gkyl_dg_vlasov_auxfields aux_inp = {.field = s->qmem, 
      .ext_field = 0, .cot_vec = 0, .alpha_geo = 0};
    // create solver
    s->slvr = gkyl_dg_updater_vlasov_new(&s->grid, &app->confBasis, &app->basis, 
      &app->local, &s->local_vel, &s->local, s->model_id, s->field_id, &aux_inp, app->use_gpu);
  }

  // acquire equation object
  s->eqn_vlasov = gkyl_dg_updater_vlasov_acquire_eqn(s->slvr);
  
  // allocate data for momentum (for use in current accumulation)
  vm_species_moment_init(app, s, &s->m1i, "M1i");
  // allocate date for density (for use in charge density accumulation)
  vm_species_moment_init(app, s, &s->m0, "M0");

  int ndm = s->info.num_diag_moments;
  // allocate data for diagnostic moments
  s->moms = gkyl_malloc(sizeof(struct vm_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    vm_species_moment_init(app, s, &s->moms[m], s->info.diag_moments[m]);

  // array for storing f^2 in each cell
  s->L2_f = mkarr(app->use_gpu, 1, s->local_ext.volume);
  // allocate data for integrated moments
  vm_species_moment_init(app, s, &s->integ_moms, "Integrated");
  if (app->use_gpu) {
    s->red_L2_f = gkyl_cu_malloc(sizeof(double));
    // GKYL_MAX_DIM = 7, allocates a 9 component double
    s->red_integ_diag = gkyl_cu_malloc(sizeof(double[2+GKYL_MAX_DIM]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments and f^2
  s->integ_L2_f = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  if (s->model_id  == GKYL_MODEL_PKPM)
    s->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, 9);
  else
    s->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  s->is_first_integ_L2_write_call = true;
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

  // set species source id
  s->source_id = s->info.source.source_id;
  
  // determine collision type to use in vlasov update
  s->collision_id = s->info.collisions.collision_id;
  if (s->collision_id == GKYL_LBO_COLLISIONS) {
    vm_species_lbo_init(app, s, &s->lbo);
  }

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int dir=0; dir<app->cdim; ++dir) {
    s->lower_bc[dir] = s->upper_bc[dir] = GKYL_SPECIES_COPY;
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
  // Certain operations fail if absorbing BCs used because absorbing BCs 
  // means the mass density is 0 in the ghost cells (divide by zero)
  s->bc_is_absorb = false;
  for (int d=0; d<app->cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (s->lower_bc[d] == GKYL_SPECIES_COPY) {
      bctype = GKYL_BC_COPY;
    }
    else if (s->lower_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype = GKYL_BC_ABSORB;
      s->bc_is_absorb = true;
    }
    else if (s->lower_bc[d] == GKYL_SPECIES_REFLECT) {
      if (s->model_id  == GKYL_MODEL_PKPM)
        bctype = GKYL_BC_PKPM_SPECIES_REFLECT;
      else
        bctype = GKYL_BC_REFLECT;
    }
    else if (s->lower_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      bctype = GKYL_BC_FIXED_FUNC;
    }

    s->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, &s->local_ext, ghost, bctype,
      app->basis_on_dev.basis, s->f->ncomp, app->cdim, app->use_gpu);
    // Upper BC updater. Copy BCs by default.
    if (s->upper_bc[d] == GKYL_SPECIES_COPY) {
      bctype = GKYL_BC_COPY;
    }
    else if (s->upper_bc[d] == GKYL_SPECIES_ABSORB) {
      bctype = GKYL_BC_ABSORB;
      s->bc_is_absorb = true;
    }
    else if (s->upper_bc[d] == GKYL_SPECIES_REFLECT) {
      if (s->model_id  == GKYL_MODEL_PKPM)
        bctype = GKYL_BC_PKPM_SPECIES_REFLECT;
      else
        bctype = GKYL_BC_REFLECT;
    }
    else if (s->upper_bc[d] == GKYL_SPECIES_FIXED_FUNC) {
      bctype = GKYL_BC_FIXED_FUNC;
    }

    s->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, &s->local_ext, ghost, bctype,
      app->basis_on_dev.basis, s->f->ncomp, app->cdim, app->use_gpu);
  }
}

void
vm_species_apply_ic(gkyl_vlasov_app *app, struct vm_species *species, double t0)
{
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj;
  if (species->model_id == GKYL_MODEL_PKPM)
    // PKPM often run with hybrid basis, so do a higher order projection for ICs
    proj = gkyl_proj_on_basis_new(&species->grid, &app->basis,
      8, 2, species->info.init, species->info.ctx);
  else
    proj = gkyl_proj_on_basis_new(&species->grid, &app->basis,
      poly_order+1, 1, species->info.init, species->info.ctx);

  // run updater; need to project onto extended range for ease of handling
  // subsequent operations over extended range such as primitive variable computations
  // This is needed to fill the corner cells as the corner cells may not be filled by
  // boundary conditions and we cannot divide by 0 anywhere or the weak divisions will fail
  gkyl_proj_on_basis_advance(proj, t0, &species->local_ext, species->f_host);
  gkyl_proj_on_basis_release(proj);    

  if (app->use_gpu) // note: f_host is same as f when not on GPUs
    gkyl_array_copy(species->f, species->f_host);

  // we are pre-computing acceleration for now as it is time-independent
  vm_species_calc_accel(app, species, t0);

  // we are pre-computing source for now as it is time-independent
  vm_species_source_calc(app, species, t0);

  // copy contents of initial conditions into buffer if specific BCs require them
  // *only works in x dimension for now*
  gkyl_bc_basic_buffer_fixed_func(species->bc_lo[0], species->bc_buffer_lo_fixed, species->f);
  gkyl_bc_basic_buffer_fixed_func(species->bc_up[0], species->bc_buffer_up_fixed, species->f);
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
vm_species_calc_pkpm_vars(gkyl_vlasov_app *app, struct vm_species *species, 
  const struct gkyl_array *fin, const struct gkyl_array *fluidin[])
{
  if (species->model_id == GKYL_MODEL_PKPM) {
    gkyl_array_clear(species->pkpm_div_ppar, 0.0);
    gkyl_array_clear(species->pkpm_p_ij, 0.0);
    gkyl_array_clear(species->pkpm_prim, 0.0);
    gkyl_array_clear(species->pkpm_prim_surf, 0.0);

    // Compute rho, p_par, & p_perp
    vm_species_moment_calc(&species->pkpm_moms, species->local_ext,
      app->local_ext, fin);
    // Compute div(p_par b_hat) for consistent pressure force in vpar acceleration
    gkyl_dg_calc_pkpm_dist_vars_div_ppar(species->calc_pkpm_dist_vars, 
        &app->local, &species->local, app->field->bvar, fin, species->pkpm_div_ppar);
    gkyl_array_scale(species->pkpm_div_ppar, species->info.mass);
    // Compute p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij
    gkyl_dg_calc_pkpm_vars_pressure(species->calc_pkpm_vars, &app->local_ext, 
      app->field->bvar, species->pkpm_moms.marr, species->pkpm_p_ij);
    // Compute primitive variables in both the volume and on surfaces
    if (species->bc_is_absorb) {
      gkyl_dg_calc_pkpm_vars_advance(species->calc_pkpm_vars,
        species->pkpm_moms.marr, fluidin[species->pkpm_fluid_index], 
        species->pkpm_div_ppar, species->cell_avg_prim, 
        species->pkpm_prim); 
      gkyl_dg_calc_pkpm_vars_surf_advance(species->calc_pkpm_vars, 
        species->pkpm_moms.marr, fluidin[species->pkpm_fluid_index], 
        species->pkpm_p_ij, species->cell_avg_prim, 
        species->pkpm_prim_surf);
    }
    else {
      gkyl_dg_calc_pkpm_vars_advance(species->calc_pkpm_vars_ext,
        species->pkpm_moms.marr, fluidin[species->pkpm_fluid_index], 
        species->pkpm_div_ppar, species->cell_avg_prim, 
        species->pkpm_prim); 
      gkyl_dg_calc_pkpm_vars_surf_advance(species->calc_pkpm_vars_ext, 
        species->pkpm_moms.marr, fluidin[species->pkpm_fluid_index], 
        species->pkpm_p_ij, species->cell_avg_prim, 
        species->pkpm_prim_surf);
    }
  }
}

void
vm_species_calc_pkpm_update_vars(gkyl_vlasov_app *app, struct vm_species *species, 
  const struct gkyl_array *fin)
{
  if (species->model_id == GKYL_MODEL_PKPM) {
    gkyl_array_clear(species->pkpm_accel, 0.0);
    gkyl_array_clear(species->pkpm_lax, 0.0);
    gkyl_dg_calc_pkpm_vars_accel(species->calc_pkpm_vars, &app->local, 
      app->field->bvar, species->pkpm_prim_surf, 
      species->pkpm_prim, species->lbo.nu_sum, 
      species->pkpm_lax, species->pkpm_accel); 

    // Calculate distrbution functions for coupling different Laguerre moments
    gkyl_dg_calc_pkpm_dist_vars_mirror_force(species->calc_pkpm_dist_vars, 
      &app->local, &species->local, 
      species->pkpm_prim, species->lbo.nu_prim_moms, species->pkpm_accel, 
      fin, species->F_k_p_1, 
      species->g_dist_source, species->F_k_m_1);
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

  if (species->collision_id == GKYL_LBO_COLLISIONS)
    vm_species_lbo_rhs(app, species, &species->lbo, fin, rhs);
  
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

      switch (species->upper_bc[d]) {
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
}

void
vm_species_calc_L2(gkyl_vlasov_app *app, double tm, const struct vm_species *species)
{
  gkyl_dg_calc_l2_range(app->basis, 0, species->L2_f, 0, species->f, species->local);
  gkyl_array_scale_range(species->L2_f, species->grid.cellVolume, species->local);
  
  double L2[1] = { 0.0 };
  if (app->use_gpu) {
    gkyl_array_reduce_range(species->red_L2_f, species->L2_f, GKYL_SUM, species->local);
    gkyl_cu_memcpy(L2, species->red_L2_f, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else { 
    gkyl_array_reduce_range(L2, species->L2_f, GKYL_SUM, species->local);
  }
  
  gkyl_dynvec_append(species->integ_L2_f, tm, L2);
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

  if (app->use_gpu)
    gkyl_array_release(s->f_host);
  
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

  gkyl_dg_eqn_release(s->eqn_vlasov);
  gkyl_dg_updater_vlasov_release(s->slvr);

  // Release arrays for different types of Vlasov equations
  if (s->field_id  == GKYL_FIELD_E_B)
    gkyl_array_release(s->qmem);

  if (s->model_id  == GKYL_MODEL_SR) {
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
  if (s->model_id == GKYL_MODEL_PKPM) {
    vm_species_moment_release(app, &s->pkpm_moms);
    vm_species_moment_release(app, &s->pkpm_moms_diag);

    gkyl_array_release(s->g_dist_source);
    gkyl_array_release(s->F_k_m_1);
    gkyl_array_release(s->F_k_p_1);
    gkyl_array_release(s->m1i_pkpm);
    gkyl_array_release(s->pkpm_div_ppar);
    gkyl_array_release(s->pkpm_prim);
    gkyl_array_release(s->pkpm_prim_surf);
    gkyl_array_release(s->pkpm_p_ij);
    gkyl_array_release(s->pkpm_lax);
    gkyl_array_release(s->cell_avg_prim);
    gkyl_array_release(s->pkpm_accel);
    gkyl_dg_calc_pkpm_vars_release(s->calc_pkpm_vars);
    gkyl_dg_calc_pkpm_vars_release(s->calc_pkpm_vars_ext);
    gkyl_dg_calc_pkpm_dist_vars_release(s->calc_pkpm_dist_vars);
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

  if (s->collision_id == GKYL_LBO_COLLISIONS)
    vm_species_lbo_release(app, &s->lbo);

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
  else
    gkyl_free(s->omegaCfl_ptr);
}
