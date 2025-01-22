#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_proj_on_basis.h>

#include <assert.h>
#include <time.h>

// initialize species grid and arrays
void
gk_neut_species_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s)
{
  int cdim = app->cdim, vdim = app->vdim+1; // neutral species are 3v
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = gk->cells[d];
    lower[d] = gk->lower[d];
    upper[d] = gk->upper[d];
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

  // Velocity space mapping.
  assert(s->info.mapc2p.mapping == 0); // mapped v-space not implemented for neutrals yet.
  s->vel_map = gkyl_velocity_map_new(s->info.mapc2p, s->grid, s->grid_vel,
    s->local, s->local_ext, s->local_vel, s->local_ext_vel, app->use_gpu);
  // allocate distribution function array for initialization and I/O
  s->f = mkarr(app->use_gpu, app->neut_basis.num_basis, s->local_ext.volume);

  s->f_host = s->f;
  if (app->use_gpu)
    s->f_host = mkarr(false, app->neut_basis.num_basis, s->local_ext.volume);

  if (app->use_gpu)
    s->omega_cfl = gkyl_cu_malloc(sizeof(double));
  else 
    s->omega_cfl = gkyl_malloc(sizeof(double));

  // Create skin/ghost ranges fir applying BCs. Only used for dynamic neutrals but included here to avoid
  // code duplication since the "ghost" array is needed.
  for (int dir=0; dir<cdim; ++dir) {
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->lower_skin[dir], &s->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &s->local_ext, ghost);
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->upper_skin[dir], &s->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &s->local_ext, ghost);
  }

  // Global skin and ghost ranges, only valid (i.e. volume>0) in ranges
  // abutting boundaries.
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&s->global_lower_skin[dir], &s->global_lower_ghost[dir],
      dir, GKYL_LOWER_EDGE, &s->global_ext, ghost); 
    gkyl_sub_range_intersect(&s->global_lower_skin[dir], &s->local_ext, &s->global_lower_skin[dir]);
    gkyl_sub_range_intersect(&s->global_lower_ghost[dir], &s->local_ext, &s->global_lower_ghost[dir]);
    gkyl_skin_ghost_ranges(&s->global_upper_skin[dir], &s->global_upper_ghost[dir],
      dir, GKYL_UPPER_EDGE, &s->local_ext, ghost);
    gkyl_sub_range_intersect(&s->global_upper_skin[dir], &s->local_ext, &s->global_upper_skin[dir]);
    gkyl_sub_range_intersect(&s->global_upper_ghost[dir], &s->local_ext, &s->global_upper_ghost[dir]);
  }

  // allocate data for density 
  gk_neut_species_moment_init(app, s, &s->m0, "M0");
  // allocate data for integrated moments
  gk_neut_species_moment_init(app, s, &s->integ_moms, "Integrated");

  // allocate data for diagnostic moments
  int ndm = s->info.num_diag_moments;
  s->moms = gkyl_malloc(sizeof(struct gk_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    gk_neut_species_moment_init(app, s, &s->moms[m], s->info.diag_moments[m]);

  if (app->use_gpu) {
    s->red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
    s->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[vdim+2]));
  } else {
    s->red_integ_diag = gkyl_malloc(sizeof(double[vdim+2]));
    s->red_integ_diag_global = gkyl_malloc(sizeof(double[vdim+2]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments 
  s->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  s->is_first_integ_write_call = true;

  // initialize projection routine for initial conditions
  gk_neut_species_projection_init(app, s, s->info.projection, &s->proj_init);
}

// initialize species object
void
gk_neut_species_new(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s)
{
  int cdim = app->cdim, vdim = app->vdim+1; // neutral species are 3v
  int pdim = cdim+vdim;
  gk_neut_species_init(gk, app, s);
  
  // allocate additional distribution function arrays for time stepping
  s->f1 = mkarr(app->use_gpu, app->neut_basis.num_basis, s->local_ext.volume);
  s->fnew = mkarr(app->use_gpu, app->neut_basis.num_basis, s->local_ext.volume);
  
  // Allocate cflrate (scalar array).
  s->cflrate = mkarr(app->use_gpu, 1, s->local_ext.volume);
  
  // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
  struct gkyl_basis surf_basis, surf_quad_basis;
  gkyl_cart_modal_serendip(&surf_basis, pdim-1, app->poly_order);
  gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, app->poly_order);
  
  // always 3v
  int alpha_surf_sz = 3*surf_basis.num_basis; 
  int sgn_alpha_surf_sz = 3*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points
  
  // allocate arrays to store fields: 
  // 1. alpha_surf (surface phase space flux)
  // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
  // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
  s->alpha_surf = mkarr(app->use_gpu, alpha_surf_sz, s->local_ext.volume);
  s->sgn_alpha_surf = mkarr(app->use_gpu, sgn_alpha_surf_sz, s->local_ext.volume);
  s->const_sgn_alpha = mk_int_arr(app->use_gpu, 3, s->local_ext.volume);
  // 4. cotangent vectors e^i = g^ij e_j 
  s->cot_vec = mkarr(app->use_gpu, 9*app->confBasis.num_basis, app->local_ext.volume);
  
  // Pre-compute alpha_surf, sgn_alpha_surf, const_sgn_alpha, and cot_vec since they are time-independent
  struct gkyl_dg_calc_vlasov_gen_geo_vars *calc_vars = gkyl_dg_calc_vlasov_gen_geo_vars_new(&s->grid, 
    &app->confBasis, &app->neut_basis, app->gk_geom, app->use_gpu);
  gkyl_dg_calc_vlasov_gen_geo_vars_alpha_surf(calc_vars, &app->local, &s->local, &s->local_ext, 
    s->alpha_surf, s->sgn_alpha_surf, s->const_sgn_alpha);
  gkyl_dg_calc_vlasov_gen_geo_vars_cot_vec(calc_vars, &app->local, s->cot_vec);

  // by default, we do not have zero-flux boundary conditions in any direction
  bool is_zero_flux[2*GKYL_MAX_DIM] = {false};

  // Determine which directions are not periodic. If any BCs are zero-flux,
  // need to set it in is_zero_flux.
  // Keep a copy of num_periodic_dir and periodic_dirs in species so we can
  // modify it in GK_IWL BCs without modifying the app's.
  s->num_periodic_dir = app->num_periodic_dir;
  for (int d=0; d<s->num_periodic_dir; ++d)
    s->periodic_dirs[d] = app->periodic_dirs[d];

  for (int d=0; d<app->cdim; ++d) s->bc_is_np[d] = true;
  for (int d=0; d<s->num_periodic_dir; ++d)
    s->bc_is_np[s->periodic_dirs[d]] = false;

  for (int dir=0; dir<app->cdim; ++dir) {
    s->lower_bc[dir].type = s->upper_bc[dir].type = GKYL_SPECIES_COPY;
    if (s->bc_is_np[dir]) {
      const struct gkyl_gyrokinetic_bcs *bc;
      if (dir == 0)
        bc = &s->info.bcx;
      else if (dir == 1)
        bc = &s->info.bcy;
      else
        bc = &s->info.bcz;

      s->lower_bc[dir] = bc->lower;
      s->upper_bc[dir] = bc->upper;
      if (s->lower_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        is_zero_flux[dir] = true;
      }
      if (s->upper_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        is_zero_flux[dir+pdim] = true;
      }
    }
  }

  s->model_id = GKYL_MODEL_GEN_GEO;
  s->react_neut = (struct gk_react) { };
  struct gkyl_dg_vlasov_auxfields aux_inp = {.field = 0, .cot_vec = s->cot_vec, 
    .alpha_surf = s->alpha_surf, .sgn_alpha_surf = s->sgn_alpha_surf, .const_sgn_alpha = s->const_sgn_alpha };
  // Set field type and model id for neutral species in GK system and create solver
  s->field_id = GKYL_FIELD_NULL;
  s->slvr = gkyl_dg_updater_vlasov_new(&s->grid, &app->confBasis, &app->neut_basis, 
    &app->local, &s->local_vel, &s->local, is_zero_flux, s->model_id, s->field_id, &aux_inp, app->use_gpu);

  // acquire equation object
  s->eqn_vlasov = gkyl_dg_updater_vlasov_acquire_eqn(s->slvr);
  
  if (s->info.react_neut.num_react) {
    gk_neut_species_react_init(app, s, s->info.react_neut, &s->react_neut);
  }

  // set species source id
  s->src = (struct gk_source) { };

  // Set cfl and bc function pointers
  s->cfl_func = gk_neut_species_omega_cfl;
  s->bc_func = gk_neut_species_bc_func;
  s->release_func = gk_neut_species_release_func;

  // Allocate buffer needed for BCs.
  long buff_sz = 0;
  for (int dir=0; dir<cdim; ++dir) {
    long vol = GKYL_MAX2(s->lower_skin[dir].volume, s->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  s->bc_buffer = mkarr(app->use_gpu, app->neut_basis.num_basis, buff_sz);
  // buffer arrays for fixed function boundary conditions on distribution function
  s->bc_buffer_lo_fixed = mkarr(app->use_gpu, app->neut_basis.num_basis, buff_sz);
  s->bc_buffer_up_fixed = mkarr(app->use_gpu, app->neut_basis.num_basis, buff_sz);

  for (int d=0; d<cdim; ++d) {
    // Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (s->lower_bc[d].type == GKYL_SPECIES_RECYCLE) {
      ;
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

      s->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev.neut_basis,
        &s->lower_skin[d], &s->lower_ghost[d], s->f->ncomp, app->cdim, app->use_gpu);

      if (s->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
        // Fill the buffer used for BCs.
        struct gk_proj gk_proj_bc_lo;
        gk_neut_species_projection_init(app, s, s->lower_bc[d].projection, &gk_proj_bc_lo);
        gk_neut_species_projection_calc(app, s, &gk_proj_bc_lo, s->f1, 0.0); // Temporarily use f1.
        gkyl_bc_basic_buffer_fixed_func(s->bc_lo[d], s->bc_buffer_lo_fixed, s->f1);
        gkyl_array_clear(s->f1, 0.0);
        gk_neut_species_projection_release(app, &gk_proj_bc_lo);
      }
    }

    if (s->upper_bc[d].type == GKYL_SPECIES_RECYCLE) {
      ;
    }
    else {
      // Upper BC updater. Copy BCs by default.
      if (s->upper_bc[d].type == GKYL_SPECIES_COPY) 
        bctype = GKYL_BC_COPY;
      else if (s->upper_bc[d].type == GKYL_SPECIES_ABSORB) 
        bctype = GKYL_BC_ABSORB;
      else if (s->upper_bc[d].type == GKYL_SPECIES_REFLECT) 
        bctype = GKYL_BC_REFLECT;
      else if (s->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC) 
        bctype = GKYL_BC_FIXED_FUNC;

      s->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.neut_basis,
        &s->upper_skin[d], &s->upper_ghost[d], s->f->ncomp, app->cdim, app->use_gpu);

      if (s->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC) {
        // Fill the buffer used for BCs.
        struct gk_proj gk_proj_bc_up;
        gk_neut_species_projection_init(app, s, s->upper_bc[d].projection, &gk_proj_bc_up);
        gk_neut_species_projection_calc(app, s, &gk_proj_bc_up, s->f1, 0.0); // Temporarily use f1.
        gkyl_bc_basic_buffer_fixed_func(s->bc_up[d], s->bc_buffer_up_fixed, s->f1);
        gkyl_array_clear(s->f1, 0.0);
        gk_neut_species_projection_release(app, &gk_proj_bc_up);
      }
    }
  }
}

void
gk_neut_species_static_new(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s)
{
  gk_neut_species_init(gk, app, s);

  // set pointers for RK methods
  s->f1 = s->f;
  s->fnew = s->f;

  // Set rhs function
  s->cfl_func = gk_neut_species_static_omega_cfl;
  s->bc_func = gk_neut_species_static_bc_func;
  s->release_func = gk_neut_species_static_release_func;
}

void
gk_neut_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_neut_species *species, double t0)
{
  gk_neut_species_projection_calc(app, species, &species->proj_init, species->f, t0);

  // we are pre-computing source for now as it is time-independent
  gk_neut_species_source_calc(app, species, &species->src, t0);
}

// Compute the RHS for species update, returning maximum stable
// time-step.
double
gk_neut_species_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  double omega_cfl = species->cfl_func(app, species, fin, rhs);
  return app->cfl/omega_cfl;
}

// Compute the RHS for species update, returning maximum stable
// time-step.
double
gk_neut_species_omega_cfl(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);
  
  gkyl_dg_updater_vlasov_advance(species->slvr, &species->local, 
    fin, species->cflrate, rhs);

  if (species->react_neut.num_react)
    gk_neut_species_react_rhs(app, species, &species->react_neut, fin, rhs);
  
  app->stat.nspecies_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);
  
  double omega_cfl_ho[1];
  if (app->use_gpu)
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    omega_cfl_ho[0] = species->omega_cfl[0];
  double omega_cfl = omega_cfl_ho[0];
  
  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return omega_cfl;
}

double
gk_neut_species_static_omega_cfl(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  return 1/DBL_MAX;
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for distribution function
void
gk_neut_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f)
{
  species->bc_func(app, species, f);
}

void
gk_neut_species_bc_func(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(species->comm, &species->local, &species->local_ext,
    num_periodic_dir, app->periodic_dirs, f); 
  
  for (int d=0; d<cdim; ++d) {
    if (species->bc_is_np[d]) {

      switch (species->lower_bc[d].type) {
        case GKYL_SPECIES_RECYCLE:
          ;
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer, f);
          break;
	case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer_lo_fixed, f);
          break;
        case GKYL_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        default:
          break;
      }

      switch (species->upper_bc[d].type) {
        case GKYL_SPECIES_RECYCLE:
          ;
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer_up_fixed, f);
          break;
        case GKYL_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        default:
          break;
      }      
    }
  }

  gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext, f);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_neut_species_static_bc_func(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f)
{
  // empty function
}

void
gk_neut_species_tm(gkyl_gyrokinetic_app *app)
{
  app->stat.species_rhs_tm = 0.0;
  for (int i=0; i<app->num_neut_species; ++i) {
    struct gkyl_dg_updater_vlasov_tm tm =
      gkyl_dg_updater_vlasov_get_tm(app->neut_species[i].slvr);
    app->stat.species_rhs_tm += tm.vlasov_tm;
  }
}

// release resources for species
void
gk_neut_species_release(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *s)
{
  // release various arrays
  gkyl_array_release(s->f);
  gk_neut_species_projection_release(app, &s->proj_init);
  gkyl_comm_release(s->comm);

  if (app->use_gpu)
    gkyl_array_release(s->f_host);

  gkyl_velocity_map_release(s->vel_map);

  // release moment data
  gk_neut_species_moment_release(app, &s->m0);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    gk_neut_species_moment_release(app, &s->moms[i]);
  gkyl_free(s->moms);
  gk_neut_species_moment_release(app, &s->integ_moms); 
  gkyl_dynvec_release(s->integ_diag);
  
  if (app->use_gpu) {
    gkyl_cu_free(s->omega_cfl);
    gkyl_cu_free(s->red_integ_diag);
    gkyl_cu_free(s->red_integ_diag_global);
  }
  else {
    gkyl_free(s->omega_cfl);
    gkyl_free(s->red_integ_diag);
    gkyl_free(s->red_integ_diag_global);
  }
  s->release_func(app, s);
}

// release resources for dynamic species
void
gk_neut_species_release_func(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *s)
{
  // release various arrays
  gkyl_array_release(s->bc_buffer);
  gkyl_array_release(s->bc_buffer_lo_fixed);
  gkyl_array_release(s->bc_buffer_up_fixed);

  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->cflrate);
  
  gkyl_array_release(s->alpha_surf);
  gkyl_array_release(s->sgn_alpha_surf);
  gkyl_array_release(s->const_sgn_alpha);
  gkyl_array_release(s->cot_vec);
  
  // release equation object and solver
  gkyl_dg_eqn_release(s->eqn_vlasov);
  gkyl_dg_updater_vlasov_release(s->slvr);

  gk_neut_species_source_release(app, &s->src);

  if (s->react_neut.num_react)
    gk_neut_species_react_release(app, &s->react_neut);

  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    if (s->lower_bc[d].type == GKYL_SPECIES_RECYCLE) 
      ;
    else 
      gkyl_bc_basic_release(s->bc_lo[d]);
    
    if (s->upper_bc[d].type == GKYL_SPECIES_RECYCLE) 
      ;
    else 
      gkyl_bc_basic_release(s->bc_up[d]);
  }
}

// release resources for species
void
gk_neut_species_static_release_func(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *s)
{
  // do nothing
}
