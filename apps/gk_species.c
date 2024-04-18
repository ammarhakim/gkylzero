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

// initialize species object
void
gk_species_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_species *s)
{
  int cdim = app->cdim, vdim = app->vdim;
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

  // determine field-type 
  s->gkfield_id = app->field->gkfield_id;
  if (s->info.no_by) {
    s->gkmodel_id = GKYL_GK_MODEL_NO_BY;
  }
  else {
    s->gkmodel_id = GKYL_GK_MODEL_GEN_GEO;
  }

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
    s->omega_cfl = gkyl_cu_malloc(sizeof(double));
  else 
    s->omega_cfl = gkyl_malloc(sizeof(double));

  // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
  struct gkyl_basis surf_basis, surf_quad_basis;
  if (app->poly_order > 1) {
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, app->poly_order);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, app->poly_order);
  }
  else {
    gkyl_cart_modal_gkhybrid(&surf_basis, cdim-1, vdim); // p=2 in vparallel
    gkyl_cart_modal_gkhybrid(&surf_quad_basis, cdim-1, vdim); 
  }
  int alpha_surf_sz = (cdim+1)*surf_basis.num_basis;
  int sgn_alpha_surf_sz = (cdim+1)*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points

  // allocate arrays to store fields: 
  // 1. alpha_surf (surface phase space flux)
  // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
  // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
  s->alpha_surf = mkarr(app->use_gpu, alpha_surf_sz, s->local_ext.volume);
  s->sgn_alpha_surf = mkarr(app->use_gpu, sgn_alpha_surf_sz, s->local_ext.volume);
  s->const_sgn_alpha = mk_int_arr(app->use_gpu, (cdim+1), s->local_ext.volume);
  // 4. EM fields: phi and (if EM GK) Apar and d/dt Apar  
  s->phi = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  s->apar = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  s->apardot = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);    

  s->calc_gk_vars = gkyl_dg_calc_gyrokinetic_vars_new(&s->grid, &app->confBasis, &app->basis, 
    s->info.charge, s->info.mass, s->gkmodel_id, app->gk_geom, app->use_gpu);

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
      if (s->lower_bc[dir].type == GKYL_SPECIES_GK_IWL || s->upper_bc[dir].type == GKYL_SPECIES_GK_IWL) {
        // Make the parallel direction periodic so that we sync the core before
        // applying sheath BCs in the SOL.
        s->periodic_dirs[s->num_periodic_dir] = app->cdim-1; // The last direction is the parallel one.
        s->num_periodic_dir += 1;
        // Check that the LCFS is the same on both BCs and that it's on a cell boundary within our grid.
        double xLCFS = s->lower_bc[dir].aux_parameter;
        assert(fabs(xLCFS-s->upper_bc[dir].aux_parameter) < 1e-14);
        // Assume the split happens at a cell boundary and within the domain.
        assert((app->grid.lower[0]<xLCFS) && (xLCFS<app->grid.upper[0]));
        double needint = (xLCFS-app->grid.lower[0])/app->grid.dx[0];
        assert(floor(fabs(needint-floor(needint))) < 1.);
      }
    }
  }

  struct gkyl_dg_gyrokinetic_auxfields aux_inp = { .alpha_surf = s->alpha_surf, 
    .sgn_alpha_surf = s->sgn_alpha_surf, .const_sgn_alpha = s->const_sgn_alpha, 
    .phi = s->phi, .apar = s->apar, .apardot = s->apardot };
  // create solver
  s->slvr = gkyl_dg_updater_gyrokinetic_new(&s->grid, &app->confBasis, &app->basis, 
    &app->local, &s->local, is_zero_flux, 
    s->info.charge, s->info.mass, s->gkmodel_id, 
    app->gk_geom, &aux_inp, app->use_gpu);

  // acquire equation object
  s->eqn_gyrokinetic = gkyl_dg_updater_gyrokinetic_acquire_eqn(s->slvr);

  // allocate data for density (for use in charge density accumulation and weak division for upar)
  gk_species_moment_init(app, s, &s->m0, "M0");
  // allocate data for integrated moments
  gk_species_moment_init(app, s, &s->integ_moms, "Integrated");

  // allocate data for diagnostic moments
  int ndm = s->info.num_diag_moments;
  s->moms = gkyl_malloc(sizeof(struct gk_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    gk_species_moment_init(app, s, &s->moms[m], s->info.diag_moments[m]);

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
  gk_species_projection_init(app, s, s->info.projection, &s->proj_init);

  // set species source id
  s->source_id = s->info.source.source_id;
  
  // determine collision type to use in gyrokinetic update
  s->collision_id = s->info.collisions.collision_id;
  if (s->collision_id == GKYL_LBO_COLLISIONS) 
    gk_species_lbo_init(app, s, &s->lbo);
  else if (s->collision_id == GKYL_BGK_COLLISIONS)
    gk_species_bgk_init(app, s, &s->bgk);

  s->has_reactions = false;
  s->has_neutral_reactions = false;
  if (s->info.react.num_react) {
    s->has_reactions = true;
    gk_species_react_init(app, s, s->info.react, &s->react, true);
  }
  if (s->info.react_neut.num_react) {
    s->has_neutral_reactions = true;
    gk_species_react_init(app, s, s->info.react_neut, &s->react_neut, false);
  }

  // determine radiation type to use in gyrokinetic update
  s->radiation_id = s->info.radiation.radiation_id;

  // initialize boundary fluxes for diagnostics and, if present,
  // ambipolar potential solve
  gk_species_bflux_init(app, s, &s->bflux); 

  // initialize diffusion if present
  s->has_diffusion = false;  
  if (s->info.diffusion.num_diff_dir) {
    s->has_diffusion = true;
    int diffusion_order = s->info.diffusion.order ? s->info.diffusion.order : 2;

    int szD = cdim*app->confBasis.num_basis;
    s->diffD = mkarr(app->use_gpu, szD, app->local_ext.volume);
    bool diff_dir[GKYL_MAX_CDIM] = {false};

    int num_diff_dir = s->info.diffusion.num_diff_dir ? s->info.diffusion.num_diff_dir : app->cdim;
    // Assuming diffusion along x only for now.
    assert(num_diff_dir == 1);
    assert(s->info.diffusion.diff_dirs[0] == 0);
    for (int d=0; d<num_diff_dir; ++d) {
      int dir = s->info.diffusion.diff_dirs[d]; 
      diff_dir[dir] = 1; 
      gkyl_array_shiftc(s->diffD, s->info.diffusion.D[d]*pow(sqrt(2),app->cdim), dir);
    }
    // Multiply diffD by g^xx*jacobgeo.
    gkyl_dg_mul_op(app->confBasis, 0, s->diffD, 0, app->gk_geom->gxxj, 0, s->diffD);

    s->diff_slvr = gkyl_dg_updater_diffusion_gyrokinetic_new(&s->grid, &app->basis, &app->confBasis, 
      false, diff_dir, diffusion_order, &app->local, is_zero_flux, app->use_gpu);
  }

  // create ranges and allocate buffers for applying periodic and non-periodic BCs
  long buff_sz = 0;
  // compute buffer size needed
  for (int dir=0; dir<cdim; ++dir) {
    // Create local lower skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->lower_skin[dir], &s->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &s->local_ext, ghost);
    // Create local upper skin and ghost ranges for distribution function
    gkyl_skin_ghost_ranges(&s->upper_skin[dir], &s->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &s->local_ext, ghost);

    long vol = GKYL_MAX2(s->lower_skin[dir].volume, s->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  s->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  // buffer arrays for fixed function boundary conditions on distribution function
  s->bc_buffer_lo_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  s->bc_buffer_up_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

  for (int d=0; d<cdim; ++d) {
    // Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (s->lower_bc[d].type == GKYL_SPECIES_GK_SHEATH) {
      s->bc_sheath_lo = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_LOWER_EDGE, app->basis_on_dev.basis, 
        &s->lower_skin[d], &s->lower_ghost[d], &s->grid, cdim, 2.0*(s->info.charge/s->info.mass), app->use_gpu);
    }
    else if (s->lower_bc[d].type == GKYL_SPECIES_GK_IWL) {
      double xLCFS = s->lower_bc[d].aux_parameter;
      // Index of the cell that abuts the xLCFS from below.
      int idxLCFS_m = (xLCFS-1e-8 - app->grid.lower[0])/app->grid.dx[0]+1;
      gkyl_range_shorten_from_below(&s->lower_skin_par_sol, &s->lower_skin[d], 0, app->grid.cells[0]-idxLCFS_m+1);
      gkyl_range_shorten_from_below(&s->lower_ghost_par_sol, &s->lower_ghost[d], 0, app->grid.cells[0]-idxLCFS_m+1);

      s->bc_sheath_lo = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_LOWER_EDGE, app->basis_on_dev.basis, 
        &s->lower_skin_par_sol, &s->lower_ghost_par_sol, &s->grid, cdim, 2.0*(s->info.charge/s->info.mass), app->use_gpu);
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

    if (s->upper_bc[d].type == GKYL_SPECIES_GK_SHEATH) {
      s->bc_sheath_up = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_UPPER_EDGE, app->basis_on_dev.basis, 
        &s->upper_skin[d], &s->upper_ghost[d], &s->grid, cdim, 2.0*(s->info.charge/s->info.mass), app->use_gpu);
    }
    else if (s->lower_bc[d].type == GKYL_SPECIES_GK_IWL) {
      double xLCFS = s->upper_bc[d].aux_parameter;
      // Index of the cell that abuts the xLCFS from below.
      int idxLCFS_m = (xLCFS-1e-8 - app->grid.lower[0])/app->grid.dx[0]+1;
      gkyl_range_shorten_from_below(&s->upper_skin_par_sol, &s->upper_skin[d], 0, app->grid.cells[0]-idxLCFS_m+1);
      gkyl_range_shorten_from_below(&s->upper_ghost_par_sol, &s->upper_ghost[d], 0, app->grid.cells[0]-idxLCFS_m+1);

      s->bc_sheath_up = gkyl_bc_sheath_gyrokinetic_new(d, GKYL_UPPER_EDGE, app->basis_on_dev.basis, 
        &s->upper_skin_par_sol, &s->upper_ghost_par_sol, &s->grid, cdim, 2.0*(s->info.charge/s->info.mass), app->use_gpu);
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

      s->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.basis,
        &s->upper_skin[d], &s->upper_ghost[d], s->f->ncomp, app->cdim, app->use_gpu);
    }
  }

  // Positivity enforcing by shifting f.
  s->enforce_positivity = false;
  if (s->info.enforce_positivity) {
    s->enforce_positivity = true;
    s->pos_shift_op = gkyl_positivity_shift_gyrokinetic_new(cdim, app->basis, s->grid, s->info.mass, app->gk_geom, app->use_gpu);
    s->ps_intmom_grid = mkarr(app->use_gpu, vdim+2, app->local_ext.volume);
    s->ps_integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
    s->is_first_ps_integ_write_call = true;
  }
}

void
gk_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_species *species, double t0)
{
  gk_species_projection_calc(app, species, &species->proj_init, species->f, t0);

  // we are pre-computing source for now as it is time-independent
  if (species->source_id)
    gk_species_source_calc(app, species, &species->src, t0);

  // copy contents of initial conditions into buffer if specific BCs require them
  // *only works in x dimension for now for cdim > 1*
  if (app->cdim>1) {
    gkyl_bc_basic_buffer_fixed_func(species->bc_lo[0], species->bc_buffer_lo_fixed, species->f);
    gkyl_bc_basic_buffer_fixed_func(species->bc_up[0], species->bc_buffer_up_fixed, species->f);
  }
}

// Compute the RHS for species update, returning maximum stable
// time-step.
double
gk_species_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkyl_array_set(species->phi, 1.0, app->field->phi_smooth);

  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Compute the surface expansion of the phase space flux
  // Note: Each cell stores the *lower* surface expansions of the 
  // phase space flux, so local_ext range needed to index the output
  // values of alpha_surf even though we only loop over local ranges
  // to avoid evaluating quantities such as geometry in ghost cells
  // where they are not defined.
  gkyl_dg_calc_gyrokinetic_vars_alpha_surf(species->calc_gk_vars, 
    &app->local, &species->local, &species->local_ext, 
    species->phi, species->alpha_surf, species->sgn_alpha_surf, species->const_sgn_alpha);

  gkyl_dg_updater_gyrokinetic_advance(species->slvr, &species->local, 
    fin, species->cflrate, rhs);

  if (species->collision_id == GKYL_LBO_COLLISIONS)
    gk_species_lbo_rhs(app, species, &species->lbo, fin, rhs);
  else if (species->collision_id == GKYL_BGK_COLLISIONS)
    gk_species_bgk_rhs(app, species, &species->bgk, fin, rhs);
  
  if (species->has_diffusion)
    gkyl_dg_updater_diffusion_gyrokinetic_advance(species->diff_slvr, &species->local, 
      species->diffD, app->gk_geom->jacobgeo_inv, fin, species->cflrate, rhs);

  if (species->radiation_id == GKYL_GK_RADIATION)
    gk_species_radiation_rhs(app, species, &species->rad, fin, rhs);

  if (species->has_reactions)
    gk_species_react_rhs(app, species, &species->react, fin, rhs);
  if (species->has_neutral_reactions)
    gk_species_react_rhs(app, species, &species->react_neut, fin, rhs);
  
  app->stat.nspecies_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);

  double omega_cfl_ho[1];
  if (app->use_gpu)
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    omega_cfl_ho[0] = species->omega_cfl[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  
  return app->cfl/omega_cfl_ho[0];
}

// Apply boundary conditions to the distribution function.
void
gk_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_species *species, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = species->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(species->comm, &species->local, &species->local_ext,
    num_periodic_dir, species->periodic_dirs, f); 
  
  for (int d=0; d<cdim; ++d) {
    if (species->bc_is_np[d]) {

      switch (species->lower_bc[d].type) {
        case GKYL_SPECIES_GK_SHEATH:
        case GKYL_SPECIES_GK_IWL:
          gkyl_bc_sheath_gyrokinetic_advance(species->bc_sheath_lo, app->field->phi_smooth, 
            app->field->phi_wall_lo, f, &app->local);
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
          break;
        default:
          break;
      }

      switch (species->upper_bc[d].type) {
        case GKYL_SPECIES_GK_SHEATH:
        case GKYL_SPECIES_GK_IWL:
          gkyl_bc_sheath_gyrokinetic_advance(species->bc_sheath_up, app->field->phi_smooth, 
            app->field->phi_wall_up, f, &app->local);
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
gk_species_coll_tm(gkyl_gyrokinetic_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      struct gkyl_dg_updater_lbo_gyrokinetic_tm tm =
        gkyl_dg_updater_lbo_gyrokinetic_get_tm(app->species[i].lbo.coll_slvr);
      app->stat.species_lbo_coll_diff_tm[i] = tm.diff_tm;
      app->stat.species_lbo_coll_drag_tm[i] = tm.drag_tm;
    }
  }
}

void
gk_species_tm(gkyl_gyrokinetic_app *app)
{
  app->stat.species_rhs_tm = 0.0;
  for (int i=0; i<app->num_species; ++i) {
    struct gkyl_dg_updater_gyrokinetic_tm tm =
      gkyl_dg_updater_gyrokinetic_get_tm(app->species[i].slvr);
    app->stat.species_rhs_tm += tm.gyrokinetic_tm;
  }
}

// release resources for species
void
gk_species_release(const gkyl_gyrokinetic_app* app, const struct gk_species *s)
{
  // release various arrays and species objects
  gkyl_array_release(s->f);
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->cflrate);
  gkyl_array_release(s->bc_buffer);
  gkyl_array_release(s->bc_buffer_lo_fixed);
  gkyl_array_release(s->bc_buffer_up_fixed);

  gk_species_projection_release(app, &s->proj_init);

  gkyl_comm_release(s->comm);

  if (app->use_gpu)
    gkyl_array_release(s->f_host);

  gkyl_array_release(s->alpha_surf);
  gkyl_array_release(s->sgn_alpha_surf);
  gkyl_array_release(s->const_sgn_alpha);
  gkyl_array_release(s->phi);
  gkyl_array_release(s->apar);
  gkyl_array_release(s->apardot);

  gkyl_dg_calc_gyrokinetic_vars_release(s->calc_gk_vars);

  // release equation object and solver
  gkyl_dg_eqn_release(s->eqn_gyrokinetic);
  gkyl_dg_updater_gyrokinetic_release(s->slvr);

  // release moment data
  gk_species_moment_release(app, &s->m0);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    gk_species_moment_release(app, &s->moms[i]);
  gkyl_free(s->moms);
  gk_species_moment_release(app, &s->integ_moms); 
  gkyl_dynvec_release(s->integ_diag);

  if (s->source_id) {
    gk_species_source_release(app, &s->src);
  }

  if (s->collision_id == GKYL_LBO_COLLISIONS)
    gk_species_lbo_release(app, &s->lbo);
  else if (s->collision_id == GKYL_BGK_COLLISIONS)
    gk_species_bgk_release(app, &s->bgk);

  if (s->has_reactions)
    gk_species_react_release(app, &s->react);
  if (s->has_neutral_reactions)
    gk_species_react_release(app, &s->react_neut);  

  if (s->radiation_id == GKYL_GK_RADIATION) 
    gk_species_radiation_release(app, &s->rad);

  gk_species_bflux_release(app, &s->bflux);

  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    if ((s->lower_bc[d].type == GKYL_SPECIES_GK_SHEATH) ||
        (s->lower_bc[d].type == GKYL_SPECIES_GK_IWL))
      gkyl_bc_sheath_gyrokinetic_release(s->bc_sheath_lo);
    else 
      gkyl_bc_basic_release(s->bc_lo[d]);
    
    if ((s->upper_bc[d].type == GKYL_SPECIES_GK_SHEATH) ||
        (s->upper_bc[d].type == GKYL_SPECIES_GK_IWL))
      gkyl_bc_sheath_gyrokinetic_release(s->bc_sheath_up);
    else 
      gkyl_bc_basic_release(s->bc_up[d]);
  }
  
  if (app->use_gpu) {
    gkyl_cu_free(s->omega_cfl);
    gkyl_cu_free(s->red_integ_diag);
    gkyl_cu_free(s->red_integ_diag_global);
  }
  else {
    gkyl_free(s->red_integ_diag);
    gkyl_free(s->red_integ_diag_global);
    gkyl_free(s->omega_cfl);
  }

  if (s->enforce_positivity) {
    gkyl_positivity_shift_gyrokinetic_release(s->pos_shift_op);
    gkyl_array_release(s->ps_intmom_grid);
    gkyl_dynvec_release(s->ps_integ_diag);
  }
}
