#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_vlasov_poisson_priv.h>
#include <gkyl_proj_on_basis.h>

#include <assert.h>
#include <time.h>

void
vp_species_init(struct gkyl_vp *vp, struct gkyl_vlasov_poisson_app *app, struct vp_species *vps)
{
  // Initialize species object.

  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = vp->cells[d];
    lower[d] = vp->lower[d];
    upper[d] = vp->upper[d];
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    // Full phase space grid.
    cells[cdim+d] = vps->info.cells[d];
    lower[cdim+d] = vps->info.lower[d];
    upper[cdim+d] = vps->info.upper[d];
    ghost[cdim+d] = 0; // No ghost-cells in velocity space.

    // Only velocity space.
    cells_vel[d] = vps->info.cells[d];
    lower_vel[d] = vps->info.lower[d];
    upper_vel[d] = vps->info.upper[d];
    ghost_vel[d] = 0; // No ghost-cells in velocity space.
  }
  // Full phase space grid.
  gkyl_rect_grid_init(&vps->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&vps->grid, ghost, &vps->global_ext, &vps->global);
  
  // Velocity space grid.
  gkyl_rect_grid_init(&vps->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&vps->grid_vel, ghost_vel, &vps->local_ext_vel, &vps->local_vel);

  // Phase-space communicator.
  vps->comm = gkyl_comm_extend_comm(app->comm, &vps->local_vel);

  // Create local and local_ext from app local range.
  struct gkyl_range local;
  // Local = conf-local X local_vel.
  gkyl_range_ten_prod(&local, &app->local, &vps->local_vel);
  gkyl_create_ranges(&local, ghost, &vps->local_ext, &vps->local);

  // Determine field-type.
  if (app->has_field)
    vps->field_id = app->field->info.field_id;
  else
    vps->field_id = GKYL_VP_FIELD_PHI;

  // Allocate distribution function arrays.
  vps->f = mkarr(app->use_gpu, app->basis.num_basis, vps->local_ext.volume);
  vps->f1 = mkarr(app->use_gpu, app->basis.num_basis, vps->local_ext.volume);
  vps->fnew = mkarr(app->use_gpu, app->basis.num_basis, vps->local_ext.volume);

  vps->f_host = vps->f;
  if (app->use_gpu)
    vps->f_host = mkarr(false, app->basis.num_basis, vps->local_ext.volume);

  // Allocate cflrate (scalar array).
  vps->cflrate = mkarr(app->use_gpu, 1, vps->local_ext.volume);

  if (app->use_gpu)
    vps->omega_cfl = gkyl_cu_malloc(sizeof(double));
  else 
    vps->omega_cfl = gkyl_malloc(sizeof(double));

  // Allocate arrays to store (q/m)*(phi,A_ext):
  vps->qbym = vps->info.charge/vps->info.mass;
  vps->qmem = mkarr(app->use_gpu, 4*app->confBasis.num_basis, app->local_ext.volume);

  // By default, we do not have zero-flux boundary conditions in any direction.
  bool is_zero_flux[2*GKYL_MAX_DIM] = {false};

  // Determine which directions are not periodic, if any directions are zero-flux,
  // need to set is_zero_flux.
  vps->num_periodic_dir = app->num_periodic_dir;
  for (int d=0; d<vps->num_periodic_dir; ++d)
    vps->periodic_dirs[d] = app->periodic_dirs[d];

  for (int d=0; d<app->cdim; ++d) vps->bc_is_np[d] = true;
  for (int d=0; d<vps->num_periodic_dir; ++d)
    vps->bc_is_np[vps->periodic_dirs[d]] = false;

  for (int dir=0; dir<app->cdim; ++dir) {
    vps->lower_bc[dir].type = vps->upper_bc[dir].type = GKYL_SPECIES_COPY;
    if (vps->bc_is_np[dir]) {
      const struct gkyl_vlasov_poisson_bcs *bc;
      if (dir == 0)
        bc = &vps->info.bcx;
      else if (dir == 1)
        bc = &vps->info.bcy;
      else
        bc = &vps->info.bcz;

      vps->lower_bc[dir] = bc->lower;
      vps->upper_bc[dir] = bc->upper;
      if (vps->lower_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        is_zero_flux[dir] = true;
      }
      if (vps->upper_bc[dir].type == GKYL_SPECIES_ZERO_FLUX) {
        is_zero_flux[dir+pdim] = true;
      }
    }
  }

  struct gkyl_dg_vlasov_poisson_auxfields aux_inp = {.field = vps->qmem, };
  // Create collisionless solver.
  vps->slvr = gkyl_dg_updater_vlasov_poisson_new(&vps->grid, &app->confBasis, &app->basis, 
    &app->local, &vps->local_vel, &vps->local, is_zero_flux, vps->model_id, vps->field_id, &aux_inp, app->use_gpu);

  // Acquire equation object.
  vps->eqn_vlasov = gkyl_dg_updater_vlasov_poisson_acquire_eqn(vps->slvr);

  // Allocate date for density.
  vp_species_moment_init(app, vps, &vps->m0, "M0");
  // Allocate data for integrated moments.
  vp_species_moment_init(app, vps, &vps->integ_moms, "Integrated");

  // Allocate data for diagnostic moments.
  int ndm = vps->info.num_diag_moments;
  vps->moms = gkyl_malloc(sizeof(struct vp_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    vp_species_moment_init(app, vps, &vps->moms[m], vps->info.diag_moments[m]);

  // Allocate arrays to reduce the integrated moments M0, M1i, M2.
  if (app->use_gpu) {
    vps->red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
    vps->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[vdim+2]));
  } else {
    vps->red_integ_diag = gkyl_malloc(sizeof(double[vdim+2]));
    vps->red_integ_diag_global = gkyl_malloc(sizeof(double[vdim+2]));
  }
  // Allocate dynamic-vector to store all-reduced integrated moments.
  vps->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  vps->is_first_integ_write_call = true;

  // Initialize projection routine for initial conditions.
  vp_species_projection_init(app, vps, vps->info.projection, &vps->proj_init);

  // Set species source id.
  vps->source_id = vps->info.source.source_id;
  
  // Determine collision type to use in Vlasov update.
  vps->collision_id = vps->info.collisions.collision_id;
  if (vps->collision_id == GKYL_LBO_COLLISIONS) {
    vp_species_lbo_init(app, vps, &vps->lbo);
  } 
  else if (vps->collision_id == GKYL_BGK_COLLISIONS) {
    vp_species_bgk_init(app, vps, &vps->bgk);
    // assert(false); // Not ready.
  }

  // Create ranges and allocate buffers for applying periodic and non-periodic BCs.
  long buff_sz = 0;
  // Compute buffer size needed.
  for (int dir=0; dir<cdim; ++dir) {
    // Create local lower skin and ghost ranges for distribution function.
    gkyl_skin_ghost_ranges(&vps->lower_skin[dir], &vps->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &vps->local_ext, ghost);
    // Create local upper skin and ghost ranges for distribution function.
    gkyl_skin_ghost_ranges(&vps->upper_skin[dir], &vps->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &vps->local_ext, ghost);

    long vol = GKYL_MAX2(vps->lower_skin[dir].volume, vps->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }

  // Initialize boundary fluxes for diagnostics and sources, for example.
  vp_species_bflux_init(app, vps, &vps->bflux);

  vps->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  // Buffer arrays for fixed function boundary conditions on distribution function.
  vps->bc_buffer_lo_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);
  vps->bc_buffer_up_fixed = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

  for (int d=0; d<cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (vps->lower_bc[d].type == GKYL_SPECIES_EMISSION) {
      vps->emit_lo = true;
      vp_species_emission_init(app, &vps->bc_emission_lo, d, GKYL_LOWER_EDGE, vps->lower_bc[d].aux_ctx,
        app->use_gpu);
    }
    else {  
      if (vps->lower_bc[d].type == GKYL_SPECIES_COPY)
        bctype = GKYL_BC_COPY;
      else if (vps->lower_bc[d].type == GKYL_SPECIES_ABSORB)
        bctype = GKYL_BC_ABSORB;
      else if (vps->lower_bc[d].type == GKYL_SPECIES_REFLECT)
        bctype = GKYL_BC_REFLECT;
      else if (vps->lower_bc[d].type == GKYL_SPECIES_FIXED_FUNC)
        bctype = GKYL_BC_FIXED_FUNC;

      vps->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev.basis,
        &vps->lower_skin[d], &vps->lower_ghost[d], vps->f->ncomp, app->cdim, app->use_gpu);
    }

    // Upper BC updater. Copy BCs by default.
    if (vps->upper_bc[d].type == GKYL_SPECIES_EMISSION) {
      vps->emit_up = true;
      vp_species_emission_init(app, &vps->bc_emission_up, d, GKYL_UPPER_EDGE, vps->upper_bc[d].aux_ctx,
        app->use_gpu);
    }
    else {
      if (vps->upper_bc[d].type == GKYL_SPECIES_COPY)
        bctype = GKYL_BC_COPY;
      else if (vps->upper_bc[d].type == GKYL_SPECIES_ABSORB)
        bctype = GKYL_BC_ABSORB;
      else if (vps->upper_bc[d].type == GKYL_SPECIES_REFLECT)
        bctype = GKYL_BC_REFLECT;
      else if (vps->upper_bc[d].type == GKYL_SPECIES_FIXED_FUNC)
        bctype = GKYL_BC_FIXED_FUNC;

      vps->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.basis,
        &vps->upper_skin[d], &vps->upper_ghost[d], vps->f->ncomp, app->cdim, app->use_gpu);
    }
  }
}

void
vp_species_apply_ic(gkyl_vlasov_poisson_app *app, struct vp_species *species, double t0)
{
  vp_species_projection_calc(app, species, &species->proj_init, species->f, t0);

  // we are pre-computing source for now as it is time-independent
  vp_species_source_calc(app, species, &species->src, t0);

  vp_species_bflux_rhs(app, species, &species->bflux, species->f, species->f);

    // Optional runtime configuration to use BGK collisions but with fixed input 
  // temperature relaxation based on the initial temperature value. 
  if (species->bgk.fixed_temp_relax) {
    vp_species_bgk_moms_fixed_temp(app, species, &species->bgk, species->f);
  }

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
vp_species_rhs(gkyl_vlasov_poisson_app *app, struct vp_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  if (app->has_field)
    gkyl_array_set_offset(species->qmem, species->qbym, app->field->phi, 0);

  gkyl_dg_updater_vlasov_poisson_advance(species->slvr, &species->local, 
    fin, species->cflrate, rhs);

  if (species->collision_id == GKYL_LBO_COLLISIONS)
    vp_species_lbo_rhs(app, species, &species->lbo, fin, rhs);
   else if (species->collision_id == GKYL_BGK_COLLISIONS)
     vp_species_bgk_rhs(app, species, &species->bgk, fin, rhs);
  
  vp_species_bflux_rhs(app, species, &species->bflux, fin, rhs);

  app->stat.nspecies_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);

  double omegaCfl_ho[1];
  if (app->use_gpu)
    gkyl_cu_memcpy(omegaCfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    omegaCfl_ho[0] = species->omega_cfl[0];
  double omegaCfl = omegaCfl_ho[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  
  return app->cfl/omegaCfl;
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for distribution function
void
vp_species_apply_bc(gkyl_vlasov_poisson_app *app, const struct vp_species *species, struct gkyl_array *f, double tcurr)
{
  struct timespec wst = gkyl_wall_clock();
  
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(species->comm, &species->local, &species->local_ext,
    num_periodic_dir, app->periodic_dirs, f); 
  
  for (int d=0; d<cdim; ++d) {
    if (species->bc_is_np[d]) {

      switch (species->lower_bc[d].type) {
        case GKYL_SPECIES_EMISSION:
          vp_species_emission_apply_bc(app, &species->bc_emission_lo, f, tcurr);
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_lo[d], species->bc_buffer_lo_fixed, f);
          break;
        default:
          break;
      }

      switch (species->upper_bc[d].type) {
        case GKYL_SPECIES_EMISSION:
          vp_species_emission_apply_bc(app, &species->bc_emission_up, f, tcurr);
          break;
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_REFLECT:
        case GKYL_SPECIES_ABSORB:
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer, f);
          break;
        case GKYL_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(species->bc_up[d], species->bc_buffer_up_fixed, f);
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
vp_species_coll_tm(gkyl_vlasov_poisson_app *app)
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
vp_species_tm(gkyl_vlasov_poisson_app *app)
{
  app->stat.species_rhs_tm = 0.0;
  for (int i=0; i<app->num_species; ++i) {
    struct gkyl_dg_updater_vlasov_poisson_tm tm =
      gkyl_dg_updater_vlasov_poisson_get_tm(app->species[i].slvr);
    app->stat.species_rhs_tm += tm.vlasov_poisson_tm;
  }
}

// release resources for species
void
vp_species_release(const gkyl_vlasov_poisson_app* app, const struct vp_species *s)
{
  // release various arrays
  gkyl_array_release(s->f);
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->cflrate);
  gkyl_array_release(s->bc_buffer);
  gkyl_array_release(s->bc_buffer_lo_fixed);
  gkyl_array_release(s->bc_buffer_up_fixed);

  vp_species_projection_release(app, &s->proj_init);

  vp_species_bflux_release(app, &s->bflux);

  gkyl_comm_release(s->comm);

  if (app->use_gpu)
    gkyl_array_release(s->f_host);

  gkyl_array_release(s->qmem);

  // release equation object and solver
  gkyl_dg_eqn_release(s->eqn_vlasov);
  gkyl_dg_updater_vlasov_poisson_release(s->slvr);

  // release moment data
  vp_species_moment_release(app, &s->m0);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    vp_species_moment_release(app, &s->moms[i]);
  gkyl_free(s->moms);
  vp_species_moment_release(app, &s->integ_moms); 
  gkyl_dynvec_release(s->integ_diag);

  if (s->source_id) {
    vp_species_source_release(app, &s->src);
  }

  if (s->collision_id == GKYL_LBO_COLLISIONS)
    vp_species_lbo_release(app, &s->lbo);
   else if (s->collision_id == GKYL_BGK_COLLISIONS)
     vp_species_bgk_release(app, &s->bgk);

  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    if (s->lower_bc[d].type == GKYL_SPECIES_EMISSION)
      vp_species_emission_release(&s->bc_emission_lo);
    else 
      gkyl_bc_basic_release(s->bc_lo[d]);
    
    if (s->upper_bc[d].type == GKYL_SPECIES_EMISSION)
      vp_species_emission_release(&s->bc_emission_up);
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
}
