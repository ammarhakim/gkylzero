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

  // allocate distribution function arrays
  s->f = mkarr(app->use_gpu, app->neut_basis.num_basis, s->local_ext.volume);
  s->f1 = mkarr(app->use_gpu, app->neut_basis.num_basis, s->local_ext.volume);
  s->fnew = mkarr(app->use_gpu, app->neut_basis.num_basis, s->local_ext.volume);

  s->f_host = s->f;
  if (app->use_gpu)
    s->f_host = mkarr(false, app->neut_basis.num_basis, s->local_ext.volume);

  // allocate cflrate (scalar array)
  s->cflrate = mkarr(app->use_gpu, 1, s->local_ext.volume);

  if (app->use_gpu)
    s->omegaCfl_ptr = gkyl_cu_malloc(sizeof(double));
  else 
    s->omegaCfl_ptr = gkyl_malloc(sizeof(double));

  // by default, we do not have zero-flux boundary conditions in any direction
  bool is_zero_flux[GKYL_MAX_DIM] = {false};

  // determine which directions are not periodic, if any directions are zero-flux, need to set is_zero_flux
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
      // Zero flux BCs can only be applied jointly on lower and upper boundary
      if (s->lower_bc[dir] == GKYL_SPECIES_ZERO_FLUX || s->upper_bc[dir] == GKYL_SPECIES_ZERO_FLUX)
        is_zero_flux[dir] = true;
    }
  }

  struct gkyl_dg_vlasov_auxfields aux_inp = {.field = 0, .cot_vec = app->gk_geom->dxdz, .alpha_geo = 0};
  // Set field type and model id for neutral species in GK system and create solver
  s->field_id = GKYL_FIELD_NULL;
  s->model_id = GKYL_MODEL_GEN_GEO;
  s->slvr = gkyl_dg_updater_vlasov_new(&s->grid, &app->confBasis, &app->neut_basis, 
    &app->local, &s->local_vel, &s->local, is_zero_flux, s->model_id, s->field_id, &aux_inp, app->use_gpu);

  // acquire equation object
  s->eqn_vlasov = gkyl_dg_updater_vlasov_acquire_eqn(s->slvr);

  // allocate date for density 
  gk_neut_species_moment_init(app, s, &s->m0, "M0");
  // allocate data for integrated moments
  gk_neut_species_moment_init(app, s, &s->integ_moms, "Integrated");

  // allocate data for diagnostic moments
  int ndm = s->info.num_diag_moments;
  s->moms = gkyl_malloc(sizeof(struct gk_neut_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    gk_neut_species_moment_init(app, s, &s->moms[m], s->info.diag_moments[m]);

  if (app->use_gpu) 
    s->red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
  // allocate dynamic-vector to store all-reduced integrated moments 
  s->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  s->is_first_integ_write_call = true;

  // set species source id
  s->source_id = s->info.source.source_id;

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
  s->bc_buffer = mkarr(app->use_gpu, app->neut_basis.num_basis, buff_sz);
  // buffer arrays for fixed function boundary conditions on distribution function
  s->bc_buffer_lo_fixed = mkarr(app->use_gpu, app->neut_basis.num_basis, buff_sz);
  s->bc_buffer_up_fixed = mkarr(app->use_gpu, app->neut_basis.num_basis, buff_sz);

  for (int d=0; d<cdim; ++d) {
    // Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (s->lower_bc[d] == GKYL_SPECIES_RECYCLE) {
      ;
    }
    else { 
      if (s->lower_bc[d] == GKYL_SPECIES_COPY) 
        bctype = GKYL_BC_COPY;
      else if (s->lower_bc[d] == GKYL_SPECIES_ABSORB) 
        bctype = GKYL_BC_ABSORB;
      else if (s->lower_bc[d] == GKYL_SPECIES_REFLECT) 
        bctype = GKYL_BC_REFLECT;
      else if (s->lower_bc[d] == GKYL_SPECIES_FIXED_FUNC) 
        bctype = GKYL_BC_FIXED_FUNC;

      s->bc_lo[d] = gkyl_bc_basic_new(d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev.neut_basis,
        &s->lower_skin[d], &s->lower_ghost[d], s->f->ncomp, app->cdim, app->use_gpu);
    }

    if (s->upper_bc[d] == GKYL_SPECIES_RECYCLE) {
      ;
    }
    else {
      // Upper BC updater. Copy BCs by default.
      if (s->upper_bc[d] == GKYL_SPECIES_COPY) 
        bctype = GKYL_BC_COPY;
      else if (s->upper_bc[d] == GKYL_SPECIES_ABSORB) 
        bctype = GKYL_BC_ABSORB;
      else if (s->upper_bc[d] == GKYL_SPECIES_REFLECT) 
        bctype = GKYL_BC_REFLECT;
      else if (s->upper_bc[d] == GKYL_SPECIES_FIXED_FUNC) 
        bctype = GKYL_BC_FIXED_FUNC;

      s->bc_up[d] = gkyl_bc_basic_new(d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev.neut_basis,
        &s->upper_skin[d], &s->upper_ghost[d], s->f->ncomp, app->cdim, app->use_gpu);
    }
  }
}

void
gk_neut_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_neut_species *species, double t0)
{
  int poly_order = app->poly_order;
  if (species->info.is_maxwellian){
    // Project n, udrift, and vt^2 based on input functions
    struct gkyl_array *m0 = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    struct gkyl_array *udrift = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    struct gkyl_array *vtsq = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      poly_order+1, 1, species->info.init_density, species->info.ctx_density);
    gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      poly_order+1, 1, species->info.init_upar, species->info.ctx_upar);
    gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      poly_order+1, 1, species->info.init_temp, species->info.ctx_temp);

    gkyl_proj_on_basis_advance(proj_m0, 0.0, &app->local_ext, m0); 
    gkyl_proj_on_basis_advance(proj_udrift, 0.0, &app->local_ext, udrift);
    gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &app->local_ext, vtsq);
    gkyl_array_scale(vtsq, 1/species->info.mass);

    // proj_maxwellian expects the primitive moments as a single array.
    struct gkyl_array *prim_moms = mkarr(false, 2*app->confBasis.num_basis, app->local_ext.volume);
    gkyl_array_set_offset(prim_moms, 1.0, udrift, 0*app->confBasis.num_basis);
    gkyl_array_set_offset(prim_moms, 1.0, vtsq  , 1*app->confBasis.num_basis);

    // Initialize Maxwellian projection object
    gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_new(&species->grid,
        &app->confBasis, &app->neut_basis, poly_order+1, app->use_gpu);

    // If on GPUs, need to copy n, udrift, and vt^2 onto device
    struct gkyl_array *prim_moms_dev, *m0_dev;
    if (app->use_gpu) {
      prim_moms_dev = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
      m0_dev = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

      gkyl_array_copy(prim_moms_dev, prim_moms);
      gkyl_array_copy(m0_dev, m0);
      gkyl_proj_maxwellian_on_basis_prim_mom(proj_max, &species->local_ext, &app->local_ext, 
        m0_dev, prim_moms_dev, species->f);
    }
    else {
      gkyl_proj_maxwellian_on_basis_prim_mom(proj_max, &species->local_ext, &app->local_ext, 
        m0_dev, prim_moms_dev, species->f);
    }
    // Now compute and scale the density to the desired density function based on input density from Maxwellian projection
    gk_neut_species_moment_calc(&species->m0, species->local_ext, app->local_ext, species->f); 

    // Rescale projected density to desired input density function
    struct gkyl_array *m0mod = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    struct gkyl_dg_bin_op_mem *mem;
    if (app->use_gpu) {
      mem = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->confBasis.num_basis);
      gkyl_dg_div_op_range(mem, app->confBasis, 0, m0mod, 0, m0_dev, 0, species->m0.marr, &app->local);
    }
    else {
      mem = gkyl_dg_bin_op_mem_new(app->local.volume, app->confBasis.num_basis);
      gkyl_dg_div_op_range(mem, app->confBasis, 0, m0mod, 0, m0, 0, species->m0.marr, &app->local);
    }
    gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->neut_basis, species->f, 
        m0mod, species->f, &app->local_ext, &species->local_ext);

    // multiply final distribution function by Jacobian
    gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->neut_basis, species->f, 
        app->gk_geom->jacobgeo, species->f, &app->local_ext, &species->local_ext);

    // Free temporary variables and projection objects
    gkyl_array_release(m0);
    gkyl_array_release(udrift); 
    gkyl_array_release(vtsq);
    gkyl_array_release(prim_moms);
    if (app->use_gpu) {
      gkyl_array_release(m0_dev);
      gkyl_array_release(prim_moms_dev);      
    }
    gkyl_proj_on_basis_release(proj_m0);
    gkyl_proj_on_basis_release(proj_udrift);
    gkyl_proj_on_basis_release(proj_vtsq);
    gkyl_array_release(m0mod); 
    gkyl_dg_bin_op_mem_release(mem);
    gkyl_proj_maxwellian_on_basis_release(proj_max);
  }
  else {
    gkyl_proj_on_basis *proj;
    proj = gkyl_proj_on_basis_new(&species->grid, &app->neut_basis,
      poly_order+1, 1, species->info.init_dist, species->info.ctx_dist);

    gkyl_proj_on_basis_advance(proj, t0, &species->local, species->f_host);
    gkyl_proj_on_basis_release(proj);    

    if (app->use_gpu) // note: f_host is same as f when not on GPUs
      gkyl_array_copy(species->f, species->f_host);

  }

  // we are pre-computing source for now as it is time-independent
  if (species->source_id)
    gk_neut_species_source_calc(app, species, t0);

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
gk_neut_species_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  gkyl_dg_updater_vlasov_advance(species->slvr, &species->local, 
    fin, species->cflrate, rhs);

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
gk_neut_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f)
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

      switch (species->lower_bc[d]) {
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
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
          break;
        case GKYL_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        default:
          break;
      }

      switch (species->upper_bc[d]) {
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
        case GKYL_SPECIES_NO_SLIP:
        case GKYL_SPECIES_WEDGE:
          assert(false);
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
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->cflrate);
  gkyl_array_release(s->bc_buffer);
  if (app->cdim > 1) {
    gkyl_array_release(s->bc_buffer_lo_fixed);
    gkyl_array_release(s->bc_buffer_up_fixed);
  }

  gkyl_comm_release(s->comm);

  if (app->use_gpu)
    gkyl_array_release(s->f_host);

  // release equation object and solver
  gkyl_dg_eqn_release(s->eqn_vlasov);
  gkyl_dg_updater_vlasov_release(s->slvr);

  // release moment data
  gk_neut_species_moment_release(app, &s->m0);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    gk_neut_species_moment_release(app, &s->moms[i]);
  gkyl_free(s->moms);
  gk_neut_species_moment_release(app, &s->integ_moms); 
  gkyl_dynvec_release(s->integ_diag);

  if (s->source_id) {
    gk_neut_species_source_release(app, &s->src);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d=0; d<app->cdim; ++d) {
    if (s->lower_bc[d] == GKYL_SPECIES_RECYCLE) 
      ;
    else 
      gkyl_bc_basic_release(s->bc_lo[d]);
    
    if (s->upper_bc[d] == GKYL_SPECIES_RECYCLE) 
      ;
    else 
      gkyl_bc_basic_release(s->bc_up[d]);
  }
  
  if (app->use_gpu) {
    gkyl_cu_free(s->omegaCfl_ptr);
    gkyl_cu_free(s->red_integ_diag);
  }
  else {
    gkyl_free(s->omegaCfl_ptr);
  }
}