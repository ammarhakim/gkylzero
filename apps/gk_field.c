#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_array_average.h>

#include <assert.h>
#include <float.h>
#include <time.h>

// initialize field object
struct gk_field* 
gk_field_new(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app)
{
  struct gk_field *f = gkyl_malloc(sizeof(struct gk_field));

  f->info = gk->field;
  f->gkfield_id = f->info.gkfield_id ? f->info.gkfield_id : GKYL_GK_FIELD_ES;

  // We add the position of the lcfs in the poisson BC structure
  // to pass it to the deflated FEM solver
  f->info.poisson_bcs.xLCFS = f->info.xLCFS;

  // Add a updater to compute the flux surface average
  if (app->cdim == 3 && f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    gkyl_rect_grid_init(&f->grid_x, 1, &app->grid.lower[0], &app->grid.upper[0], &app->grid.cells[0]);
    gkyl_cart_modal_serendip(&f->confBasis_x, 1, app->poly_order);
    int ghost_x[] = {1};
    gkyl_create_grid_ranges(&f->grid_x, ghost_x, &f->local_x_ext, &f->local_x);
    int dim_fs_avg[] = {0,1,1};
    struct gkyl_array_average_inp input_fs_avg = {
      .grid = &app->grid,
      .basis = app->confBasis,
      .basis_avg = f->confBasis_x,
      .local = &app->local,
      .local_avg = &f->local_x,
      .local_avg_ext = &f->local_x_ext,
      .weight = app->gk_geom->jacobgeo,
      .avg_dim = dim_fs_avg,
      .use_gpu = app->use_gpu
    };
    f->up_fs_avg = gkyl_array_average_new(&input_fs_avg);
    f->fs_avg = mkarr(app->use_gpu, f->confBasis_x.num_basis, f->local_x_ext.volume);
    // gkyl_comm_allreduce
    // // Do a MPI sum reduction of the weight to have correct average when multi MPI proc
    // struct gkyl_array *Intw, *Intw_global;
    // // Get a pointer to the averaged weights
    // Intw = gkyl_array_average_acquire_weight_avg(f->up_fs_avg);
    // // Get the averaging volume
    // double avg_volume = gkyl_array_average_get_avg_volume(f->up_fs_avg);
    // // Multiply the averaged weights by the averaging volume
    // gkyl_array_scale_range(Intw, avg_volume, *f->local_x);
    // gkyl_array_accumulate()

    // gkyl_array_release(Intw);

  }
  // allocate an array to store the temperature (used in target corner bias)
  f->temp_elc = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  // allocate arrays for charge density
  f->rho_c = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  f->rho_c_global_dg = mkarr(app->use_gpu, app->confBasis.num_basis, app->global_ext.volume);
  f->rho_c_global_smooth = mkarr(app->use_gpu, app->confBasis.num_basis, app->global_ext.volume);

  // allocate arrays for electrostatic potential
  // global phi (only used in 1x simulations)
  f->phi_fem = mkarr(app->use_gpu, app->confBasis.num_basis, app->global_ext.volume);
  // local phi (assuming domain decomposition is *only* in z right now)
  f->phi_smooth = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    f->apar_fem = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    f->apardot_fem = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  }

  f->init_phi_pol = false;
  if (f->info.polarization_potential) {
    // Project the initial potential onto a p+1 tensor basis and compute the polarization
    // density to use use by species in calculating the initial ion density.
    f->init_phi_pol = true;
    struct gkyl_basis phi_pol_basis;
    gkyl_cart_modal_tensor(&phi_pol_basis, app->cdim, app->poly_order+1);

    f->phi_pol = mkarr(app->use_gpu, phi_pol_basis.num_basis, app->local_ext.volume);
    struct gkyl_array *phi_pol_ho = app->use_gpu? mkarr(false, f->phi_pol->ncomp, f->phi_pol->size)
                                            : gkyl_array_acquire(f->phi_pol);

    struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_new(&app->grid, &phi_pol_basis,
      1, f->info.polarization_potential, f->info.polarization_potential_ctx);
    gkyl_eval_on_nodes_advance(phi_pol_proj, 0.0, &app->local, phi_pol_ho);
    gkyl_array_copy(f->phi_pol, phi_pol_ho);
    
    gkyl_eval_on_nodes_release(phi_pol_proj);
    gkyl_array_release(phi_pol_ho);
  }

  // Create global subrange we'll copy the field solver solution from (into local).
  int intersect = gkyl_sub_range_intersect(&f->global_sub_range, &app->global, &app->local);

  // detect if this process contains an edge in the z dimension by comparing the local and global indices
  // for target corner BC
  int ndim = app->grid.ndim;
  f->info.poisson_bcs.contains_lower_z_edge = f->global_sub_range.lower[ndim-1] == app->global.lower[ndim-1];
  f->info.poisson_bcs.contains_upper_z_edge = f->global_sub_range.upper[ndim-1] == app->global.upper[ndim-1];

  if (f->gkfield_id == GKYL_GK_FIELD_BOLTZMANN || f->gkfield_id == GKYL_GK_FIELD_ADIABATIC)
    assert(app->cdim == 1); // Not yet implemented for cdim>1.

  f->weight = 0;
  f->epsilon = 0;
  f->kSq = 0;  // not currently used by fem_perp_poisson
  double polarization_weight = 0.0; 
  double es_energy_fac_1d_adiabatic = 0.0; 
  if (f->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
    polarization_weight = 1.0; 
    f->ambi_pot = gkyl_ambi_bolt_potential_new(&app->grid, &app->confBasis, 
      f->info.electron_mass, f->info.electron_charge, f->info.electron_temp, app->use_gpu);
    for (int j=0; j<app->cdim; ++j) {
      f->sheath_vals[2*j] = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
      f->sheath_vals[2*j+1] = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
    }
  } else {

    struct gkyl_array* bmag_mid_host = app->use_gpu? mkarr(false, 1, 1) : gkyl_array_acquire(app->gk_geom->bmag_mid);
    gkyl_array_copy(bmag_mid_host, app->gk_geom->bmag_mid);
    double *bmag_mid_ptr = gkyl_array_fetch(bmag_mid_host, 0);
    double polarization_bmag = f->info.polarization_bmag ? f->info.polarization_bmag : bmag_mid_ptr[0];
    gkyl_array_release(bmag_mid_host);
    // Linearized polarization density
    for (int i=0; i<app->num_species; ++i) {
      struct gk_species *s = &app->species[i];
      polarization_weight += s->info.polarization_density*s->info.mass/(polarization_bmag*polarization_bmag);
    }
    if (app->cdim == 1) {
      // Need to set weight to kperpsq*polarizationWeight for use in potential smoothing.
      f->weight = mkarr(false, app->confBasis.num_basis, app->global_ext.volume); // fem_parproj expects weight on host
 
      // Gather jacobgeo for smoothing in z.
      struct gkyl_array *jacobgeo_global = mkarr(app->use_gpu, app->confBasis.num_basis, app->global_ext.volume);
      gkyl_comm_array_allgather(app->comm, &app->local, &app->global, app->gk_geom->jacobgeo, jacobgeo_global);
      gkyl_array_copy(f->weight, jacobgeo_global);

      gkyl_array_scale(f->weight, polarization_weight);
      gkyl_array_scale(f->weight, f->info.kperpSq);


      if (f->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
        // Add the contribution from adiabatic electrons (in principle any
        // species can be adiabatic, which we can add suport for later).
        double n_s0 = f->info.electron_density;
        double q_s = f->info.electron_charge;
        double T_s = f->info.electron_temp;
        double quasineut_contr = q_s*n_s0*q_s/T_s;
        struct gkyl_array *weight_adiab = mkarr(false, app->confBasis.num_basis, app->global_ext.volume); // fem_parproj expects weight on host
        gkyl_array_copy(weight_adiab, jacobgeo_global);
        gkyl_array_scale(weight_adiab, quasineut_contr);
        gkyl_array_accumulate(f->weight, 1., weight_adiab);
        gkyl_array_release(weight_adiab);

        es_energy_fac_1d_adiabatic = 0.5*quasineut_contr; 
      }

      gkyl_array_release(jacobgeo_global);
    }
    else if (app->cdim > 1) {
      // set whatever epsilon we need
      // initialize a the weight to be used by deflated_fem_poisson
      f->epsilon = mkarr(app->use_gpu, (2*(app->cdim-1)-1)*app->confBasis.num_basis, app->local_ext.volume);
      gkyl_array_set_offset(f->epsilon, polarization_weight, app->gk_geom->gxxj, 0*app->confBasis.num_basis);
      if (app->cdim > 2) {
        gkyl_array_set_offset(f->epsilon, polarization_weight, app->gk_geom->gxyj, 1*app->confBasis.num_basis);
        gkyl_array_set_offset(f->epsilon, polarization_weight, app->gk_geom->gyyj, 2*app->confBasis.num_basis);
      }
      // deflated Poisson solve is performed on range assuming decomposition is *only* in z right now
      // need sub range of global range corresponding to where we are in z to properly index global charge density
      f->deflated_fem_poisson = gkyl_deflated_fem_poisson_new(app->grid, app->basis_on_dev.confBasis, app->confBasis,
        app->local, f->global_sub_range, f->epsilon, f->info.poisson_bcs, app->use_gpu);
    }
  }

  // Potential smoothing (in z) updater
  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    enum gkyl_fem_parproj_bc_type fem_parproj_bc_core, fem_parproj_bc_sol;
    if (app->cdim == 2) {
      fem_parproj_bc_core = GKYL_FEM_PARPROJ_PERIODIC;
      fem_parproj_bc_sol = GKYL_FEM_PARPROJ_NONE;
    } else {
      fem_parproj_bc_core = GKYL_FEM_PARPROJ_DIRICHLET;
      fem_parproj_bc_sol = GKYL_FEM_PARPROJ_NONE;
    }
    // construct core and SOL ranges.
    double xLCFS = f->info.xLCFS;
    // Index of the cell that abuts the xLCFS from below.
    int idxLCFS_m = (xLCFS-1e-8 - app->grid.lower[0])/app->grid.dx[0]+1;
    gkyl_range_shorten_from_below(&f->global_sol, &app->global, 0, app->grid.cells[0]-idxLCFS_m+1);
    gkyl_range_shorten_from_below(&f->global_ext_sol, &app->global_ext, 0, app->grid.cells[0]-idxLCFS_m+1);
    gkyl_range_shorten_from_above(&f->global_core, &app->global, 0, idxLCFS_m+1);
    gkyl_range_shorten_from_above(&f->global_ext_core, &app->global_ext, 0, idxLCFS_m+1);
    f->fem_parproj_core = gkyl_fem_parproj_new(&f->global_core, &f->global_ext_core, 
      &app->confBasis, fem_parproj_bc_core, f->weight, app->use_gpu);
    f->fem_parproj_sol = gkyl_fem_parproj_new(&f->global_sol, &f->global_ext_sol, 
      &app->confBasis, fem_parproj_bc_sol, f->weight, app->use_gpu);
  } 
  else {
    f->fem_parproj = gkyl_fem_parproj_new(&app->global, &app->global_ext, 
      &app->confBasis, f->info.fem_parbc, f->weight, app->use_gpu);
  }

  f->phi_host = f->phi_smooth;  
  if (app->use_gpu) {
    f->phi_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    f->em_energy_red = gkyl_cu_malloc(sizeof(double[1]));
    f->em_energy_red_global = gkyl_cu_malloc(sizeof(double[1]));
  } else {
    f->em_energy_red = gkyl_malloc(sizeof(double[1]));
    f->em_energy_red_global = gkyl_malloc(sizeof(double[1]));
  }

  f->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  f->is_first_energy_write_call = true;
  // Factors for ES energy. 
  f->es_energy_fac = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  f->es_energy_fac_1d = 0.0;
  if (app->cdim==1) {
    if (f->gkfield_id == GKYL_GK_FIELD_BOLTZMANN)
      f->es_energy_fac_1d = polarization_weight;
    else
      f->es_energy_fac_1d = polarization_weight*f->info.kperpSq + es_energy_fac_1d_adiabatic;

    f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->confBasis, 
      1, GKYL_ARRAY_INTEGRATE_OP_SQ, app->use_gpu);
  }
  else {
    gkyl_array_shiftc(f->es_energy_fac, sqrt(pow(2,app->cdim)), 0); // Sets es_energy_fac=1.
    gkyl_array_scale(f->es_energy_fac, 0.5*polarization_weight);

    f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->confBasis, 
      1, GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ, app->use_gpu);
  }

  // setup biased lower wall (same size as electrostatic potential), by default is 0.0
  f->phi_wall_lo = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  f->has_phi_wall_lo = false;
  f->phi_wall_lo_evolve = false;
  if (f->info.phi_wall_lo) {
    f->has_phi_wall_lo = true;
    if (f->info.phi_wall_lo_evolve)
      f->phi_wall_lo_evolve = f->info.phi_wall_lo_evolve;

    f->phi_wall_lo_host = f->phi_wall_lo;
    if (app->use_gpu) 
      f->phi_wall_lo_host = mkarr(false, f->phi_wall_lo->ncomp, f->phi_wall_lo->size);

    f->phi_wall_lo_proj = gkyl_eval_on_nodes_new(&app->grid, &app->confBasis, 
      1, f->info.phi_wall_lo, f->info.phi_wall_lo_ctx);

    // Compute phi_wall_lo at t = 0
    gkyl_eval_on_nodes_advance(f->phi_wall_lo_proj, 0.0, &app->local_ext, f->phi_wall_lo_host);
    if (app->use_gpu) // note: phi_wall_lo_host is same as phi_wall_lo when not on GPUs
      gkyl_array_copy(f->phi_wall_lo, f->phi_wall_lo_host);
  }

  // setup biased upper wall (same size as electrostatic potential), by default is 0.0
  f->phi_wall_up = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  f->has_phi_wall_up = false;
  f->phi_wall_up_evolve = false;
  if (f->info.phi_wall_up) {
    f->has_phi_wall_up = true;
    if (f->info.phi_wall_up_evolve)
      f->phi_wall_up_evolve = f->info.phi_wall_up_evolve;

    f->phi_wall_up_host = f->phi_wall_up;
    if (app->use_gpu) 
      f->phi_wall_up_host = mkarr(false, f->phi_wall_up->ncomp, f->phi_wall_up->size);

    f->phi_wall_up_proj = gkyl_eval_on_nodes_new(&app->grid, &app->confBasis,
      1, f->info.phi_wall_up, f->info.phi_wall_up_ctx);

    // Compute phi_wall_up at t = 0
    gkyl_eval_on_nodes_advance(f->phi_wall_up_proj, 0.0, &app->local_ext, f->phi_wall_up_host);
    if (app->use_gpu) // note: phi_wall_up_host is same as phi_wall_up when not on GPUs
      gkyl_array_copy(f->phi_wall_up, f->phi_wall_up_host);
  }

  // If we are in a 3x simulation and IWL, we add TS BC and SSFG updaters
  if(app->cdim ==3){
    if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
      gk_field_add_TSBC_and_SSFG_updaters(app,f);
    }
  }
  
  return f;
}

void
gk_field_add_TSBC_and_SSFG_updaters(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  // We take the first species to copy the function for the TS BC
  struct gk_species *gks = &app->species[0];
  // Get the z BC info from the first species in our app
  const struct gkyl_gyrokinetic_bcs *bcz = &gks->info.bcz;
  // we consider only dir=2 (z-BC)
  int zdir = 2; 

  // Create local and local_ext from app local range with ghosts
  int ghost[] = {1, 1, 1};
  gkyl_create_ranges(&app->local, ghost, &f->local_ext, &f->local);

  //-------- Define sub range of ghost and skin cells that spans only the core
  double xLCFS = f->info.xLCFS;
  // Index of the cell that abuts the xLCFS from below.
  int idxLCFS_m = (xLCFS-1e-8 - app->grid.lower[0])/app->grid.dx[0]+1;

  // Create a core local range, extended in the BC dir.
  int ndim = app->cdim;
  int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
  for (int i=0; i<ndim; i++) {
    lower_bcdir_ext[i] = f->local.lower[i];
    upper_bcdir_ext[i] = f->local.upper[i];
  }
  upper_bcdir_ext[0] = idxLCFS_m;
  lower_bcdir_ext[zdir] = f->local_ext.lower[zdir];
  upper_bcdir_ext[zdir] = f->local_ext.upper[zdir];
  gkyl_sub_range_init(&f->local_par_ext_core, &f->local_ext, lower_bcdir_ext, upper_bcdir_ext);

  //TS BC updater for up to low TS for the lower edge
  //this sets ghost_L = T_LU(ghost_L)
  struct gkyl_bc_twistshift_inp T_LU_lo = {
    .bc_dir = zdir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = GKYL_LOWER_EDGE,
    .cdim = app->cdim,
    .bcdir_ext_update_r = f->local_par_ext_core,
    .num_ghost = ghost, // one ghost per config direction
    .basis = app->confBasis,
    .grid = app->grid,
    .shift_func = bcz->lower.aux_profile,
    .shift_func_ctx = bcz->lower.aux_ctx,
    .use_gpu = app->use_gpu,
  };

  //TS BC updater for up to low TS for the lower edge
  //this sets ghost_L = T_UL(ghost_L)
  struct gkyl_bc_twistshift_inp T_UL_lo = {
    .bc_dir = zdir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = GKYL_LOWER_EDGE,
    .cdim = app->cdim,
    .bcdir_ext_update_r = f->local_par_ext_core,
    .num_ghost = ghost,
    .basis = app->confBasis,
    .grid = app->grid,
    .shift_func = bcz->upper.aux_profile,
    .shift_func_ctx = bcz->upper.aux_ctx,
    .use_gpu = app->use_gpu,
  };

  //TS BC updater for low to up TS for the upper edge
  //this sets ghost_U = T_UL(ghost_U)
  struct gkyl_bc_twistshift_inp T_UL_up = {
    .bc_dir = zdir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = GKYL_UPPER_EDGE,
    .cdim = app->cdim,
    .bcdir_ext_update_r = f->local_par_ext_core,
    .num_ghost = ghost,
    .basis = app->confBasis,
    .grid = app->grid,
    .shift_func = bcz->upper.aux_profile,
    .shift_func_ctx = bcz->upper.aux_ctx,
    .use_gpu = app->use_gpu,
  };

  //TS BC updater for up to low TS for the lower edge
  //this sets ghost_U = T_LU(ghost_U)
  struct gkyl_bc_twistshift_inp T_LU_up = {
    .bc_dir = zdir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = GKYL_UPPER_EDGE,
    .cdim = app->cdim,
    .bcdir_ext_update_r = f->local_par_ext_core,
    .num_ghost = ghost,
    .basis = app->confBasis,
    .grid = app->grid,
    .shift_func = bcz->lower.aux_profile,
    .shift_func_ctx = bcz->lower.aux_ctx,
    .use_gpu = app->use_gpu,
  };

  // Add the forward TS updater to f
  f->bc_T_LU_lo = gkyl_bc_twistshift_new(&T_LU_lo);
  f->bc_T_UL_up = gkyl_bc_twistshift_new(&T_UL_up);
  // Add the backward TS updater to f
  f->bc_T_UL_lo = gkyl_bc_twistshift_new(&T_UL_lo);
  f->bc_T_LU_up = gkyl_bc_twistshift_new(&T_LU_up);

  //------------ Skin surface from ghost updater
  // The par_ext ranges have only extension of ghosts in the parallel direciton, 
  // hence their ghost array is not {1,1,1} but
  int ghost_par[] = {0, 0, 1};
  // create lower and upper skin and ghost ranges for the z BC in the core region
  gkyl_skin_ghost_ranges( &f->lower_skin_core, &f->lower_ghost_core, zdir, 
                          GKYL_LOWER_EDGE, &f->local_par_ext_core, ghost_par);
  gkyl_skin_ghost_ranges( &f->upper_skin_core, &f->upper_ghost_core, zdir, 
                          GKYL_UPPER_EDGE, &f->local_par_ext_core, ghost_par);
  // add the SSFG updater for lower and upper application
  f->ssfg_lo = gkyl_skin_surf_from_ghost_new(zdir,GKYL_LOWER_EDGE,
                app->confBasis,&f->lower_skin_core,&f->lower_ghost_core,app->use_gpu);
  f->ssfg_up = gkyl_skin_surf_from_ghost_new(zdir,GKYL_UPPER_EDGE,
                app->confBasis,&f->upper_skin_core,&f->upper_ghost_core,app->use_gpu);

  //We allocate the array to store the ghost values of phi_smooth but switching upper and lower
  f->aux_array = mkarr(app->use_gpu, f->phi_smooth->ncomp, f->phi_smooth->size);

  // Create a BC_REFLECT updater.
  f->bc_reflect_lo = gkyl_bc_basic_new(zdir, GKYL_LOWER_EDGE, GKYL_BC_REFLECT, app->basis_on_dev.confBasis,
      &f->lower_skin_core, &f->lower_ghost_core, f->aux_array->ncomp, app->cdim, app->use_gpu);
  f->bc_reflect_up = gkyl_bc_basic_new(zdir, GKYL_UPPER_EDGE, GKYL_BC_REFLECT, app->basis_on_dev.confBasis,
      &f->upper_skin_core, &f->upper_ghost_core, f->aux_array->ncomp, app->cdim, app->use_gpu);
  // These bc_basic updaters need a buffer, allocate it.
  long buff_sz = GKYL_MAX2(f->lower_skin_core.volume, f->upper_skin_core.volume);
  f->bc_buffer = mkarr(app->use_gpu, app->confBasis.num_basis, buff_sz);

  // Finally we need a config space communicator to sync the inner cell data and
  // avoid applying TS BC inside the domain
  f->comm_conf = gkyl_comm_extend_comm(app->comm, &f->local);
}

void
gk_field_calc_phi_wall(gkyl_gyrokinetic_app *app, struct gk_field *field, double tm)
{
  if (field->has_phi_wall_lo && field->phi_wall_lo_evolve) {
    gkyl_eval_on_nodes_advance(field->phi_wall_lo_proj, tm, &app->local_ext, field->phi_wall_lo_host);
    if (app->use_gpu) // note: phi_wall_lo_host is same as phi_wall_lo when not on GPUs
      gkyl_array_copy(field->phi_wall_lo, field->phi_wall_lo_host);
  }
  if (field->has_phi_wall_up && field->phi_wall_up_evolve) {
    gkyl_eval_on_nodes_advance(field->phi_wall_up_proj, tm, &app->local_ext, field->phi_wall_up_host);
    if (app->use_gpu) // note: phi_wall_up_host is same as phi_wall_up when not on GPUs
      gkyl_array_copy(field->phi_wall_up, field->phi_wall_up_host);
  }
}

void
gk_field_accumulate_rho_c(gkyl_gyrokinetic_app *app, struct gk_field *field, 
  const struct gkyl_array *fin[])
{
  gkyl_array_clear(field->rho_c, 0.0);
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];

    gk_species_moment_calc(&s->m0, s->local, app->local, fin[i]);
    if (field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
      // For Boltzmann electrons, we only need ion density, not charge density.
      gkyl_array_accumulate_range(field->rho_c, 1.0, s->m0.marr, &app->local);
    } else {
      gkyl_array_accumulate_range(field->rho_c, s->info.charge, s->m0.marr, &app->local);
      if (field->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
        // Add the background (electron) charge density.
        double n_s0 = field->info.electron_density;
        double q_s = field->info.electron_charge;
        gkyl_array_shiftc_range(field->rho_c, q_s*n_s0*sqrt(2.), 0, &app->local);
      }
    }
  } 
}

void
gk_field_calc_ambi_pot_sheath_vals(gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];

    // Assumes symmetric sheath BCs for now only in 1D
    gkyl_ambi_bolt_potential_sheath_calc(field->ambi_pot, GKYL_LOWER_EDGE, 
      &app->lower_skin[0], &app->lower_ghost[0], app->gk_geom->jacobgeo_inv, 
      s->bflux.gammai[0].marr, field->rho_c, field->sheath_vals[0]);
    gkyl_ambi_bolt_potential_sheath_calc(field->ambi_pot, GKYL_UPPER_EDGE, 
      &app->upper_skin[0], &app->upper_ghost[0], app->gk_geom->jacobgeo_inv, 
      s->bflux.gammai[1].marr, field->rho_c, field->sheath_vals[1]);

    // Broadcast the sheath values from skin processes to other processes.
    int comm_sz[1];
    gkyl_comm_get_size(app->comm, comm_sz);
    gkyl_comm_array_bcast(app->comm, field->sheath_vals[0], field->sheath_vals[0], 0);
    gkyl_comm_array_bcast(app->comm, field->sheath_vals[1], field->sheath_vals[1], comm_sz[0]-1);
    gkyl_array_accumulate(field->sheath_vals[0], 1., field->sheath_vals[1]);
    gkyl_array_scale(field->sheath_vals[0], 0.5);
  } 
}

void 
gk_field_calc_target_corner_bias(gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *fin[]){

  bool constant_TCBC = true;

  if (constant_TCBC){ // phi = 0 time constant target corner bias
    field->target_corner_bias = 0;

  } else { // time dependent target corner bias
    // -- This is the method for setting phi_TC = phi_fs(xLCFS)
    // We update the flux surface average of phi as well
    // gkyl_array_average_advance(field->up_fs_avg, field->phi_smooth, field->fs_avg);

    // -- This is the method for setting phi_TC = lambda Te_fs(xLCFS)/e
    // find the electrons
    struct gk_species *s_elc, *s_ion;
    int elc_idx = 0;
    int ion_idx = 0;
    for (int i=0; i<app->num_species; ++i){
      if (strcmp("elc", app->species[i].info.name) == 0){
        s_elc = &app->species[i];
        elc_idx = i;
      }
      if (strcmp("ion", app->species[i].info.name) == 0){
        s_ion = &app->species[i];
        ion_idx = i;
      }
    }
    gk_species_moment_calc(&s_elc->maxwellian_moments, s_elc->local, app->local, fin[elc_idx]);
    gkyl_array_set_offset(field->temp_elc, s_elc->info.mass, s_elc->maxwellian_moments.marr, 2*app->confBasis.num_basis);
    gkyl_array_average_advance(field->up_fs_avg, field->temp_elc, field->fs_avg);

    // get the flux-surface averaged temp at the LCFS interface
    // (for now we take the average of the cell right to the LCFS)
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &field->local_x);
    // Index of the cell that abuts the xLCFS from below.
    int idxLCFS_m = (field->info.xLCFS-1e-8 - app->grid.lower[0])/app->grid.dx[0]+1;
    // Get the avg value on host
    while (gkyl_range_iter_next(&iter)) {
      long lidx_avg = gkyl_range_idx(&field->local_x, iter.idx);
      if (iter.idx[0] == idxLCFS_m) {
        const double *avg_value = gkyl_array_cfetch(field->fs_avg, lidx_avg);
        if (app->use_gpu)
          gkyl_cu_memcpy(&field->fs_avg_LCFS, avg_value, sizeof(double), GKYL_CU_MEMCPY_D2H);
        else
          memcpy(&field->fs_avg_LCFS, avg_value, sizeof(double));
        // Multiply by sqrt(2)/2 to get the average value of the cell from the 0th DG coeffitient
        field->fs_avg_LCFS = sqrt(2.0)/2.0 * field->fs_avg_LCFS;
      }
    }
    // Reduce the avg over all MPI processes
    int err_check = gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &field->fs_avg_LCFS, &field->fs_LCFS_global);
    assert(err_check == false);
    int nproc;
    err_check = gkyl_comm_get_size(app->comm, &nproc);
    assert(err_check == false);
    field->fs_LCFS_global /= nproc;

    // fs avg phi version
    // field->target_corner_bias = field->fs_LCFS_global;

    // For the temperature method we compute now lambda Te/e
    // double me = s_elc->info.mass;
    // double qe = s_elc->info.charge;
    // double mi = s_ion->info.mass;
    // double qi = s_ion->info.charge;
    // double Ti0 = 100*qe;//100 eV (this is hard encoded)
    // double Te0 = 100*qe;//100 eV
    // double lambda = sqrt(log(qi*mi/(qe*me) / (2.*M_PI*(1. + Ti0/Te0))));
    // We hard encode it for now
    double lambda = 2.2325777683556205;
    field->target_corner_bias = lambda * field->fs_LCFS_global / 1.602176634e-19; // this is lambda Te / e
  }
}

// Compute the electrostatic potential
void
gk_field_rhs(gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  struct timespec wst = gkyl_wall_clock();
  if (field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) { 
    gkyl_ambi_bolt_potential_phi_calc(field->ambi_pot, &app->local, &app->local_ext,
      app->gk_geom->jacobgeo_inv, field->rho_c, field->sheath_vals[0], field->phi_smooth);

    // gather potential into global array for smoothing in z
    gkyl_comm_array_allgather(app->comm, &app->local, &app->global, field->phi_smooth, field->rho_c_global_dg);

    gkyl_fem_parproj_set_rhs(field->fem_parproj, field->rho_c_global_dg, field->rho_c_global_dg);
    gkyl_fem_parproj_solve(field->fem_parproj, field->phi_fem);
    // copy globally smoothed potential to local potential per process for update
    gkyl_array_copy_range_to_range(field->phi_smooth, field->phi_fem, &app->local, &field->global_sub_range);
  }
  else {
    // gather charge density into global array for smoothing in z
    gkyl_comm_array_allgather(app->comm, &app->local, &app->global, field->rho_c, field->rho_c_global_dg);
    if (app->cdim == 1) {
      gkyl_fem_parproj_set_rhs(field->fem_parproj, field->rho_c_global_dg, field->rho_c_global_dg);
      gkyl_fem_parproj_solve(field->fem_parproj, field->phi_fem);
      // copy globally smoothed potential to local potential per process for update
      gkyl_array_copy_range_to_range(field->phi_smooth, field->phi_fem, &app->local, &field->global_sub_range);
    }
    else if (app->cdim > 1) {
      // input is rho_c_global_dg, globally smoothed in z, and then output should be in *local* phi_smooth
      if (field->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
        gkyl_fem_parproj_set_rhs(field->fem_parproj_core, field->rho_c_global_dg, field->rho_c_global_dg);
        gkyl_fem_parproj_solve(field->fem_parproj_core, field->rho_c_global_smooth);
        gkyl_fem_parproj_set_rhs(field->fem_parproj_sol, field->rho_c_global_dg, field->rho_c_global_dg);
        gkyl_fem_parproj_solve(field->fem_parproj_sol, field->rho_c_global_smooth);
      }
      else {
        gkyl_fem_parproj_set_rhs(field->fem_parproj, field->rho_c_global_dg, field->rho_c_global_dg);
        gkyl_fem_parproj_solve(field->fem_parproj, field->rho_c_global_smooth);
      }

      // printf("target corner b ias = %g\n", field->target_corner_bias[0]);
      gk_field_calc_target_corner_bias(app, field, NULL);
      gkyl_deflated_fem_poisson_advance(field->deflated_fem_poisson, field->rho_c_global_smooth, field->phi_smooth, field->target_corner_bias);

      /*
      * If we are in a 3x simulation with IWL we apply TS BC to the upper and lower edges
      * so that we can enforce that the surface value of the global skin cells are matching the
      * twist-and-shift BC exactly. (this should enforce periodicity of y-avg phi)
      */
      if (field->gkfield_id == GKYL_GK_FIELD_ES_IWL && app->cdim == 3)   
        gk_field_apply_bc(app, field, field->phi_smooth);

    }
  }
  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_apply_bc(const gkyl_gyrokinetic_app *app, const struct gk_field *field, struct gkyl_array *finout)
{
  //1. Apply the periodicity to fill the ghost cells
  int num_periodic_dir = 1; // we need only periodicity in z
  int periodic_dirs[] = {2}; // z direction
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext,
    num_periodic_dir, periodic_dirs, finout); 
  // finout ghost cells: | f_up | ----- | f_lo |

  //=== We now perform the phase 1a, i.e. 
  // f_L = T_LU(f_up)
  // f_U = f_up

  //2. call the TS BC updater to update the z ghosts with TS
  gkyl_bc_twistshift_advance(field->bc_T_LU_lo, finout, finout);
  // finout is now | T_LU(f_up) | ----- | f_lo |

  //3. Synchronize the array between the MPI processes
  gkyl_comm_array_sync(field->comm_conf, &field->local, &field->local_ext, finout);
/* Note:
  * Here the TS BC has been applied blindly i.e. regardless if the process is at the 
  * border of the domain. Technically, only the processes with the global skin cells 
  * should apply it. We resolve this with an array sync so that the ghost
  * cells of inner skin cells are overwritten by the neighbor process skin value.
  * This method creates a lighter code and is the same to the one used for the TS BC
  * of the distribution function in gk_species.c
  */

  //4. Copy the ghost surface value to the skin surface value (SSFG)
  gkyl_skin_surf_from_ghost_advance(field->ssfg_lo, finout);
  // finout is now | T_LU(f_up) | \---- | f_lo |
  // "\" denote that the skin cell has been adapted to the edge (only the lower)
  // The upper skin cell remains unchanged
}

void
gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *field)
{
  gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth, 
    app->grid.cellVolume, field->es_energy_fac, &app->local, field->em_energy_red);

  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red, field->em_energy_red_global);

  double energy_global[1] = { 0.0 };
  if (app->use_gpu)
    gkyl_cu_memcpy(energy_global, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  else
    energy_global[0] = field->em_energy_red_global[0];
  
  if (app->cdim == 1)
    energy_global[0] *= field->es_energy_fac_1d; 

  gkyl_dynvec_append(field->integ_energy, tm, energy_global);
}

// release resources for field
void
gk_field_release(const gkyl_gyrokinetic_app* app, struct gk_field *f)
{
  gkyl_array_release(f->rho_c);
  gkyl_array_release(f->rho_c_global_dg);
  gkyl_array_release(f->rho_c_global_smooth);
  gkyl_array_release(f->phi_fem);
  gkyl_array_release(f->phi_smooth);

  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    gkyl_array_release(f->apar_fem);
    gkyl_array_release(f->apardot_fem);
  }

  if (f->init_phi_pol) {
    gkyl_array_release(f->phi_pol);
  }

  gkyl_array_release(f->es_energy_fac);
  if (f->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
    gkyl_ambi_bolt_potential_release(f->ambi_pot);
    for (int i=0; i<2*app->cdim; ++i) 
      gkyl_array_release(f->sheath_vals[i]);
  } 
  else {
    if (app->cdim == 1) {
      gkyl_array_release(f->weight);
    }
    else if (app->cdim > 1) {
      gkyl_array_release(f->epsilon);
      gkyl_deflated_fem_poisson_release(f->deflated_fem_poisson);
    }
  }

  if (app->cdim > 1 && f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    gkyl_fem_parproj_release(f->fem_parproj_core);
    gkyl_fem_parproj_release(f->fem_parproj_sol);
  } else {
    gkyl_fem_parproj_release(f->fem_parproj);
  }

  // Release TS BS and SSFG updater
  if (app->cdim == 3 && f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    gkyl_bc_twistshift_release(f->bc_T_UL_up);
    gkyl_bc_twistshift_release(f->bc_T_UL_lo);
    gkyl_bc_twistshift_release(f->bc_T_LU_up);
    gkyl_bc_twistshift_release(f->bc_T_LU_lo);
    gkyl_skin_surf_from_ghost_release(f->ssfg_up);
    gkyl_skin_surf_from_ghost_release(f->ssfg_lo);
    gkyl_comm_release(f->comm_conf);
    gkyl_array_release(f->aux_array);
    gkyl_array_release(f->bc_buffer);
    gkyl_bc_basic_release(f->bc_reflect_lo);
    gkyl_bc_basic_release(f->bc_reflect_up);
    gkyl_array_average_release(f->up_fs_avg);
    gkyl_array_release(f->fs_avg);
    gkyl_array_release(f->temp_elc);
  }

  gkyl_dynvec_release(f->integ_energy);
  gkyl_array_integrate_release(f->calc_em_energy);
  if (app->use_gpu) {
    gkyl_array_release(f->phi_host);
    gkyl_cu_free(f->em_energy_red);
    gkyl_cu_free(f->em_energy_red_global);
  } else {
    gkyl_free(f->em_energy_red);
    gkyl_free(f->em_energy_red_global);
  }

  gkyl_array_release(f->phi_wall_lo);
  gkyl_array_release(f->phi_wall_up);
  if (f->has_phi_wall_lo) {
    gkyl_eval_on_nodes_release(f->phi_wall_lo_proj);
    if (app->use_gpu) 
      gkyl_array_release(f->phi_wall_lo_host);
  }
  if (f->has_phi_wall_up) {
    gkyl_eval_on_nodes_release(f->phi_wall_up_proj);
    if (app->use_gpu) 
      gkyl_array_release(f->phi_wall_up_host);
  }

  gkyl_free(f);
}

