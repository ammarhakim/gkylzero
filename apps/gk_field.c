#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_comm_io.h>

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

  f->calc_init_field = !f->info.zero_init_field;
  f->update_field = !f->info.is_static;
  // The combination update_field=true, calc_init_field=false is not allowed.
  assert(!(f->update_field && (!f->calc_init_field)));

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

  if (f->gkfield_id == GKYL_GK_FIELD_BOLTZMANN || f->gkfield_id == GKYL_GK_FIELD_ADIABATIC)
    assert(app->cdim == 1); // Not yet implemented for cdim>1.

  f->epsilon = 0;
  struct gkyl_array *epsilon_global = 0;
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

    // Allocate array for the polarization weight times geometric coefficients.
    f->epsilon = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->confBasis.num_basis, app->local_ext.volume);

    double polarization_bmag = f->info.polarization_bmag ? f->info.polarization_bmag : app->bmag_ref;
    // Linearized polarization density
    for (int i=0; i<app->num_species; ++i) {
      struct gk_species *s = &app->species[i];
      polarization_weight += s->info.polarization_density*s->info.mass/pow(polarization_bmag,2);
    }
    if (app->cdim == 1) {
      // Need to set weight to kperpsq*polarizationWeight for use in potential smoothing.
      gkyl_array_copy(f->epsilon, app->gk_geom->jacobgeo);
      gkyl_array_scale(f->epsilon, polarization_weight);
      gkyl_array_scale(f->epsilon, f->info.kperpSq);


      if (f->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
        // Add the contribution from adiabatic electrons (in principle any
        // species can be adiabatic, which we can add suport for later).
        double n_s0 = f->info.electron_density;
        double q_s = f->info.electron_charge;
        double T_s = f->info.electron_temp;
        double quasineut_contr = q_s*n_s0*q_s/T_s;
        
        struct gkyl_array *epsilon_adiab = mkarr(app->use_gpu, f->epsilon->ncomp, f->epsilon->size);
        gkyl_array_copy(epsilon_adiab, app->gk_geom->jacobgeo);
        gkyl_array_scale(epsilon_adiab, quasineut_contr);
        gkyl_array_accumulate(f->epsilon, 1., epsilon_adiab);
        gkyl_array_release(epsilon_adiab);

        es_energy_fac_1d_adiabatic = 0.5*quasineut_contr; 
      }

      // Gather gather epsilon for (global) smoothing in z.
      epsilon_global = mkarr(app->use_gpu, f->epsilon->ncomp, app->global_ext.volume);
      gkyl_comm_array_allgather(app->comm, &app->local, &app->global, f->epsilon, epsilon_global);
    }
    else if (app->cdim > 1) {
      // Initialize the polarization weight.
      gkyl_array_set_offset(f->epsilon, polarization_weight, app->gk_geom->gxxj, 0*app->confBasis.num_basis);
      if (app->cdim > 2) {
        gkyl_array_set_offset(f->epsilon, polarization_weight, app->gk_geom->gxyj, 1*app->confBasis.num_basis);
        gkyl_array_set_offset(f->epsilon, polarization_weight, app->gk_geom->gyyj, 2*app->confBasis.num_basis);
      }
      // Initialize the Poisson solver.
//      f->deflated_fem_poisson = gkyl_deflated_fem_poisson_new(app->grid, app->basis_on_dev.confBasis, app->confBasis,
//        app->local, f->global_sub_range, f->epsilon, f->info.poisson_bcs, app->use_gpu);
      f->fem_poisson = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->confBasis,
        &f->info.poisson_bcs, f->epsilon, NULL, app->use_gpu);
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

    f->fem_parproj_core = gkyl_fem_parproj_new(&f->global_core, &app->confBasis,
      fem_parproj_bc_core, 0, 0, app->use_gpu);
    f->fem_parproj_sol = gkyl_fem_parproj_new(&f->global_sol, &app->confBasis,
      fem_parproj_bc_sol, 0, 0, app->use_gpu);
  } 
  else {
    f->fem_parproj = gkyl_fem_parproj_new(&app->global, &app->confBasis,
      f->info.fem_parbc, epsilon_global, 0, app->use_gpu);
  }

  if (epsilon_global)
    gkyl_array_release(epsilon_global);

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
  f->es_energy_fac = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->confBasis.num_basis, app->local_ext.volume);
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
    gkyl_array_set(f->es_energy_fac, 0.5, f->epsilon);

    f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->confBasis, 
      1, GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ, app->use_gpu);
  }

  if (app->fdot_diagnostics) {
    if (app->use_gpu) {
      f->em_energy_red_new = gkyl_cu_malloc(sizeof(double[1]));
      f->em_energy_red_old = gkyl_cu_malloc(sizeof(double[1]));
    } else {
      f->em_energy_red_new = gkyl_malloc(sizeof(double[1]));
      f->em_energy_red_old = gkyl_malloc(sizeof(double[1]));
    }
    f->integ_energy_dot = gkyl_dynvec_new(GKYL_DOUBLE, 1);
    f->is_first_energy_dot_write_call = true;
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

  return f;
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
      s->bflux_solver.gammai[0].marr, field->rho_c, field->sheath_vals[0]);
    gkyl_ambi_bolt_potential_sheath_calc(field->ambi_pot, GKYL_UPPER_EDGE, 
      &app->upper_skin[0], &app->upper_ghost[0], app->gk_geom->jacobgeo_inv, 
      s->bflux_solver.gammai[1].marr, field->rho_c, field->sheath_vals[1]);

    // Broadcast the sheath values from skin processes to other processes.
    int comm_sz[1];
    gkyl_comm_get_size(app->comm, comm_sz);
    gkyl_comm_array_bcast(app->comm, field->sheath_vals[0], field->sheath_vals[0], 0);
    gkyl_comm_array_bcast(app->comm, field->sheath_vals[1], field->sheath_vals[1], comm_sz[0]-1);
    gkyl_array_accumulate(field->sheath_vals[0], 1., field->sheath_vals[1]);
    gkyl_array_scale(field->sheath_vals[0], 0.5);
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
      // Smooth the charge density. Input is rho_c_global_dg, globally smoothed in z,
      // and then output should be in *local* phi_smooth.
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
//      gkyl_deflated_fem_poisson_advance(field->deflated_fem_poisson, field->rho_c_global_smooth, field->phi_smooth);

      // Solve the Poisson equation.
      gkyl_fem_poisson_perp_set_rhs(field->fem_poisson, field->rho_c_global_smooth);
      gkyl_fem_poisson_perp_solve(field->fem_poisson, field->phi_smooth);

      // Smooth the potential.
      gkyl_fem_parproj_set_rhs(field->fem_parproj, field->phi_smooth, field->phi_smooth);
      gkyl_fem_parproj_solve(field->fem_parproj, field->phi_smooth);
    }
  }
  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);
}

static struct gkyl_app_restart_status
header_from_file(gkyl_gyrokinetic_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = { .io_status = 2 };
  
  FILE *fp = 0;
  with_file(fp, fname, "r") {
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    rstat.io_status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);

    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      if (hdr.etype != GKYL_DOUBLE)
        rstat.io_status = GKYL_ARRAY_RIO_DATA_MISMATCH;
    }

    struct gyrokinetic_output_meta meta =
      gk_meta_from_mpack( &(struct gkyl_msgpack_data) {
          .meta = hdr.meta,
          .meta_sz = hdr.meta_size
        }
      );

    rstat.frame = meta.frame;
    rstat.stime = meta.stime;

    gkyl_grid_sub_array_header_release(&hdr);
  }
  
  return rstat;
}

void
gk_field_file_import_init(struct gkyl_gyrokinetic_app *app, struct gkyl_gyrokinetic_ic_import inp)
{
  // Import initial condition from a file.
  struct gkyl_app_restart_status rstat = header_from_file(app, inp.file_name);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    struct gkyl_app_restart_status rstat;
    rstat.io_status = gkyl_comm_array_read(app->comm, &app->grid, &app->local, app->field->phi_host, inp.file_name);
    gkyl_array_copy(app->field->phi_smooth, app->field->phi_host);
  }
  else
    assert(false);
}

void
gk_field_project_init(struct gkyl_gyrokinetic_app *app)
{
  // Project the initial field.
  struct gkyl_eval_on_nodes *phi_proj = gkyl_eval_on_nodes_new(&app->grid, &app->confBasis,
    1, app->field->info.init_field_profile, app->field->info.init_field_profile_ctx);
  gkyl_eval_on_nodes_advance(phi_proj, 0.0, &app->local, app->field->phi_host);
  gkyl_eval_on_nodes_release(phi_proj);
  gkyl_array_copy(app->field->phi_smooth, app->field->phi_host);
}

void
gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *field)
{
  gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth, 
    1.0, field->es_energy_fac, &app->local, &app->local, field->em_energy_red);

  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red, field->em_energy_red_global);

  double energy_global[1] = { 0.0 };
  if (app->use_gpu)
    gkyl_cu_memcpy(energy_global, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  else
    energy_global[0] = field->em_energy_red_global[0];
  
  if (app->cdim == 1)
    energy_global[0] *= field->es_energy_fac_1d; 

  gkyl_dynvec_append(field->integ_energy, tm, energy_global);

  if (app->fdot_diagnostics) {
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red_old, field->em_energy_red_global);
    double energy_dot_global_old[1] = { 0.0 };
    if (app->use_gpu)
      gkyl_cu_memcpy(energy_dot_global_old, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    else
      energy_dot_global_old[0] = field->em_energy_red_global[0];
    if (app->cdim == 1)
      energy_dot_global_old[0] *= field->es_energy_fac_1d; 

    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red_new, field->em_energy_red_global);
    double energy_dot_global_new[1] = { 0.0 };
    if (app->use_gpu)
      gkyl_cu_memcpy(energy_dot_global_new, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    else
      energy_dot_global_new[0] = field->em_energy_red_global[0];
    if (app->cdim == 1)
      energy_dot_global_new[0] *= field->es_energy_fac_1d; 

    double energy_dot_global[1] = { 0.0 };
    energy_dot_global[0] = energy_dot_global_new[0] - energy_dot_global_old[0];

    gkyl_dynvec_append(field->integ_energy_dot, tm, energy_dot_global);
  }
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
    gkyl_array_release(f->epsilon);
    if (app->cdim > 1) {
//      gkyl_deflated_fem_poisson_release(f->deflated_fem_poisson);
      gkyl_fem_poisson_perp_release(f->fem_poisson);
    }
  }

  if (app->cdim > 1 && f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    gkyl_fem_parproj_release(f->fem_parproj_core);
    gkyl_fem_parproj_release(f->fem_parproj_sol);
  } else {
    gkyl_fem_parproj_release(f->fem_parproj);
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

  if (app->fdot_diagnostics) {
    gkyl_dynvec_release(f->integ_energy_dot);
    if (app->use_gpu) {
      gkyl_cu_free(f->em_energy_red_new);
      gkyl_cu_free(f->em_energy_red_old);
    } else {
      gkyl_free(f->em_energy_red_new);
      gkyl_free(f->em_energy_red_old);
    }
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

