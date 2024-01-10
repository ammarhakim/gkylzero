#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>

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

  // allocate arrays for charge density
  f->rho_c = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  f->rho_c_smooth = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  // allocate arrays for electrostatic potential
  f->phi_fem = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  f->phi_smooth = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    f->apar_fem = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    f->apardot_fem = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  }

  f->weight = 0;
  f->epsilon = 0;
  f->kSq = 0;  // not currently used by fem_perp_poisson
  double polarization_weight = 0.0; 
  if (f->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
    polarization_weight = 1.0; 
    f->ambi_pot = gkyl_ambi_bolt_potential_new(&app->grid, &app->confBasis, 
      f->info.electron_mass, f->info.electron_charge, f->info.electron_temp, app->use_gpu);
    for (int j=0; j<app->cdim; ++j) {
      f->sheath_vals[2*j] = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
      f->sheath_vals[2*j+1] = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
    }
  }
  else {
    // Linearized polarization density
    for (int i=0; i<app->num_species; ++i) {
      struct gk_species *s = &app->species[i];
      polarization_weight += s->info.polarization_density*s->info.mass/(f->info.bmag_fac*f->info.bmag_fac);
    }
    if (app->cdim == 1) {
      // in 1D case need to set weight to kperpsq*polarizationWeight for use in potential smoothing
      f->weight = mkarr(false, app->confBasis.num_basis, app->local_ext.volume); // fem_parproj expects weight on host
      gkyl_array_shiftc(f->weight, sqrt(2.0), 0); // Sets weight=1.
      gkyl_array_scale(f->weight, polarization_weight);
      gkyl_array_scale(f->weight, f->info.kperp2);
    }
    else {
      // allocate arrays for Poisson solver and Poisson solver 
      f->epsilon = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
      gkyl_array_set_offset(f->epsilon, polarization_weight, app->gk_geom->gxxj, 0*app->confBasis.num_basis);
      gkyl_array_set_offset(f->epsilon, polarization_weight, app->gk_geom->gxyj, 1*app->confBasis.num_basis);
      gkyl_array_set_offset(f->epsilon, polarization_weight, app->gk_geom->gyyj, 2*app->confBasis.num_basis);

      f->fem_poisson_perp = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->confBasis, 
        &f->info.poisson_bcs, f->epsilon, f->kSq, app->use_gpu);
    }
  }
  // Potential smoothing (in z) updater
  f->fem_parproj = gkyl_fem_parproj_new(&app->local, &app->local_ext, 
    &app->confBasis, f->info.fem_parbc, f->weight, app->use_gpu);

  f->phi_host = f->phi_smooth;  
  if (app->use_gpu) {
    f->phi_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    f->em_energy_red = gkyl_cu_malloc(sizeof(double[1]));
  }

  f->es_energy_fac = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  gkyl_array_shiftc(f->es_energy_fac, sqrt(2.0), 0); // Sets es_energy_fac=1.
  gkyl_array_scale(f->es_energy_fac, 0.5*polarization_weight);
  if (app->cdim == 1){
    gkyl_array_scale(f->es_energy_fac, 0.5*f->info.kperp2);
  }

  f->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  f->is_first_energy_write_call = true;
  if (app->cdim==1) {
    f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->confBasis, 
      1, GKYL_ARRAY_INTEGRATE_OP_SQ, app->use_gpu);
  }
  else {
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
      f->phi_wall_lo_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);

    f->phi_wall_lo_proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis, app->confBasis.poly_order+1,
      1, f->info.phi_wall_lo, &f->info.phi_wall_lo_ctx);

    // Compute phi_wall_lo at t = 0
    gkyl_proj_on_basis_advance(f->phi_wall_lo_proj, 0.0, &app->local_ext, f->phi_wall_lo_host);
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
      f->phi_wall_up_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);

    f->phi_wall_up_proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis, app->confBasis.poly_order+1,
      1, f->info.phi_wall_up, &f->info.phi_wall_up_ctx);

    // Compute phi_wall_up at t = 0
    gkyl_proj_on_basis_advance(f->phi_wall_up_proj, 0.0, &app->local_ext, f->phi_wall_up_host);
    if (app->use_gpu) // note: phi_wall_up_host is same as phi_wall_up when not on GPUs
      gkyl_array_copy(f->phi_wall_up, f->phi_wall_up_host);
  }

  return f;
}

void
gk_field_calc_phi_wall(gkyl_gyrokinetic_app *app, struct gk_field *field, double tm)
{
  if (field->has_phi_wall_lo && field->phi_wall_lo_evolve) {
    gkyl_proj_on_basis_advance(field->phi_wall_lo_proj, tm, &app->local_ext, field->phi_wall_lo_host);
    if (app->use_gpu) // note: phi_wall_lo_host is same as phi_wall_lo when not on GPUs
      gkyl_array_copy(field->phi_wall_lo, field->phi_wall_lo_host);
  }
  if (field->has_phi_wall_up && field->phi_wall_up_evolve) {
    gkyl_proj_on_basis_advance(field->phi_wall_up_proj, tm, &app->local_ext, field->phi_wall_up_host);
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
    // if adiabatic electrons, we only need ion density, not charge density
    if (field->gkfield_id == GKYL_GK_FIELD_ADIABATIC) 
      gkyl_array_accumulate_range(field->rho_c, 1.0, s->m0.marr, &app->local);
    else
      gkyl_array_accumulate_range(field->rho_c, s->info.charge, s->m0.marr, &app->local);
  } 
}

void
gk_field_calc_ambi_pot_sheath_vals(gkyl_gyrokinetic_app *app, struct gk_field *field, 
  const struct gkyl_array *fin[], struct gkyl_array *rhs[])
{
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];

    gk_species_bflux_rhs(app, s, &s->bflux, fin[i], rhs[i]);

    // Assumes symmetric sheath BCs for now only in 1D
    gkyl_ambi_bolt_potential_sheath_calc(field->ambi_pot, GKYL_LOWER_EDGE, 
      &app->lower_skin[0], &app->lower_ghost[0], app->gk_geom->jacobgeo, 
      s->bflux.gammai[0].marr, field->rho_c, field->sheath_vals[0]);
    gkyl_ambi_bolt_potential_sheath_calc(field->ambi_pot, GKYL_UPPER_EDGE, 
      &app->upper_skin[0], &app->upper_ghost[0], app->gk_geom->jacobgeo, 
      s->bflux.gammai[1].marr, field->rho_c, field->sheath_vals[1]);
  } 
}

// Compute the electrostatic potential
void
gk_field_rhs(gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  struct timespec wst = gkyl_wall_clock();
  if (field->gkfield_id == GKYL_GK_FIELD_ADIABATIC) { 
    // This is not currently right. There's some subtlety in how the sheath_vals are stored
    gkyl_ambi_bolt_potential_phi_calc(field->ambi_pot,
      &app->local, &app->local_ext,
      app->gk_geom->jacobgeo_inv, field->rho_c,
      field->sheath_vals[0], field->phi_fem);
    gkyl_fem_parproj_set_rhs(field->fem_parproj, field->phi_fem, field->phi_fem);
    gkyl_fem_parproj_solve(field->fem_parproj, field->phi_smooth);
  }
  else {
    if (app->cdim==1) {
      gkyl_fem_parproj_set_rhs(field->fem_parproj, field->rho_c, 0);
      gkyl_fem_parproj_solve(field->fem_parproj, field->phi_smooth);
    }
    else {
      gkyl_fem_parproj_set_rhs(field->fem_parproj, field->rho_c, 0);
      gkyl_fem_parproj_solve(field->fem_parproj, field->rho_c_smooth);

      gkyl_fem_poisson_perp_set_rhs(field->fem_poisson_perp, field->rho_c_smooth);
      gkyl_fem_poisson_perp_solve(field->fem_poisson_perp, field->phi_fem);

      gkyl_fem_parproj_set_rhs(field->fem_parproj, field->phi_fem, field->phi_fem);
      gkyl_fem_parproj_solve(field->fem_parproj, field->phi_smooth);
    }
  }
  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *field)
{
  double energy[1] = { 0.0 };
  if (app->use_gpu) {
    gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth, 
      app->grid.cellVolume, field->es_energy_fac, &app->local, field->em_energy_red);
    gkyl_cu_memcpy(energy, field->em_energy_red, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth, 
      app->grid.cellVolume, field->es_energy_fac, &app->local, energy);
  } 

  double energy_global[6] = { 0.0 };
  gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 6, energy, energy_global);
  
  gkyl_dynvec_append(field->integ_energy, tm, energy_global);
}

// release resources for field
void
gk_field_release(const gkyl_gyrokinetic_app* app, struct gk_field *f)
{
  gkyl_array_release(f->rho_c);
  gkyl_array_release(f->rho_c_smooth);
  gkyl_array_release(f->phi_fem);
  gkyl_array_release(f->phi_smooth);

  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    gkyl_array_release(f->apar_fem);
    gkyl_array_release(f->apardot_fem);
  }

  if (f->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
    gkyl_ambi_bolt_potential_release(f->ambi_pot);
    for (int i=0; i<2*app->cdim; ++i) 
      gkyl_array_release(f->sheath_vals[i]);
  } 
  else {
    if (app->cdim == 1) {
      gkyl_array_release(f->weight);
    }
    else {
      gkyl_array_release(f->epsilon);
      gkyl_fem_poisson_perp_release(f->fem_poisson_perp);
    }
  }
  gkyl_fem_parproj_release(f->fem_parproj);
  
  gkyl_dynvec_release(f->integ_energy);
  gkyl_array_integrate_release(f->calc_em_energy);
  if (app->use_gpu) {
    gkyl_array_release(f->phi_host);
    gkyl_cu_free(f->em_energy_red);
  }

  gkyl_array_release(f->phi_wall_lo);
  gkyl_array_release(f->phi_wall_up);
  if (f->has_phi_wall_lo && app->use_gpu) 
    gkyl_array_release(f->phi_wall_lo);
  if (f->has_phi_wall_up && app->use_gpu) 
    gkyl_array_release(f->phi_wall_up);

  gkyl_free(f);
}

