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

  // allocate arrays for Poisson smoothing and solver
  f->weight = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_shiftc(f->weight, sqrt(2.0), 0); // Sets weight=1.

  f->epsilon = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
  // Linearized polarization density
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    gkyl_array_accumulate_range(f->epsilon, 
      s->info.mass/(f->info.bmag_fac*f->info.bmag_fac), s->polarization_density, &app->local_ext);
  }

  f->kSq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);

  f->fem_parproj = gkyl_fem_parproj_new(&app->local, &app->local_ext, 
    &app->confBasis, f->info.fem_parbc, f->weight, app->use_gpu);
  f->fem_poisson_perp = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->confBasis, 
    &f->info.poisson_bcs, f->epsilon, f->kSq, app->use_gpu);

  f->phi_host = f->phi_smooth;  
  if (app->use_gpu) {
    f->phi_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    f->em_energy_red = gkyl_cu_malloc(sizeof(double[1]));
  }

  f->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  f->is_first_energy_write_call = true;
  f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->confBasis, 
    1, GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ, app->use_gpu);

  // Initialize wall potential
  int Nbasis_surf = app->confBasis.num_basis/(app->confBasis.poly_order + 1); // *only valid for tensor bases for cdim > 1*
  f->phi_wall = mkarr(app->use_gpu, Nbasis_surf, app->local_ext.volume);

  return f;
}

void
gk_field_accumulate_rho_c(gkyl_gyrokinetic_app *app, struct gk_field *field, 
  const struct gkyl_array *fin[])
{
  gkyl_array_clear(field->rho_c, 0.0);
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];

    gk_species_moment_calc(&s->m0, s->local, app->local, fin[i]);
    gkyl_array_accumulate_range(field->rho_c, 1.0, s->m0.marr, &app->local);
  } 
}

// Compute the electrostatic potential
void
gk_field_rhs(gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_fem_parproj_set_rhs(field->fem_parproj, field->rho_c, 0);
  gkyl_fem_parproj_solve(field->fem_parproj, field->rho_c_smooth);

  gkyl_fem_poisson_perp_set_rhs(field->fem_poisson_perp, field->rho_c_smooth);
  gkyl_fem_poisson_perp_solve(field->fem_poisson_perp, field->phi_fem);

  gkyl_fem_parproj_set_rhs(field->fem_parproj, field->phi_fem, field->phi_fem);
  gkyl_fem_parproj_solve(field->fem_parproj, field->phi_smooth);

  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *field)
{
  double energy[1] = { 0.0 };
  if (app->use_gpu) {
    gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth, 
      app->grid.cellVolume, field->weight, &app->local, field->em_energy_red);
    gkyl_cu_memcpy(energy, field->em_energy_red, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth, 
      app->grid.cellVolume, field->weight, &app->local, energy);
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

  gkyl_array_release(f->weight);
  gkyl_array_release(f->epsilon);
  gkyl_array_release(f->kSq);
  gkyl_fem_parproj_release(f->fem_parproj);
  gkyl_fem_poisson_perp_release(f->fem_poisson_perp);
  
  gkyl_dynvec_release(f->integ_energy);
  gkyl_array_integrate_release(f->calc_em_energy);
  if (app->use_gpu) {
    gkyl_array_release(f->phi_host);
    gkyl_cu_free(f->em_energy_red);
  }

  gkyl_array_release(f->phi_wall);

  gkyl_free(f);
}

