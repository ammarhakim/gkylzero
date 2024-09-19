#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_poisson_priv.h>

#include <assert.h>
#include <float.h>
#include <time.h>

void
vp_field_calc_ext_em(gkyl_vlasov_poisson_app *app, struct vp_field *field, double tm)
{
  gkyl_eval_on_nodes_advance(field->ext_em_proj, tm, &app->local, field->ext_em_host);
  if (app->use_gpu) {
    gkyl_array_copy(field->ext_em, field->ext_em_host);
  }
}

struct vp_field*
vp_field_new(struct gkyl_vp *vp, struct gkyl_vlasov_poisson_app *app)
{
  // Initialize field object.

  struct vp_field *vpf = gkyl_malloc(sizeof(struct vp_field));

  vpf->info = vp->field;

  vpf->field_id = vpf->info.field_id ? vpf->info.field_id : GKYL_VP_FIELD_PHI;

  // Allocate arrays for charge density.
  vpf->rho_c        = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  vpf->rho_c_global = mkarr(app->use_gpu, app->confBasis.num_basis, app->global_ext.volume);

  // Allocate arrays for electrostatic potential.
  vpf->phi        = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  vpf->phi_global = mkarr(app->use_gpu, app->confBasis.num_basis, app->global_ext.volume);

  // Create global subrange we'll copy the field solver solution from (into local).
  int intersect = gkyl_sub_range_intersect(&vpf->global_sub_range, &app->global, &app->local);

  // Set the permittivity in the Poisson equation.
  vpf->epsilon = mkarr(app->use_gpu, app->confBasis.num_basis, app->global_ext.volume);
  gkyl_array_clear(vpf->epsilon, 0.0);
  gkyl_array_shiftc(vpf->epsilon, vpf->info.permittivity*pow(sqrt(2.0),app->cdim), 0);

  // Create Poisson solver.
  vpf->fem_poisson = gkyl_fem_poisson_new(&app->global, &app->grid, app->confBasis,
    &vpf->info.poisson_bcs, vpf->epsilon, NULL, true, app->use_gpu);

  // Initialize external potentials (always used by implicit fluid sources, so always initialize) 
  vpf->ext_em = mkarr(app->use_gpu, 4*app->confBasis.num_basis, app->local_ext.volume);
  vpf->has_ext_em = vpf->ext_em_evolve = false;
  if (vpf->info.ext_em) {
    vpf->has_ext_em = true;
    if (vpf->info.ext_em_evolve) {
      vpf->ext_em_evolve = vpf->info.ext_em_evolve;
    }

    vpf->ext_em_host = app->use_gpu? mkarr(false, vpf->ext_em->ncomp, vpf->ext_em->size)
                                   : gkyl_array_acquire(vpf->ext_em);

    vpf->ext_em_proj = gkyl_eval_on_nodes_new(&app->grid, &app->confBasis,
      4, vpf->info.ext_em, vpf->info.ext_em_ctx);

    // Compute the external potentials.
    vp_field_calc_ext_em(app, vpf, 0.0);
  }

  vpf->phi_host = vpf->phi;
  if (app->use_gpu) {
    vpf->phi_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    vpf->es_energy_red = gkyl_cu_malloc(sizeof(double[1]));
    vpf->es_energy_red_global = gkyl_cu_malloc(sizeof(double[1]));
  } else {
    vpf->es_energy_red = gkyl_malloc(sizeof(double[1]));
    vpf->es_energy_red_global = gkyl_malloc(sizeof(double[1]));
  }

  vpf->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 1);

  vpf->es_energy_fac = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_shiftc(vpf->es_energy_fac, pow(sqrt(2.0),app->cdim), 0); // Sets es_energy_fac=1.

  vpf->calc_es_energy = gkyl_array_integrate_new(&app->grid, &app->confBasis,
    1, GKYL_ARRAY_INTEGRATE_OP_GRAD_SQ, app->use_gpu);
  vpf->is_first_energy_write_call = true;

  return vpf;
}

void
vp_field_accumulate_rho_c(gkyl_vlasov_poisson_app *app, struct vp_field *field,
  const struct gkyl_array *fin[])
{
  // Calcualte the charge density.

  gkyl_array_clear(field->rho_c, 0.0);

  for (int i=0; i<app->num_species; ++i) {
    struct vp_species *vps = &app->species[i];

    vp_species_moment_calc(&vps->m0, vps->local, app->local, fin[i]);

    gkyl_array_accumulate_range(field->rho_c, vps->info.charge, vps->m0.marr, &app->local);
  }
}

void
vp_field_rhs(gkyl_vlasov_poisson_app *app, struct vp_field *field)
{
  // Compute the electrostatic potential.

  struct timespec wst = gkyl_wall_clock();
  // Gather charge density into global array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, field->rho_c, field->rho_c_global);

  // Solve the Poisson problem.
  gkyl_fem_poisson_set_rhs(field->fem_poisson, field->rho_c_global);
  gkyl_fem_poisson_solve(field->fem_poisson, field->phi_global);

  // Copy the portion of global potential corresponding to this MPI pcross to the local potential.
  gkyl_array_copy_range_to_range(field->phi, field->phi_global, &app->local, &field->global_sub_range);
  
  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);
}

void
vp_field_calc_energy(gkyl_vlasov_poisson_app *app, double tm, const struct vp_field *field)
{
  gkyl_array_integrate_advance(field->calc_es_energy, field->phi,
    app->grid.cellVolume, field->es_energy_fac, &app->local, field->es_energy_red);

  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->es_energy_red, field->es_energy_red_global);

  double energy_global[1] = { 0.0 };
  if (app->use_gpu)
    gkyl_cu_memcpy(energy_global, field->es_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  else
    energy_global[0] = field->es_energy_red_global[0];

  gkyl_dynvec_append(field->integ_energy, tm, energy_global);
}
 
void
vp_field_release(const gkyl_vlasov_poisson_app* app, struct vp_field *vpf)
{
  // Release resources for Vlasov-Poisson field.

  gkyl_dynvec_release(vpf->integ_energy);
  gkyl_array_integrate_release(vpf->calc_es_energy);
  gkyl_array_release(vpf->es_energy_fac);

  if (app->use_gpu) {
    gkyl_array_release(vpf->phi_host);
    gkyl_cu_free(vpf->es_energy_red);
    gkyl_cu_free(vpf->es_energy_red_global);
  } else {
    gkyl_free(vpf->es_energy_red);
    gkyl_free(vpf->es_energy_red_global);
  }

  gkyl_fem_poisson_release(vpf->fem_poisson);

  gkyl_array_release(vpf->epsilon);

  gkyl_array_release(vpf->ext_em);
  if (vpf->has_ext_em) {
    gkyl_array_release(vpf->ext_em_host);
    gkyl_eval_on_nodes_release(vpf->ext_em_proj);
  }
  gkyl_array_release(vpf->phi);
  gkyl_array_release(vpf->phi_global);

  gkyl_array_release(vpf->rho_c_global);
  gkyl_array_release(vpf->rho_c);

  free(vpf);
}
