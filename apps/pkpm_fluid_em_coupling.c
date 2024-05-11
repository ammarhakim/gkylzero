#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_pkpm_priv.h>

// initialize fluid-EM coupling object for PKPM system
struct pkpm_fluid_em_coupling* 
pkpm_fluid_em_coupling_init(struct gkyl_pkpm_app *app)
{
  struct pkpm_fluid_em_coupling *pkpm_em = gkyl_malloc(sizeof(struct pkpm_fluid_em_coupling));

  int num_species = app->num_species;
  double qbym[GKYL_MAX_SPECIES] = {0.0};
  for (int i=0; i<num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    qbym[i] = s->info.charge/s->info.mass;
  }
  // Initialize solver
  pkpm_em->slvr = gkyl_dg_calc_pkpm_em_coupling_new(&app->confBasis, &app->local, 
    num_species, qbym, app->field->info.epsilon0, app->use_gpu);

  return pkpm_em;  
}

void
pkpm_fluid_em_coupling_update(struct gkyl_pkpm_app *app, struct pkpm_fluid_em_coupling *pkpm_em, 
  double tcurr, double dt)
{
  int num_species = app->num_species;

  struct gkyl_array *fluids[GKYL_MAX_SPECIES];
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES];
  const struct gkyl_array *vlasov_pkpm_moms[GKYL_MAX_SPECIES];
  const struct gkyl_array *pkpm_u[GKYL_MAX_SPECIES];

  for (int i=0; i<num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    fluids[i] = s->fluid;
    app_accels[i] = s->app_accel; 
    // Compute applied accelerations if present and time-dependent.
    // Note: applied accelerations use proj_on_basis 
    // so does copy to GPU every call if app->use_gpu = true.
    if (s->app_accel_evolve) {
      pkpm_species_calc_app_accel(app, s, tcurr);
    }
    // Compute the PKPM moments from the kinetic equation at the current time 
    pkpm_species_moment_calc(&s->pkpm_moms, s->local, app->local, s->f);
    vlasov_pkpm_moms[i] = s->pkpm_moms.marr;
    // Compute the flow velocity at the current time
    gkyl_dg_calc_pkpm_vars_u(s->calc_pkpm_vars,
      vlasov_pkpm_moms[i], fluids[i], 
      s->cell_avg_prim, s->pkpm_u);
    pkpm_u[i] = s->pkpm_u;
  }
  // Compute external EM field or applied currents if present and time-dependent.
  // Note: external EM field and applied currents use proj_on_basis 
  // so does copy to GPU every call if app->use_gpu = true.
  if (app->has_field) {
    if (app->field->ext_em_evolve) {
      pkpm_field_calc_ext_em(app, app->field, tcurr);
    }
    if (app->field->app_current_evolve) {
      pkpm_field_calc_app_current(app, app->field, tcurr); 
    }
  }

  gkyl_dg_calc_pkpm_em_coupling_advance(pkpm_em->slvr, dt, 
    app_accels, app->field->ext_em, app->field->app_current, 
    vlasov_pkpm_moms, pkpm_u, fluids, app->field->em);

  for (int i=0; i<num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    pkpm_fluid_species_apply_bc(app, s, fluids[i]);
  }
  pkpm_field_apply_bc(app, app->field, app->field->em);
}

void
pkpm_fluid_em_coupling_release(struct gkyl_pkpm_app *app, struct pkpm_fluid_em_coupling *pkpm_em)
{
  gkyl_dg_calc_pkpm_em_coupling_release(pkpm_em->slvr);
  gkyl_free(pkpm_em);
}
