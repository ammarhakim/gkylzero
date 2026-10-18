#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_pkpm_priv.h>

// initialize implicit fluid-EM coupling object for PKPM system
struct pkpm_fluid_em_coupling *
pkpm_fluid_em_coupling_init(struct gkyl_pkpm_app *app)
{
  struct pkpm_fluid_em_coupling *pkpm_em = gkyl_malloc(sizeof(struct pkpm_fluid_em_coupling));

  int num_species = app->num_species;
  double qbym[GKYL_MAX_SPECIES] = {0.0};
  for (int i = 0; i < num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    qbym[i] = s->info.charge / s->info.mass;
  }

  // Initialize solver
  pkpm_em->slvr = gkyl_dg_calc_pkpm_em_coupling_new(
    &app->confBasis, &app->local, num_species, qbym, app->field->info.epsilon0,
    app->field->info.is_static, app->use_gpu
  );

  return pkpm_em;
}

void
pkpm_fluid_em_coupling_update(
  struct gkyl_pkpm_app *app, struct pkpm_fluid_em_coupling *pkpm_em, double tcurr, double dt
)
{
  int num_species = app->num_species;

  struct gkyl_array *fluids[GKYL_MAX_SPECIES];
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES];

  // The nodal implicit solve advances the fluid state [rho, M, E] (momentum by the solve,
  // energy by the work of the field) and the electric field in place. It works with the
  // PHYSICAL state: the evolved fluid and field carry the position-map Jacobian, which is
  // divided out before and restored after (identity map: copies).
  for (int i = 0; i < num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    gkyl_vlasov_position_map_divide_jacobpos_conf(app->pos_map, &app->local, s->fluid->ncomp, s->fluid, s->fluid_no_J);
    fluids[i] = s->fluid_no_J;
  }
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local, app->field->em->ncomp, app->field->em, app->field->em_no_J
  );

  // Compute applied acceleration, external EM fields, and applied currents if present.
  // Note: these quantities use proj_on_basis so we copy to GPU every call if app->use_gpu = true.
  if (app->field->ext_em_evolve) {
    pkpm_field_calc_ext_em(app, app->field, tcurr);
  }
  if (app->field->app_current_evolve) {
    pkpm_field_calc_app_current(app, app->field, tcurr);
  }
  for (int i = 0; i < num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    if (s->app_accel_evolve) {
      pkpm_species_calc_app_accel(app, s, tcurr);
    }
    app_accels[i] = s->app_accel;
  }

  gkyl_dg_calc_pkpm_em_coupling_advance(
    pkpm_em->slvr, dt, app_accels, app->field->ext_em, app->field->app_current_no_J, fluids,
    app->field->em_no_J
  );

  for (int i = 0; i < num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    gkyl_vlasov_position_map_rescale_jacobpos_conf(app->pos_map, &app->local, s->fluid_no_J, s->fluid);
    pkpm_species_apply_positivity(app, s, s->fluid);
    pkpm_fluid_species_apply_bc(app, s, s->fluid);
  }
  gkyl_vlasov_position_map_rescale_jacobpos_conf(
    app->pos_map, &app->local, app->field->em_no_J, app->field->em
  );
  pkpm_field_apply_bc(app, app->field, app->field->em);
}

void
pkpm_fluid_em_coupling_release(struct gkyl_pkpm_app *app, struct pkpm_fluid_em_coupling *pkpm_em)
{
  gkyl_dg_calc_pkpm_em_coupling_release(pkpm_em->slvr);
  gkyl_free(pkpm_em);
}
