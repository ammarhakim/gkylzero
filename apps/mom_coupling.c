#include <gkyl_moment_priv.h>

// initialize source solver: this should be called after all species
// and fields are initialized
void
moment_coupling_init(const struct gkyl_moment_app *app, struct moment_coupling *src)
{
  struct gkyl_moment_em_coupling_inp src_inp = {
    .grid = &app->grid,
    .nfluids = app->num_species,
    // if there is a field, need to update electric field too, otherwise just updating fluid
    .epsilon0 = app->field.epsilon0 ? app->field.epsilon0 : 0.0, 
  };

  for (int i=0; i<app->num_species; ++i)
    src_inp.param[i] = (struct gkyl_moment_em_coupling_data) {
      .type = app->species[i].eqn_type,
      .charge = app->species[i].charge,
      .mass = app->species[i].mass,
      .k0 = app->species[i].k0
    };

  src_inp.has_collision = app->has_collision;
  src_inp.gas_gamma = app->gas_gamma;
  for (int s=0; s<app->num_species; ++s)
    for (int r=0; r<app->num_species; ++r)
      src_inp.nu_base[s][r] = app->nu_base[s][r];

  src_inp.has_user_sources = app->has_user_sources;
  src_inp.k_Boltzmann = app->k_Boltzmann;

  // create updater to solve for sources
  src->slvr = gkyl_moment_em_coupling_new(src_inp);
}

// update sources: 'nstrang' is 0 for the first Strang step and 1 for
// the second step
void
moment_coupling_update(gkyl_moment_app *app, struct moment_coupling *src,
  int nstrang, double tcurr, double dt)
{
  int sidx[] = { 0, app->ndim };
  struct gkyl_array *fluids[GKYL_MAX_SPECIES];
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES];
  const struct gkyl_array *user_sources[GKYL_MAX_SPECIES];
  
  for (int i=0; i<app->num_species; ++i) {
    fluids[i] = app->species[i].f[sidx[nstrang]];
    
    if (app->species[i].proj_app_accel)
      gkyl_fv_proj_advance(app->species[i].proj_app_accel, tcurr, &app->local, app->species[i].app_accel);
    app_accels[i] = app->species[i].app_accel;

    if (app->species[i].proj_user_source
        && !(app->user_sources_set_only_once && app->user_sources_is_set))
    {
      gkyl_fv_proj_advance(app->species[i].proj_user_source, tcurr,
                           &app->local, app->species[i].user_source);
    }
    user_sources[i] = app->species[i].user_source;
  }
  app->user_sources_is_set = true;

  if (app->field.proj_app_current)
    gkyl_fv_proj_advance(app->field.proj_app_current, tcurr, &app->local, app->field.app_current);
  
  if (app->field.proj_ext_em) {

    if (!app->field.was_ext_em_computed)
      gkyl_fv_proj_advance(app->field.proj_ext_em, tcurr, &app->local, app->field.ext_em);

    if (app->field.is_ext_em_static)
      app->field.was_ext_em_computed = true;
    else
      app->field.was_ext_em_computed = false;
  }

  gkyl_moment_em_coupling_advance(src->slvr, dt, &app->local,
    fluids, app_accels,
    app->field.f[sidx[nstrang]], app->field.app_current, app->field.ext_em,
    user_sources);

  for (int i=0; i<app->num_species; ++i)
    moment_species_apply_bc(app, tcurr, &app->species[i], fluids[i]);

  moment_field_apply_bc(app, tcurr, &app->field, app->field.f[sidx[nstrang]]);
}

// free sources
void
moment_coupling_release(const struct moment_coupling *src)
{
  gkyl_moment_em_coupling_release(src->slvr);
}

