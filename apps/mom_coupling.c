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
    .mu0 = app->field.mu0 ? app->field.mu0 : 0.0, 
    // linear ramping function for slowing turning on applied accelerations, E fields, or currents
    .t_ramp_E = app->field.t_ramp_ext_em ? app->field.t_ramp_ext_em : 0.0,
    .t_ramp_curr = app->field.t_ramp_curr ? app->field.t_ramp_curr : 0.0,
  };

  for (int i=0; i<app->num_species; ++i)
    src_inp.param[i] = (struct gkyl_moment_em_coupling_data) {
      .type = app->species[i].eqn_type,
      .charge = app->species[i].charge,
      .mass = app->species[i].mass,
      // If gradient-based closure is present, k0=0.0 in source solve to avoid applying local closure
      .k0 = app->species[i].has_grad_closure ? 0.0 : app->species[i].k0,
    };

  src_inp.has_collision = app->has_collision;
  for (int s=0; s<app->num_species; ++s)
    for (int r=0; r<app->num_species; ++r)
      src_inp.nu_base[s][r] = app->nu_base[s][r];

  src_inp.has_nT_sources = false;
  for (int i=0; i<app->num_species; ++i)
    if (app->species[i].proj_nT_source)
      src_inp.has_nT_sources = true;

  // check for relativistic-species
  bool use_rel = 0;
  for (int n=0;  n<app->num_species; ++n) 
    use_rel = use_rel || (app->species[n].eqn_type == GKYL_EQN_COLDFLUID_SR) || app->field.use_explicit_em_coupling;

  // save the use-rel bool
  src_inp.use_rel = use_rel;

  // check for explicit em-coupling
  src_inp.use_explicit_em_coupling = 0;
  if (app->field.use_explicit_em_coupling)
    src_inp.use_explicit_em_coupling = 1;

  // create updater to solve for sources
  src->slvr = gkyl_moment_em_coupling_new(src_inp);

  for (int n=0; n<app->num_species; ++n) {
    int meqn = app->species[n].num_equations;
    src->pr_rhs[n] = mkarr(false, meqn, app->local_ext.volume);
    src->non_ideal_cflrate[n] = mkarr(false, 1, app->local_ext.volume);
  }

  int ghost[3] = { 1, 1, 1 };
  // create non-ideal local extended range from local range
  // has one additional cell in each direction because non-ideal variables are stored at cell vertices
  // gkyl_create_ranges(&app->local, ghost, &src->non_ideal_local_ext, &src->non_ideal_local);
  gkyl_create_vertex_ranges(&app->local, ghost, &src->non_ideal_local_ext, &src->non_ideal_local);

  // check if gradient-closure is present
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure) {
      struct gkyl_ten_moment_grad_closure_inp grad_closure_inp = {
        .grid = &app->grid,
        .k0 = app->species[i].k0,
	.cfl = app->cfl,
	.comm = app->comm,
	.mass = app->species[i].mass
      };
      src->grad_closure_slvr[i] = gkyl_ten_moment_grad_closure_new(grad_closure_inp);
    }
  }
}

// update sources: 'nstrang' is 0 for the first Strang step and 1 for
// the second step
struct gkyl_update_status
moment_coupling_update(gkyl_moment_app *app, struct moment_coupling *src,
  int nstrang, double tcurr, double dt)
{
  int sidx[] = { 0, app->ndim };
  struct gkyl_array *fluids[GKYL_MAX_SPECIES];
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES];
  const struct gkyl_array *pr_rhs_const[GKYL_MAX_SPECIES];
  const struct gkyl_array *nT_sources[GKYL_MAX_SPECIES];

  double dt_suggested = DBL_MAX;
  struct gkyl_ten_moment_grad_closure_status stat;

  for (int i=0; i<app->num_species; ++i) {
    fluids[i] = app->species[i].f[sidx[nstrang]];

    if (app->species[i].proj_app_accel) {
      
      if (!app->species[i].was_app_accel_computed)
        gkyl_fv_proj_advance(app->species[i].proj_app_accel, tcurr, &app->local, app->species[i].app_accel);
      
      if (app->species[i].is_app_accel_static)
        app->species[i].was_app_accel_computed = true;
      else
        app->species[i].was_app_accel_computed = false;
    }
    app_accels[i] = app->species[i].app_accel;

    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure) {
      // non-ideal variables defined on an extended range with one additional "cell" in each direction
      // this additional cell accounts for the fact that non-ideal variables are stored at cell vertices
      stat = gkyl_ten_moment_grad_closure_advance(src->grad_closure_slvr[i],
        &src->non_ideal_local, &app->local,
        app->species[i].f[sidx[nstrang]], app->field.f[sidx[nstrang]],
        src->non_ideal_cflrate[i], dt, src->pr_rhs[i]);
      
      if (!stat.success)
        return (struct gkyl_update_status) {
          .success = false,
          .dt_suggested = stat.dt_suggested
        };
    
      dt_suggested = fmin(dt_suggested, stat.dt_suggested);
    }
  }

  if (app->field.proj_app_current)
    gkyl_fv_proj_advance(app->field.proj_app_current, tcurr, &app->local, app->field.app_current);
  if ((app->field.proj_app_current) && (app->field.use_explicit_em_coupling)){
    // TEMP: 2x on dt
    gkyl_fv_proj_advance(app->field.proj_app_current, tcurr + dt*2, &app->local, app->field.app_current1);
    gkyl_fv_proj_advance(app->field.proj_app_current, tcurr + 2*dt/2.0, &app->local, app->field.app_current2);
  }

  if (app->field.proj_ext_em) {

    if (!app->field.was_ext_em_computed)
      gkyl_fv_proj_advance(app->field.proj_ext_em, tcurr, &app->local, app->field.ext_em);

    if (app->field.is_ext_em_static)
      app->field.was_ext_em_computed = true;
    else
      app->field.was_ext_em_computed = false;
  }

  // Get the RHS pointer for accumulation during source update
  for (int i=0; i<app->num_species; ++i) {
    pr_rhs_const[i] = src->pr_rhs[i];
  }

  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].proj_nT_source
        && !(app->species[i].nT_source_set_only_once
             && app->species[i].nT_source_is_set))
    {
      gkyl_fv_proj_advance(app->species[i].proj_nT_source, tcurr,
                           &app->local, app->species[i].nT_source);
    }
    nT_sources[i] = app->species[i].nT_source;
    app->species[i].nT_source_is_set = true;
  }

  if (app->field.use_explicit_em_coupling)
    gkyl_moment_em_coupling_explicit_advance(src->slvr, tcurr, dt, &app->local,
      fluids, app_accels, pr_rhs_const, 
      app->field.f[sidx[nstrang]], app->field.app_current, app->field.app_current1,
      app->field.app_current2, app->field.ext_em, 
      nT_sources, app->field.proj_app_current,nstrang);
  else
    gkyl_moment_em_coupling_implicit_advance(src->slvr, tcurr, dt, &app->local,
      fluids, app_accels, pr_rhs_const, 
      app->field.f[sidx[nstrang]], app->field.app_current, app->field.ext_em, 
      nT_sources);

  for (int i=0; i<app->num_species; ++i) {
    moment_species_apply_bc(app, tcurr, &app->species[i], fluids[i]);
  }
  if (app->has_field) {
    moment_field_apply_bc(app, tcurr, &app->field, app->field.f[sidx[nstrang]]);
  }
  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested
  };
}

// free sources
void
moment_coupling_release(const struct gkyl_moment_app *app, const struct moment_coupling *src)
{
  gkyl_moment_em_coupling_release(src->slvr);
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_release(src->pr_rhs[i]);
    gkyl_array_release(src->non_ideal_cflrate[i]);
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure)
      gkyl_ten_moment_grad_closure_release(src->grad_closure_slvr[i]);
  }
}

