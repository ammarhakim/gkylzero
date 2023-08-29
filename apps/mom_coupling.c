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
    // linear ramping function for slowing turning on applied accelerations, E fields, or currents
    .t_ramp_E = app->field.t_ramp_E ? app->field.t_ramp_E : 0.0,
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

  // create updater to solve for sources
  src->slvr = gkyl_moment_em_coupling_new(src_inp);

  for (int n=0; n<app->num_species; ++n) {
    int meqn = app->species[n].num_equations;
    src->rhs[n] = mkarr(false, meqn, app->local_ext.volume);
    src->non_ideal_cflrate[n] = mkarr(false, 1, app->local_ext.volume); 
  }

  // create grid and ranges for non-ideal variables (grid is in computational space)
  // this grid is the grid of node values
  int ghost[3] = { 2, 2, 2 };
  double non_ideal_lower[3] = {0.0};
  double non_ideal_upper[3] = {0.0};
  int non_ideal_cells[3] = {0};
  // non-ideal terms (e.g., heat flux tensor) grid has one "extra" cell and is half a grid cell larger past the lower and upper domain
  for (int d=0; d<app->ndim; ++d) {
    non_ideal_lower[d] = app->grid.lower[d] - (app->grid.upper[d]-app->grid.lower[d])/(2.0* (double) app->grid.cells[d]);
    non_ideal_upper[d] = app->grid.upper[d] + (app->grid.upper[d]-app->grid.lower[d])/(2.0* (double) app->grid.cells[d]);
    non_ideal_cells[d] = app->grid.cells[d] + 1;
  }
  gkyl_rect_grid_init(&src->non_ideal_grid, app->ndim, non_ideal_lower, non_ideal_upper, non_ideal_cells);
  gkyl_create_grid_ranges(&app->grid, ghost, &src->non_ideal_local_ext, &src->non_ideal_local);
  // In Gradient-closure case, non-ideal variables are 10 heat flux tensor components
  for (int n=0;  n<app->num_species; ++n) 
    src->non_ideal_vars[n] = mkarr(false, 10, src->non_ideal_local_ext.volume);

  // check if gradient-closure is present
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure) {
      struct gkyl_ten_moment_grad_closure_inp grad_closure_inp = {
        .grid = &app->grid,
        .k0 = app->species[i].k0,
      };
      src->grad_closure_slvr[i] = gkyl_ten_moment_grad_closure_new(grad_closure_inp);      
    }
  }
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
  const struct gkyl_array *rhs_const[GKYL_MAX_SPECIES];
  
  for (int i=0; i<app->num_species; ++i) {
    fluids[i] = app->species[i].f[sidx[nstrang]];

    if (app->species[i].proj_app_accel)
      gkyl_fv_proj_advance(app->species[i].proj_app_accel, tcurr, &app->local, app->species[i].app_accel);
    app_accels[i] = app->species[i].app_accel;

    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure) {
      gkyl_ten_moment_grad_closure_advance(src->grad_closure_slvr[i], 
        &src->non_ideal_local, &app->local, 
        app->species[i].f[sidx[nstrang]], app->field.f[sidx[nstrang]], 
        src->non_ideal_cflrate[i], src->non_ideal_vars[i], src->rhs[i]);
    }
  }
  
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

  // Get the RHS pointer for accumulation during source update
  for (int i=0; i<app->num_species; ++i)
    rhs_const[i] = src->rhs[i];

  gkyl_moment_em_coupling_advance(src->slvr, tcurr, dt, &app->local,
    fluids, app_accels, rhs_const, 
    app->field.f[sidx[nstrang]], app->field.app_current, app->field.ext_em);

  for (int i=0; i<app->num_species; ++i)
    moment_species_apply_bc(app, tcurr, &app->species[i], fluids[i]);

  moment_field_apply_bc(app, tcurr, &app->field, app->field.f[sidx[nstrang]]);
}

// free sources
void
moment_coupling_release(const struct gkyl_moment_app *app, const struct moment_coupling *src)
{
  gkyl_moment_em_coupling_release(src->slvr);
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_release(src->rhs[i]);
    gkyl_array_release(src->non_ideal_cflrate[i]);
    gkyl_array_release(src->non_ideal_vars[i]);
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure)
      gkyl_ten_moment_grad_closure_release(src->grad_closure_slvr[i]);
  }
}

