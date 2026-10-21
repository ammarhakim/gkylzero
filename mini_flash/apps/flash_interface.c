// Functions required for coupling with FLASH
//
#include <gkyl_flash_multib_priv.h>
#include <gkyl_flash_priv.h>

/**
 * Pass array of fluid velocity [ux, uy, uz] to set the velocity of the moving
 * reference frame for the kinetic solver.
 * The input array is finite volume data living on a grid with 2x the configuration space
 * resolution of the kinetic grid in the app.
 * 
 * @param app 
 * @param fluid_mom_in 
 */
void
gkyl_flash_app_set_fluid_vel(gkyl_flash_app *app, struct gkyl_array *fluid_vel_in)
{
  // Each species uses the same bulk velocity
  for (int i=0; i<app->num_species; ++i) {
    struct flash_species *s = &app->species[i];

    gkyl_flash_interp_vars_vel_to_dg(app->interp_vars, 
      &app->local_kinetic_ext, &app->local_fluid_ext, fluid_vel_in, s->fluid_rhou);

    // Ensure pkpm_moms are computed before scaling velocity by density
    flash_species_calc_pkpm_moms(app, &app->species[i], app->species[i].f);

    // Multiply by density because the PKPM updater looks for momentum density instead of velocity
    for (int d=0; d<3; ++d)
      gkyl_dg_mul_op(&app->confBasis, d, s->fluid_rhou, d, s->fluid_rhou, 0, s->pkpm_moms);
  }
}

void
gkyl_flash_multib_app_set_fluid_vel(gkyl_flash_multib_app *app, int blockID, struct gkyl_array *fluid_vel_in)
{
  gkyl_flash_app_set_fluid_vel(app->block_apps[blockID], fluid_vel_in);
}

/**
 * Pass array of moments [n, ux, uy, uz, T/m] to the app to set initial conditions.
 * The input array is finite volume data living on a grid with 2x the configuration space
 * resolution of the kinetic grid in the app.
 * 
 * 
 * @param app 
 * @param init_dist_moms 
 */
void
gkyl_flash_app_set_init_dist_moms(gkyl_flash_app *app, struct gkyl_array *init_dist_moms)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_flash_interp_vars_init_dist_moms_to_dg(app->interp_vars,
      &app->local_kinetic_ext, &app->local_fluid_ext, init_dist_moms, app->species[i].init_dist_moms);
  }
}

/**
 * Initialize the distribution function as a Maxwellian matching the moment array
 * of finite volume data [n, ux, uy, uz, T/m]. 
 * 
 * @param app 
 * @param init_dist_moms 
 */
void
gkyl_flash_app_init_maxwellian_from_moms(gkyl_flash_app *app, struct gkyl_array *init_dist_moms)
{
  // Pass the moment array to the app internals
  gkyl_flash_app_set_init_dist_moms(app, init_dist_moms);

  // Initialize the distribution function
  for (int i=0; i<app->num_species; ++i) {
    flash_species_calc_distf_from_moms(app, &app->species[i], app->species[i].init_dist_moms, 0.0);
    
    flash_species_apply_bc(app, &app->species[i], app->species[i].f);
  }
}

/**
 * Initialize the distribution function of a single-block app
 * as a Maxwellian matching the moment array of finite volume data [n, ux, uy, uz, T/m].
 * 
 * @param app 
 * @param blockID
 * @param init_dist_moms 
 */
void
gkyl_flash_multib_app_init_maxwellian_from_moms_block(gkyl_flash_multib_app *app, int blockID, struct gkyl_array *init_dist_moms)
{
  struct gkyl_flash_app *sbapp = app->block_apps[blockID];

  // Pass the moment array to the app internal
  gkyl_flash_app_set_init_dist_moms(sbapp, init_dist_moms);

  // Initialize the distribution function
  for (int i=0; i<sbapp->num_species; ++i) {
    flash_species_calc_distf_from_moms(sbapp, &sbapp->species[i], 
      sbapp->species[i].init_dist_moms, 0.0);
  
    flash_species_apply_bc(sbapp, &sbapp->species[i], sbapp->species[i].f);
  }

  // Apply boundary conditions
  struct gkyl_array *distf[app->num_species * app->num_local_blocks];
  struct gkyl_array *fluid[app->num_species * app->num_local_blocks];
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_flash_app *sbapp = app->block_apps[b];
    int lin_idx = b*app->num_species;
    for (int i=0; i<app->num_species; ++i) {
      distf[lin_idx + i] = sbapp->species[i].f;
      fluid[lin_idx + i] = sbapp->fluid_species[i].fluid;
    }
  }

  gkyl_flash_multib_app_apply_bc(app, app->tcurr, distf, fluid);
}

void
gkyl_flash_app_get_kinetic_fluxes(gkyl_flash_app *app, struct gkyl_array *coupling_fluxes)
{
  for (int i=0; i<app->num_species; ++i) {
    // Verify that array exists
    if (coupling_fluxes->elemsz == 0) {
      char msg[128];
      snprintf(msg, 128, 
        "Output coupling flux array is empty. Check that the proper number of arrays have been initialized. Exiting...");
      gkyl_exit(msg);
    }

    // Compute PKPM moments, [rho, p_par, p_perp].
    flash_species_calc_pkpm_moms(app, &app->species[i], app->species[i].f);

    // Compute coupling moments, self-consistent pressure
    // force (div(p_par b_hat)), and primitive variables
    flash_species_calc_pkpm_vars(app, &app->species[i], app->species[i].f);

    // Compute remaining PKPM variables: acceleration variables in kinetic equation and 
    // source distribution functions for Laguerre couplings.
    flash_species_calc_pkpm_update_vars(app, &app->species[i], app->species[i].f);

    // Compute fluxes
    flash_species_calc_coupling_fluxes(app, &app->species[i], app->species[i].f);

    // Interpolate to FV
    gkyl_flash_interp_vars_coupling_fluxes_to_fv(app->interp_vars, 
      &app->local_kinetic_ext, &app->local_fluid_ext, 
      app->species[i].flash_coupling_fluxes, coupling_fluxes);
  }
}

void
gkyl_flash_multib_app_get_kinetic_fluxes(gkyl_flash_multib_app *app, int blockID, struct gkyl_array *coupling_fluxes)
{
  gkyl_flash_app_get_kinetic_fluxes(app->block_apps[blockID], coupling_fluxes);
}

void 
gkyl_flash_app_set_vel_space_extents(gkyl_flash_app *app, double lower_vel[GKYL_MAX_VDIM], double upper_vel[GKYL_MAX_VDIM])
{
  for (int i=0; i<app->num_species; ++i) {
    flash_species_set_vel_space_extents(app, &app->species[i], lower_vel, upper_vel);
  }
}

void
gkyl_flash_multib_app_set_vel_space_extents(gkyl_flash_multib_app *app, int blockID, double lower_vel[GKYL_MAX_VDIM], double upper_vel[GKYL_MAX_VDIM])
{
  gkyl_flash_app_set_vel_space_extents(app->block_apps[blockID], lower_vel, upper_vel);
}


struct gkyl_range
gkyl_flash_app_get_fluid_conf_ext_range(gkyl_flash_app *app) 
{
  return app->local_fluid_ext;
}

struct gkyl_range
gkyl_flash_app_get_fluid_conf_range(gkyl_flash_app *app) 
{
  return app->local_fluid;
}

struct gkyl_range
gkyl_flash_multib_app_get_fluid_conf_ext_range(gkyl_flash_multib_app *app, int blockID) 
{
  return app->block_apps[blockID]->local_fluid_ext;
}

struct gkyl_range
gkyl_flash_multib_app_get_fluid_conf_range(gkyl_flash_multib_app *app, int blockID) 
{
  return app->block_apps[blockID]->local_fluid;
}

struct gkyl_update_status
gkyl_flash_kinetic_update_ssp_rk3(gkyl_flash_app *app, double dt0)
{
  int ns = app->num_species;

  const struct gkyl_array *fin[ns];
  struct gkyl_array *fout[ns];
  struct gkyl_update_status st = { .success = true };

  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  double dtmin = DBL_MAX;

  while (state != RK_COMPLETE ) {
    switch (state) {
      case RK_STAGE_1:
        do {
          for (int i=0; i<app->num_species; ++i) {
            fin[i] = app->species[i].f;
            fout[i] = app->species[i].f1;
          }
          flash_vlasov_forward_euler(app, tcurr, dt, fin, fout, &st);
          dt = st.dt_actual;

          state = RK_STAGE_2;

        } while(0);
        break;

      case RK_STAGE_2:
        do {
          for (int i=0; i<ns; ++i) {
            fin[i] = app->species[i].f1;
            fout[i] = app->species[i].fnew;
          }
          flash_vlasov_forward_euler(app, tcurr+dt, dt, fin, fout, &st);

          // Verify that the timestep is consistent with previous stages.
          // If not, take the smaller timestep and retry.
          if (st.dt_actual < dt) {
            double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;

            app->stat.nfail += 1;

            dt = st.dt_actual;
            state = RK_STAGE_1;
          }
          else {
            // No issues with the timestep, continue onto next stage.
            for (int i=0; i<ns; ++i) {

              array_combine(app->species[i].f1, 
                3.0/4.0, app->species[i].f, 
                1.0/4.0, app->species[i].fnew, &app->species[i].local_ext);
            }

            state = RK_STAGE_3;
          }
        } while(0);
        break;

      case RK_STAGE_3:
        do {
          for (int i=0; i<ns; ++i) {
            fin[i] = app->species[i].f1;
            fout[i] = app->species[i].fnew;
          }
          flash_vlasov_forward_euler(app, tcurr+dt/2, dt, fin, fout, &st);
          
          // Verify that the timestep is consistent with previous stages.
          // If not, take the smaller timestep and retry.
          if (st.dt_actual < dt) {
            double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;

            app->stat.nfail += 1;

            dt = st.dt_actual;
            state = RK_STAGE_1;
          }
          else {
            // Complete update to next timestep.
            for (int i=0; i<ns; ++i) {
              array_combine(app->species[i].f1, 
                1.0/3.0, app->species[i].f, 
                2.0/3.0, app->species[i].fnew, &app->species[i].local_ext);
              gkyl_array_copy_range(app->species[i].f, 
                app->species[i].f1, &app->species[i].local_ext);
            }
            state = RK_COMPLETE;
          }
        } while(0);
        break;
      case RK_COMPLETE:
        break;
    }
  }
  return st;
};

struct gkyl_update_status
gkyl_flash_app_kinetic_update(gkyl_flash_app *app, double dt)
{
  app->stat.nup += 1;

  struct timespec wst = gkyl_wall_clock();
  struct gkyl_update_status status;

  int ns = app->num_species;
  const struct gkyl_array *fin[ns];
  struct gkyl_array *fout[ns];

  for (int i=0; i<ns; ++i) {
    fin[i] = app->species[i].f;
    fout[i] = app->species[i].f1;
    flash_vlasov_forward_euler(app, app->tcurr, dt, fin, fout, &status);

    // gkyl_array_copy_range(app->species[i].f, fout[i], &app->species[i].local_ext);
    gkyl_array_copy_range(app->species[i].f, fout[i], &app->species[i].local);
  }

  app->tcurr += status.dt_actual;

  app->stat.total_tm += gkyl_time_diff_now_sec(wst);

  return status;
}

void
gkyl_flash_multib_app_kinetic_update(gkyl_flash_multib_app *app, double dt0)
{
  app->stat.nup += 1;

  struct timespec wst = gkyl_wall_clock();
  struct gkyl_update_status status;

  int ns = app->num_species;
  const struct gkyl_array *fin[ns];
  struct gkyl_array *fout[ns];

  // Compute maximum timestep
  // double dt = fmin(dt0, gkyl_flash_multib_app_compute_dt_max(app));
  double dt = dt0;

  // Call forward Euler update
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_flash_app *sbapp = app->block_apps[b];

    gkyl_flash_app_kinetic_update(sbapp, dt);
  }

  app->tcurr += status.dt_actual;

  app->stat.total_tm += gkyl_time_diff_now_sec(wst);

  // Apply boundary conditions
  struct gkyl_array *distf[app->num_species * app->num_local_blocks];
  struct gkyl_array *fluid[app->num_species * app->num_local_blocks];
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_flash_app *sbapp = app->block_apps[b];
    int lin_idx = b*app->num_species;
    for (int i=0; i<app->num_species; ++i) {
      distf[lin_idx + i] = sbapp->species[i].f;
      fluid[lin_idx + i] = sbapp->fluid_species[i].fluid;
    }
  }

  gkyl_flash_multib_app_apply_bc(app, app->tcurr, distf, fluid);

}