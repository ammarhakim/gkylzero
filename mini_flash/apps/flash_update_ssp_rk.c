#include <gkyl_flash_priv.h>
#include <gkyl_flash.h>

// Take a single time step with the 3-stage SSP-RK method with suggested
// timestep dt. This function will never take a timestep larger than dt,
// but may take one smaller. 
// The actual dt is returned in the timestep.
struct gkyl_update_status
gkyl_flash_update_ssp_rk3(gkyl_flash_app *app, double dt0)
{
  int ns = app->num_species;

  const struct gkyl_array *fin[ns];
  struct gkyl_array *fout[ns];
  const struct gkyl_array *fv_fluidin[ns];
  struct gkyl_array *fv_fluidout[ns];
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
            fv_fluidin[i] = app->fluid_species[i].fluid;
            fv_fluidout[i] = app->fluid_species[i].fluid1;
          }
          flash_forward_euler(app, tcurr, dt, fin, fv_fluidin, fout, fv_fluidout, &st);
          dt = st.dt_actual;

          state = RK_STAGE_2;

        } while(0);
        break;

      case RK_STAGE_2:
        do {
          for (int i=0; i<ns; ++i) {
            fin[i] = app->species[i].f1;
            fout[i] = app->species[i].fnew;
            fv_fluidin[i] = app->fluid_species[i].fluid1;
            fv_fluidout[i] = app->fluid_species[i].fluidnew;
          }
          flash_forward_euler(app, tcurr+dt, dt, fin, fv_fluidin, fout, fv_fluidout, &st);

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
              array_combine(app->fluid_species[i].fluid1, 
                3.0/4.0, app->fluid_species[i].fluid, 
                1.0/4.0, app->fluid_species[i].fluidnew, &app->fluid_species[i].local_ext);
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
            fv_fluidin[i] = app->fluid_species[i].fluid1;
            fv_fluidout[i] = app->fluid_species[i].fluidnew;
          }
          flash_forward_euler(app, tcurr+dt/2, dt, fin, fv_fluidin, fout, fv_fluidout, &st);
          
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
              array_combine(app->fluid_species[i].fluid1, 
                1.0/3.0, app->fluid_species[i].fluid, 
                2.0/3.0, app->fluid_species[i].fluidnew, &app->fluid_species[i].local_ext);
              gkyl_array_copy_range(app->fluid_species[i].fluid, 
                app->fluid_species[i].fluid1, &app->fluid_species[i].local_ext);
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
