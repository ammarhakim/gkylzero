#include <gkyl_vlasov_priv.h>

// internal function that takes a single time-step using a single-step
// Strang-split scheme
struct gkyl_update_status
vlasov_update_strang_split(gkyl_vlasov_app* app, double dt0)
{
  int ns = app->num_species;  
  int nfs = app->num_fluid_species;  

  const struct gkyl_array *fin[ns];
  struct gkyl_array *fout[ns];
  const struct gkyl_array *fluidin[nfs];
  struct gkyl_array *fluidout[nfs];
  struct gkyl_update_status st = { .success = true };
  
  // time-stepper states
  enum {
    UPDATE_DONE = 0,
    PRE_UPDATE,
    FIRST_COUPLING_UPDATE,
    RK_STAGE_1, 
    RK_STAGE_2, 
    RK_STAGE_3,
    SECOND_COUPLING_UPDATE,
    UPDATE_REDO,
  } state = PRE_UPDATE;

  double tcurr = app->tcurr, dt = dt0;
  while (state != UPDATE_DONE) {
    switch (state) {
      case PRE_UPDATE:
        do {
          state = FIRST_COUPLING_UPDATE; // next state
            
          // copy old solution in case we need to redo this step
          // Note: only need to copy fluid and EM solution since
          // kinetic species is done completely explicitly 
          // (old solution is still in app->species->f unless update successful).
          for (int i=0; i<nfs; ++i) {
            gkyl_array_copy(app->fluid_species[i].fluid_dup, app->fluid_species[i].fluid);
          }
          if (app->has_field) {
            gkyl_array_copy(app->field->em_dup, app->field->em);
          }
        } while(0);
        break;
          
      case FIRST_COUPLING_UPDATE:
        do {
          state = RK_STAGE_1; // next state

          struct timespec fl_em1_tm = gkyl_wall_clock();
          vm_fluid_em_coupling_update(app, app->fl_em, tcurr, dt/2.0);
          app->stat.fl_em_tm += gkyl_time_diff_now_sec(fl_em1_tm);
        } while(0);
        break;

      case RK_STAGE_1:
        do {
          struct timespec rk3_s1_tm = gkyl_wall_clock();

          for (int i=0; i<ns; ++i) {
            fin[i] = app->species[i].f;
            fout[i] = app->species[i].f1;
          }
          for (int i=0; i<nfs; ++i) {
            fluidin[i] = app->fluid_species[i].fluid;
            fluidout[i] = app->fluid_species[i].fluid1;
          }
          vlasov_forward_euler(app, tcurr, dt, fin, fluidin, app->has_field ? app->field->em : 0,
            fout, fluidout, app->has_field ? app->field->em1 : 0,
            &st
          );
          // Limit fluid and EM solutions if desired (done after update as post-hoc fix)
          for (int i=0; i<nfs; ++i) {
            vm_fluid_species_limiter(app, &app->fluid_species[i], fluidout[i]);
          }
          if (app->has_field) {
            vm_field_limiter(app, app->field, app->field->em1);
          }
          dt = st.dt_actual;
          state = RK_STAGE_2;

          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s1_tm);
        } while(0);
        break;

      case RK_STAGE_2:
        do {
          struct timespec rk3_s2_tm = gkyl_wall_clock();

          for (int i=0; i<ns; ++i) {
            fin[i] = app->species[i].f1;
            fout[i] = app->species[i].fnew;
          }
          for (int i=0; i<nfs; ++i) {
            fluidin[i] = app->fluid_species[i].fluid1;
            fluidout[i] = app->fluid_species[i].fluidnew;
          }
          vlasov_forward_euler(app, tcurr+dt, dt, fin, fluidin, app->has_field ? app->field->em1 : 0,
            fout, fluidout, app->has_field ? app->field->emnew : 0,
            &st
          );
          // Limit fluid and EM solutions if desired (done after update as post-hoc fix)
          for (int i=0; i<nfs; ++i) {
            vm_fluid_species_limiter(app, &app->fluid_species[i], fluidout[i]);
          }
          if (app->has_field) {
            vm_field_limiter(app, app->field, app->field->emnew);
          }
          if (st.dt_actual < dt) {
            // collect stats
            double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
            app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0],
              dt_rel_diff);
            app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1],
              dt_rel_diff);
            app->stat.nstage_2_fail += 1;

            dt = st.dt_actual;
            state = UPDATE_REDO; // Need to re-do update, fetch copy of old time step
          } 
          else {
            for (int i=0; i<ns; ++i)
              array_combine(app->species[i].f1,
                3.0/4.0, app->species[i].f, 1.0/4.0, app->species[i].fnew, &app->species[i].local_ext);
            for (int i=0; i<nfs; ++i)
              array_combine(app->fluid_species[i].fluid1,
                3.0/4.0, app->fluid_species[i].fluid, 1.0/4.0, app->fluid_species[i].fluidnew, &app->local_ext);
            if (app->has_field)
              array_combine(app->field->em1,
                3.0/4.0, app->field->em, 1.0/4.0, app->field->emnew, &app->local_ext);

            state = RK_STAGE_3;
          }

          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s2_tm);
        } while(0);
        break;

      case RK_STAGE_3:
        do {
          struct timespec rk3_s3_tm = gkyl_wall_clock();

          for (int i=0; i<ns; ++i) {
            fin[i] = app->species[i].f1;
            fout[i] = app->species[i].fnew;
          }
          for (int i=0; i<nfs; ++i) {
            fluidin[i] = app->fluid_species[i].fluid1;
            fluidout[i] = app->fluid_species[i].fluidnew;
          }
          vlasov_forward_euler(app, tcurr+dt/2, dt, fin, fluidin, app->has_field ? app->field->em1 : 0,
            fout, fluidout, app->has_field ? app->field->emnew : 0,
            &st
          );
          // Limit fluid and EM solutions if desired (done after update as post-hoc fix)
          for (int i=0; i<nfs; ++i) {
            vm_fluid_species_limiter(app, &app->fluid_species[i], fluidout[i]);
          }
          if (app->has_field) {
            vm_field_limiter(app, app->field, app->field->emnew);
          }
          if (st.dt_actual < dt) {
            // collect stats
            double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
            app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0],
              dt_rel_diff);
            app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1],
              dt_rel_diff);
            app->stat.nstage_3_fail += 1;

            dt = st.dt_actual;
            state = UPDATE_REDO; // Need to re-do update, fetch copy of old time step

            app->stat.nstage_2_fail += 1;
          }
          else {
            for (int i=0; i<ns; ++i) {
              array_combine(app->species[i].f1,
                1.0/3.0, app->species[i].f, 2.0/3.0, app->species[i].fnew, &app->species[i].local_ext);
              gkyl_array_copy_range(app->species[i].f, app->species[i].f1, &app->species[i].local_ext);
            }
            for (int i=0; i<nfs; ++i) {
              array_combine(app->fluid_species[i].fluid1,
                1.0/3.0, app->fluid_species[i].fluid, 2.0/3.0, app->fluid_species[i].fluidnew, &app->local_ext);
              gkyl_array_copy_range(app->fluid_species[i].fluid, app->fluid_species[i].fluid1, &app->local_ext);
            }
            if (app->has_field) {
              array_combine(app->field->em1,
                1.0/3.0, app->field->em, 2.0/3.0, app->field->emnew, &app->local_ext);
              gkyl_array_copy_range(app->field->em, app->field->em1, &app->local_ext);
            }

            state = SECOND_COUPLING_UPDATE;
          }

          app->stat.rk3_tm += gkyl_time_diff_now_sec(rk3_s3_tm);
        } while(0);
        break;

      case SECOND_COUPLING_UPDATE:
        do {
          state = UPDATE_DONE; // next state

          struct timespec fl_em2_tm = gkyl_wall_clock();
          vm_fluid_em_coupling_update(app, app->fl_em, tcurr, dt/2.0);
          app->stat.fl_em_tm += gkyl_time_diff_now_sec(fl_em2_tm);
        } while(0);
        break;

      case UPDATE_REDO:
        do {
          state = PRE_UPDATE; // start all-over again
            
          // restore solution and retake step
          // Note: only need to restore fluid and EM solution since
          // kinetic species is done completely explicitly 
          // (old solution is still in app->species->f unless update successful).
          for (int i=0; i<nfs; ++i) {
            gkyl_array_copy(app->fluid_species[i].fluid, app->fluid_species[i].fluid_dup);
          }
          if (app->has_field) {
            gkyl_array_copy(app->field->em, app->field->em_dup);
          }
        } while(0);          
        break;

      case UPDATE_DONE: // unreachable code! (suppresses warning)
        break;
    }
  }

  return st;
}