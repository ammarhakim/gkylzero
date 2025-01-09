#include <gkyl_gyrokinetic_multib_priv.h>

static void
gyrokinetic_multib_forward_euler(struct gkyl_gyrokinetic_multib_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], struct gkyl_array *fout[], 
  const struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[], 
  struct gkyl_update_status *st)
{
  // Take a forward Euler step with the suggested time-step dt. This may
  // not be the actual time-step taken. However, the function will never
  // take a time-step larger than dt even if it is allowed by
  // stability. The actual time-step and dt_suggested are returned in
  // the status object.
  app->stat.nfeuler += 1;

  double dtmin = DBL_MAX;

  // Compute the time rate of change of the distributions, df/dt.
  for (int b=0; b<app->num_local_blocks; ++b) {
    int li_charged = b * app->num_species;
    int li_neut = b * app->num_neut_species;
    gyrokinetic_rhs(app->singleb_apps[b], tcurr, dt, &fin[li_charged], &fout[li_charged],
      &fin_neut[li_neut], &fout_neut[li_neut], st);
    dtmin = fmin(dtmin, st->dt_actual);
  }

  // Compute minimum time-step across all processors.
  double dtmin_local = dtmin, dtmin_global;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dtmin_local, &dtmin_global);
  st->dt_actual = dtmin_global;

  // Complete update of distribution functions.
  double dta = st->dt_actual;
  for (int b=0; b<app->num_local_blocks; ++b) {
    int li_charged = b * app->num_species;
    int li_neut = b * app->num_neut_species;
    for (int i=0; i<app->num_species; ++i) {
      gkyl_array_accumulate(gkyl_array_scale(fout[li_charged+i], dta), 1.0, fin[li_charged+i]);
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      if (!app->singleb_apps[b]->neut_species[i].info.is_static) {
        gkyl_array_accumulate(gkyl_array_scale(fout_neut[li_neut+i], dta), 1.0, fin_neut[li_neut+i]);
      }
    }
  }

}

struct gkyl_update_status
gyrokinetic_multib_update_ssp_rk3(struct gkyl_gyrokinetic_multib_app* app, double dt0)
{
  // Take time-step using the RK3 method. Also sets the status object
  // which has the actual and suggested dts used. These can be different
  // from the actual time-step.
  const struct gkyl_array *fin[app->num_species * app->num_local_blocks];
  struct gkyl_array *fout[app->num_species * app->num_local_blocks];
  const struct gkyl_array *fin_neut[app->num_neut_species * app->num_local_blocks];
  struct gkyl_array *fout_neut[app->num_neut_species * app->num_local_blocks];
  struct gkyl_update_status st = { .success = true };

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        for (int b=0; b<app->num_local_blocks; ++b) {
          struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
          int li_charged = b * app->num_species;
          int li_neut = b * app->num_neut_species;
          for (int i=0; i<app->num_species; ++i) {
            fin[li_charged+i] = sbapp->species[i].f;
            fout[li_charged+i] = sbapp->species[i].f1;
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            fin_neut[li_neut+i] = sbapp->neut_species[i].f;
	    fout_neut[li_neut+i] = sbapp->neut_species[i].f1;
          }
        }

        gyrokinetic_multib_forward_euler(app, tcurr, dt, fin, fout, fin_neut, fout_neut, &st);
        // Compute the fields and apply BCs.
        gyrokinetic_multib_calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

        dt = st.dt_actual;
        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int b=0; b<app->num_local_blocks; ++b) {
          struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
          int li_charged = b * app->num_species;
          int li_neut = b * app->num_neut_species;
          for (int i=0; i<app->num_species; ++i) {
            fin[li_charged+i] = sbapp->species[i].f1;
            fout[li_charged+i] = sbapp->species[i].fnew;
          }
          for (int i=0; i<app->num_neut_species; ++i) {
	    fin_neut[li_neut+i] = sbapp->neut_species[i].f1;
	    fout_neut[li_neut+i] = sbapp->neut_species[i].fnew;
          }
        }

        gyrokinetic_multib_forward_euler(app, tcurr+dt, dt, fin, fout, fin_neut, fout_neut, &st);

        if (st.dt_actual < dt) {

          // Recalculate the field.
          for (int b=0; b<app->num_local_blocks; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            int li_charged = b * app->num_species;
            for (int i=0; i<app->num_species; ++i)
              fin[li_charged+i] = sbapp->species[i].f;
          }
          gyrokinetic_multib_calc_field(app, tcurr, fin);

          // collect stats
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_2_fail += 1;

          dt = st.dt_actual;
          state = RK_STAGE_1; // restart from stage 1

        } 
        else {
          for (int b=0; b<app->num_local_blocks; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            for (int i=0; i<app->num_species; ++i) {
              array_combine(sbapp->species[i].f1,
                3.0/4.0, sbapp->species[i].f, 1.0/4.0, sbapp->species[i].fnew, &sbapp->species[i].local_ext);
            }
            for (int i=0; i<app->num_neut_species; ++i) {
              if (!sbapp->neut_species[i].info.is_static) {
                array_combine(sbapp->neut_species[i].f1,
                  3.0/4.0, sbapp->neut_species[i].f, 1.0/4.0, sbapp->neut_species[i].fnew, &sbapp->neut_species[i].local_ext);
              }
            }
          }

          for (int b=0; b<app->num_local_blocks; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            int li_charged = b * app->num_species;
            int li_neut = b * app->num_neut_species;
            // Compute the fields and apply BCs.
            for (int i=0; i<app->num_species; ++i) {
              fout[li_charged+i] = sbapp->species[i].f1;
            }
            for (int i=0; i<app->num_neut_species; ++i) {
              fout_neut[li_neut+i] = sbapp->neut_species[i].f1;
            }
          }
          gyrokinetic_multib_calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int b=0; b<app->num_local_blocks; ++b) {
          struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
          int li_charged = b * app->num_species;
          int li_neut = b * app->num_neut_species;
          for (int i=0; i<app->num_species; ++i) {
            fin[li_charged+i] = sbapp->species[i].f1;
            fout[li_charged+i] = sbapp->species[i].fnew;
          }
          for (int i=0; i<app->num_neut_species; ++i) {
	    fin_neut[li_neut+i] = sbapp->neut_species[i].f1;
	    fout_neut[li_neut+i] = sbapp->neut_species[i].fnew;
          }
        }

        gyrokinetic_multib_forward_euler(app, tcurr+dt/2, dt, fin, fout, fin_neut, fout_neut, &st);

        if (st.dt_actual < dt) {
          // Recalculate the field.
          for (int b=0; b<app->num_local_blocks; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            int li_charged = b * app->num_species;
            for (int i=0; i<app->num_species; ++i)
              fin[li_charged+i] = sbapp->species[i].f;
          }
          gyrokinetic_multib_calc_field(app, tcurr, fin);

          // collect stats
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_3_fail += 1;

          dt = st.dt_actual;
          state = RK_STAGE_1; // restart from stage 1

          app->stat.nstage_2_fail += 1;
        }
        else {
          for (int b=0; b<app->num_local_blocks; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            for (int i=0; i<app->num_species; ++i) {
              array_combine(sbapp->species[i].f1,
                1.0/3.0, sbapp->species[i].f, 2.0/3.0, sbapp->species[i].fnew, &sbapp->species[i].local_ext);
              gkyl_array_copy_range(sbapp->species[i].f, sbapp->species[i].f1, &sbapp->species[i].local_ext);
            }
            for (int i=0; i<app->num_neut_species; ++i) {
              if (!sbapp->neut_species[i].info.is_static) {
                array_combine(sbapp->neut_species[i].f1,
                  1.0/3.0, sbapp->neut_species[i].f, 2.0/3.0, sbapp->neut_species[i].fnew, &sbapp->neut_species[i].local_ext);
                gkyl_array_copy_range(sbapp->neut_species[i].f, sbapp->neut_species[i].f1, &sbapp->neut_species[i].local_ext);
              }
            }
          }

//          if (app->enforce_positivity) {
//            // Apply positivity shift if requested.
//            int elc_idx = -1;
//            gkyl_array_clear(app->ps_delta_m0_ions, 0.0);
//            for (int i=0; i<app->num_species; ++i) {
//              struct gk_species *gks = &app->species[i];
//
//              // Copy f so we can calculate the moments of the change later. 
//              gkyl_array_set(gks->fnew, -1.0, gks->f);
//
//              // Shift each species.
//              gkyl_positivity_shift_gyrokinetic_advance(gks->pos_shift_op, &app->local, &gks->local,
//                gks->f, gks->m0.marr, gks->ps_delta_m0);
//
//              // Accumulate the shift density of all ions:
//              if (gks->info.charge > 0.0)
//                gkyl_array_accumulate(app->ps_delta_m0_ions, 1.0, gks->ps_delta_m0);
//              else if (gks->info.charge < 0.0) 
//                elc_idx = i;
//            }
//
//            // Rescale each species to enforce quasineutrality.
//            for (int i=0; i<app->num_species; ++i) {
//              struct gk_species *gks = &app->species[i];
//              if (gks->info.charge > 0.0) {
//                struct gk_species *gkelc = &app->species[elc_idx];
//                gkyl_positivity_shift_gyrokinetic_quasineutrality_scale(gks->pos_shift_op, &app->local, &gks->local,
//                  gks->ps_delta_m0, app->ps_delta_m0_ions, gkelc->ps_delta_m0, gks->m0.marr, gks->f);
//              }
//              else {
//                gkyl_positivity_shift_gyrokinetic_quasineutrality_scale(gks->pos_shift_op, &app->local, &gks->local,
//                  gks->ps_delta_m0, gks->ps_delta_m0, app->ps_delta_m0_ions, gks->m0.marr, gks->f);
//              }
//
//              gkyl_array_accumulate(gks->fnew, 1.0, gks->f);
//            }
//          }

          // Compute the fields and apply BCs
          for (int b=0; b<app->num_local_blocks; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            int li_charged = b * app->num_species;
            int li_neut = b * app->num_neut_species;
            for (int i=0; i<app->num_species; ++i) {
              fout[li_charged+i] = sbapp->species[i].f;
            }
            for (int i=0; i<app->num_neut_species; ++i) {
              fout_neut[li_neut+i] = sbapp->neut_species[i].f;
            }
          }
          gyrokinetic_multib_calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

          state = RK_COMPLETE;
        }
        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }

  return st;
}


