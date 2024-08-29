#include <gkyl_gyrokinetic_multib_priv.h>
// Take time-step using the RK3 method. Also sets the status object
// which has the actual and suggested dts used. These can be different
// from the actual time-step.
static struct gkyl_update_status
gyrokinetic_multib_update_ssp_rk3(struct gkyl_gyrokinetic_multib_app* mbapp, double dt0)
{
  const struct gkyl_array **fin[mbapp->num_local_blocks];
  struct gkyl_array **fout[mbapp->num_local_blocks];
  const struct gkyl_array **fin_neut[mbapp->num_local_blocks];
  struct gkyl_array **fout_neut[mbapp->num_local_blocks];
  struct gkyl_update_status st = { .success = true };
  struct gkyl_update_status sb_st[mbapp->num_local_blocks];

  for (int i =0; i<mbapp->num_local_blocks; i++) {
    fin[i] = gkyl_malloc(sizeof(struct gkyl_array*)*mbapp->num_species);
    fout[i] = gkyl_malloc(sizeof(struct gkyl_array*)*mbapp->num_species);
    fin_neut[i] = gkyl_malloc(sizeof(struct gkyl_array*)*mbapp->num_neut_species);
    fout_neut[i] = gkyl_malloc(sizeof(struct gkyl_array*)*mbapp->num_neut_species);
    sb_st[i].success = true;
  }

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = mbapp->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        for (int lbidx = 0; lbidx < mbapp->num_local_blocks; lbidx++) {
          struct gkyl_gyrokinetic_app *app = mbapp->singleb_apps[lbidx];
          for (int i=0; i<app->num_species; ++i) {
            fin[lbidx][i] = app->species[i].f;
            fout[lbidx][i] = app->species[i].f1;
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            fin_neut[lbidx][i] = app->neut_species[i].f;
            if (!app->neut_species[i].info.is_static) {
              fout_neut[lbidx][i] = app->neut_species[i].f1;
            }
          }
        }

        // Forward euler will loop over blocks and calculate the rhs, 
        // do a reduction using the big communicator to get dt
        // Then do loop over local blocks and do the accumulation
        gyrokinetic_multib_forward_euler(mbapp, tcurr, dt, fin, fout, fin_neut, fout_neut, &st, sb_st);

        // TO DO Compute the fields and apply BCs.

        dt = st.dt_actual;
        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int lbidx = 0; lbidx < mbapp->num_local_blocks; lbidx++) {
          struct gkyl_gyrokinetic_app *app = mbapp->singleb_apps[lbidx];
          for (int i=0; i<app->num_species; ++i) {
            fin[lbidx][i] = app->species[i].f1;
            fout[lbidx][i] = app->species[i].fnew;
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            if (!app->neut_species[i].info.is_static) {
              fin_neut[lbidx][i] = app->neut_species[i].f1;
              fout_neut[lbidx][i] = app->neut_species[i].fnew;
            }
            else {
              fin_neut[lbidx][i] = app->neut_species[i].f;
            }
          }
        }

        gyrokinetic_multib_forward_euler(mbapp, tcurr+dt, dt, fin, fout, fin_neut, fout_neut, &st, sb_st);

        if (st.dt_actual < dt) {

          // TO DO: Recalculate the field.
          for (int lbidx = 0; lbidx < mbapp->num_local_blocks; lbidx++) {
            struct gkyl_gyrokinetic_app *app = mbapp->singleb_apps[lbidx];
            for (int i=0; i<app->num_species; ++i)
              fin[lbidx][i] = app->species[i].f;
          }
          //calc_field(app, tcurr, fin);

          // collect stats
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          mbapp->stat.stage_2_dt_diff[0] = fmin(mbapp->stat.stage_2_dt_diff[0],
            dt_rel_diff);
          mbapp->stat.stage_2_dt_diff[1] = fmax(mbapp->stat.stage_2_dt_diff[1],
            dt_rel_diff);
          mbapp->stat.nstage_2_fail += 1;

          dt = st.dt_actual;
          state = RK_STAGE_1; // restart from stage 1

        } 
        else {
          for (int lbidx = 0; lbidx < mbapp->num_local_blocks; lbidx++) {
            struct gkyl_gyrokinetic_app *app = mbapp->singleb_apps[lbidx];
            for (int i=0; i<app->num_species; ++i) {
              array_combine(app->species[i].f1,
                3.0/4.0, app->species[i].f, 1.0/4.0, app->species[i].fnew, &app->species[i].local_ext);
            }
            for (int i=0; i<app->num_neut_species; ++i) {
              if (!app->neut_species[i].info.is_static) {
                array_combine(app->neut_species[i].f1,
                  3.0/4.0, app->neut_species[i].f, 1.0/4.0, app->neut_species[i].fnew, &app->neut_species[i].local_ext);
              }
            }
          }

          // TO DO: Compute the fields and apply BCs.
          for (int lbidx = 0; lbidx < mbapp->num_local_blocks; lbidx++) {
            struct gkyl_gyrokinetic_app *app = mbapp->singleb_apps[lbidx];
            for (int i=0; i<app->num_species; ++i) {
              fout[lbidx][i] = app->species[i].f1;
            }
            for (int i=0; i<app->num_neut_species; ++i) {
              fout_neut[lbidx][i] = app->neut_species[i].f1;
            }
          }
          //calc_field_and_apply_bc(mbapp, tcurr, fout, fout_neut);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int lbidx = 0; lbidx < mbapp->num_local_blocks; lbidx++) {
          struct gkyl_gyrokinetic_app *app = mbapp->singleb_apps[lbidx];
          for (int i=0; i<app->num_species; ++i) {
            fin[lbidx][i] = app->species[i].f1;
            fout[lbidx][i] = app->species[i].fnew;
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            if (!app->neut_species[i].info.is_static) {
              fin_neut[lbidx][i] = app->neut_species[i].f1;
              fout_neut[lbidx][i] = app->neut_species[i].fnew;
            }
            else {
              fin_neut[lbidx][i] = app->neut_species[i].f;
            }          
          }
        }

        gyrokinetic_multib_forward_euler(mbapp, tcurr+dt/2, dt, fin, fout, fin_neut, fout_neut, &st, sb_st);

        if (st.dt_actual < dt) {
          // TO DO: Recalculate the field.
          for (int lbidx = 0; lbidx < mbapp->num_local_blocks; lbidx++) {
            struct gkyl_gyrokinetic_app *app = mbapp->singleb_apps[lbidx];
            for (int i=0; i<app->num_species; ++i)
              fin[lbidx][i] = app->species[i].f;
          }
          //calc_field(mbapp, tcurr, fin);

          // collect stats
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          mbapp->stat.stage_3_dt_diff[0] = fmin(mbapp->stat.stage_3_dt_diff[0],
            dt_rel_diff);
          mbapp->stat.stage_3_dt_diff[1] = fmax(mbapp->stat.stage_3_dt_diff[1],
            dt_rel_diff);
          mbapp->stat.nstage_3_fail += 1;

          dt = st.dt_actual;
          state = RK_STAGE_1; // restart from stage 1

          mbapp->stat.nstage_2_fail += 1;
        }
        else {
          for (int lbidx = 0; lbidx < mbapp->num_local_blocks; lbidx++) {
            struct gkyl_gyrokinetic_app *app = mbapp->singleb_apps[lbidx];
            for (int i=0; i<app->num_species; ++i) {
              array_combine(app->species[i].f1,
                1.0/3.0, app->species[i].f, 2.0/3.0, app->species[i].fnew, &app->species[i].local_ext);
              gkyl_array_copy_range(app->species[i].f, app->species[i].f1, &app->species[i].local_ext);
            }
            for (int i=0; i<app->num_neut_species; ++i) {
              if (!app->neut_species[i].info.is_static) {
                array_combine(app->neut_species[i].f1,
                  1.0/3.0, app->neut_species[i].f, 2.0/3.0, app->neut_species[i].fnew, &app->neut_species[i].local_ext);
                gkyl_array_copy_range(app->neut_species[i].f, app->neut_species[i].f1, &app->neut_species[i].local_ext);
              }
            }
          }

          // Implement positivity shift if requested.
          for (int lbidx = 0; lbidx < mbapp->num_local_blocks; lbidx++) {
            struct gkyl_gyrokinetic_app *app = mbapp->singleb_apps[lbidx];
            for (int i=0; i<app->num_species; ++i) {
              struct gk_species *gks = &app->species[i];
              if (gks->enforce_positivity)
                gkyl_positivity_shift_gyrokinetic_advance(gks->pos_shift_op, &gks->local, &app->local, gks->f, gks->ps_intmom_grid);
            }
          }

          // TO DO: Compute the fields and apply BCs
          for (int lbidx = 0; lbidx < mbapp->num_local_blocks; lbidx++) {
            struct gkyl_gyrokinetic_app *app = mbapp->singleb_apps[lbidx];
            for (int i=0; i<app->num_species; ++i) {
              fout[lbidx][i] = app->species[i].f;
            }
            for (int i=0; i<app->num_neut_species; ++i) {
              fout_neut[lbidx][i] = app->neut_species[i].f;
            }
          }
          //calc_field_and_apply_bc(mbapp, tcurr, fout, fout_neut);

          state = RK_COMPLETE;
        }
        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }

  return st;
}

