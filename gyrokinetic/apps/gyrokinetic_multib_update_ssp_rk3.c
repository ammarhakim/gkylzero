#include <gkyl_gyrokinetic_multib_priv.h>

static void
gyrokinetic_multib_forward_euler(struct gkyl_gyrokinetic_multib_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], struct gkyl_array *fout[], 
  const struct gkyl_array **bflux_in[], struct gkyl_array **bflux_out[], 
  const struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[], 
  const struct gkyl_array **bflux_in_neut[], struct gkyl_array **bflux_out_neut[], 
  struct gkyl_update_status *st)
{
  struct timespec wst_fe = gkyl_wall_clock();
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
      &bflux_out[li_charged], &fin_neut[li_neut], &fout_neut[li_neut], &bflux_out_neut[li_neut], st);
    dtmin = fmin(dtmin, st->dt_actual);
  }

  struct timespec wtm = gkyl_wall_clock();
  // Compute minimum time-step across all processors.
  double dtmin_local = dtmin, dtmin_global;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dtmin_local, &dtmin_global);
  st->dt_actual = dtmin_global;
  app->stat.dfdt_dt_reduce_tm += gkyl_time_diff_now_sec(wtm);

  struct timespec wst = gkyl_wall_clock();
  // Complete update of distribution functions.
  double dta = st->dt_actual;
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int li_charged = b * app->num_species;
    int li_neut = b * app->num_neut_species;
    for (int i=0; i<app->num_species; ++i) {
      struct gk_species *gks = &sbapp->species[i];
      gk_species_step_f(gks, fout[li_charged+i], dta, fin[li_charged+i]);
      gk_species_bflux_step_f(sbapp, &gks->bflux, bflux_out[li_charged+i], dta, bflux_in[li_charged+i]);
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      struct gk_neut_species *gkns = &sbapp->neut_species[i];
      gk_neut_species_step_f(gkns, fout_neut[li_charged+i], dta, fin_neut[li_charged+i]);
      gk_neut_species_bflux_step_f(sbapp, &gkns->bflux, bflux_out_neut[li_neut+i], dta, bflux_in_neut[li_neut+i]);
    }
  }

  app->stat.fwd_euler_step_f_tm += gkyl_time_diff_now_sec(wst);
  app->stat.fwd_euler_tm += gkyl_time_diff_now_sec(wst_fe);
}

struct gkyl_update_status
gyrokinetic_multib_update_ssp_rk3(struct gkyl_gyrokinetic_multib_app* app, double dt0)
{
  // Take time-step using the RK3 method. Also sets the status object
  // which has the actual and suggested dts used. These can be different
  // from the actual time-step.
  int ns_charged = app->num_species;
  int ns_neut = app->num_neut_species;
  int nblocks_local = app->num_local_blocks;

  const struct gkyl_array *fin[ns_charged * nblocks_local];
  struct gkyl_array *fout[ns_charged * nblocks_local];
  const struct gkyl_array **bflux_in[ns_charged * nblocks_local];
  struct gkyl_array **bflux_out[ns_charged * nblocks_local];

  const struct gkyl_array *fin_neut[ns_neut * nblocks_local];
  struct gkyl_array *fout_neut[ns_neut * nblocks_local];
  const struct gkyl_array **bflux_in_neut[ns_neut * nblocks_local];
  struct gkyl_array **bflux_out_neut[ns_neut * nblocks_local];

  struct gkyl_update_status st = { .success = true };

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        for (int b=0; b<nblocks_local; ++b) {
          struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
          int li_charged = b * ns_charged;
          int li_neut = b * ns_neut;
          for (int i=0; i<ns_charged; ++i) {
            struct gk_species *gks = &sbapp->species[i];
            fin[li_charged+i] = gks->f;
            fout[li_charged+i] = gks->f1;
            // Boundary fluxes.
            gk_species_bflux_clear(sbapp, &gks->bflux, gks->bflux.f, 0.0);
            gk_species_bflux_clear(sbapp, &gks->bflux, gks->bflux.f1, 0.0);
            bflux_in[li_charged+i] = (const struct gkyl_array **)gks->bflux.f;
            bflux_out[li_charged+i] = gks->bflux.f1;
          }
          for (int i=0; i<ns_neut; ++i) {
            struct gk_neut_species *gkns = &sbapp->neut_species[i];
            fin_neut[li_neut+i] = gkns->f;
	    fout_neut[li_neut+i] = gkns->f1;
            // Boundary fluxes.
            gk_neut_species_bflux_clear(sbapp, &gkns->bflux, gkns->bflux.f, 0.0);
            gk_neut_species_bflux_clear(sbapp, &gkns->bflux, gkns->bflux.f1, 0.0);
            bflux_in_neut[li_neut+i] = (const struct gkyl_array **)gkns->bflux.f;
            bflux_out_neut[li_neut+i] = gkns->bflux.f1;
          }
        }

        gyrokinetic_multib_forward_euler(app, tcurr, dt, fin, fout, bflux_in, bflux_out,
          fin_neut, fout_neut, bflux_in_neut, bflux_out_neut, &st);
        dt = st.dt_actual;

        for (int b=0; b<nblocks_local; ++b) {
          struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
          for (int i=0; i<ns_charged; ++i) {
            struct gk_species *gks = &sbapp->species[i];
            // Compute moment of f_old to later compute moment of df/dt.
            // Do it before the fields are updated, but after dt is calculated.
            gk_species_calc_int_mom_dt(sbapp, gks, dt, gks->fdot_mom_old);
          }

          // Compute field energy divided by dt for energy balance diagnostics.
          gk_field_calc_energy_dt(sbapp, sbapp->field, dt, sbapp->field->em_energy_red_old);
        }

        // Compute the fields and apply BCs.
        gyrokinetic_multib_calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int b=0; b<nblocks_local; ++b) {
          struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
          int li_charged = b * ns_charged;
          int li_neut = b * ns_neut;
          for (int i=0; i<ns_charged; ++i) {
            struct gk_species *gks = &sbapp->species[i];
            fin[li_charged+i] = gks->f1;
            fout[li_charged+i] = gks->fnew;
            // Boundary fluxes.
            bflux_in[li_charged+i] = (const struct gkyl_array **)gks->bflux.f1;
            bflux_out[li_charged+i] = gks->bflux.fnew;
          }
          for (int i=0; i<ns_neut; ++i) {
            struct gk_neut_species *gkns = &sbapp->neut_species[i];
	    fin_neut[li_neut+i] = gkns->f1;
	    fout_neut[li_neut+i] = gkns->fnew;
            // Boundary fluxes.
            bflux_in_neut[li_neut+i] = (const struct gkyl_array **)gkns->bflux.f1;
            bflux_out_neut[li_neut+i] = gkns->bflux.fnew;
          }
        }

        gyrokinetic_multib_forward_euler(app, tcurr+dt, dt, fin, fout, bflux_in, bflux_out,
          fin_neut, fout_neut, bflux_in_neut, bflux_out_neut, &st);

        if (st.dt_actual < dt) {

          // Recalculate the field.
          for (int b=0; b<nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            int li_charged = b * ns_charged;
            for (int i=0; i<ns_charged; ++i)
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
          struct timespec wst = gkyl_wall_clock();
          for (int b=0; b<nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            for (int i=0; i<ns_charged; ++i) {
	      struct gk_species *gks = &sbapp->species[i]; 
	      gk_species_combine(gks, gks->f1, 3.0/4.0, gks->f, 1.0/4.0, gks->fnew, &gks->local_ext);
              gk_species_bflux_combine(sbapp, &gks->bflux, gks->bflux.f1,
                3.0/4.0, gks->bflux.f, 1.0/4.0, gks->bflux.fnew);
	    }
            for (int i=0; i<ns_neut; ++i) {
	      struct gk_neut_species *gkns = &sbapp->neut_species[i]; 
	      gk_neut_species_combine(gkns, gkns->f1, 3.0/4.0, gkns->f, 1.0/4.0, gkns->fnew, &gkns->local_ext);
              gk_neut_species_bflux_combine(sbapp, &gkns->bflux, gkns->bflux.f1,
                3.0/4.0, gkns->bflux.f, 1.0/4.0, gkns->bflux.fnew);
            }
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

          for (int b=0; b<nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            int li_charged = b * ns_charged;
            int li_neut = b * ns_neut;
            // Compute the fields and apply BCs.
            for (int i=0; i<ns_charged; ++i) {
              fout[li_charged+i] = sbapp->species[i].f1;
            }
            for (int i=0; i<ns_neut; ++i) {
              fout_neut[li_neut+i] = sbapp->neut_species[i].f1;
            }
          }
          gyrokinetic_multib_calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int b=0; b<nblocks_local; ++b) {
          struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
          int li_charged = b * ns_charged;
          int li_neut = b * ns_neut;
          for (int i=0; i<ns_charged; ++i) {
	    struct gk_species *gks = &sbapp->species[i]; 
            fin[li_charged+i] = gks->f1;
            fout[li_charged+i] = gks->fnew;
            // Boundary fluxes.
            bflux_in[li_charged+i] = (const struct gkyl_array **)gks->bflux.f1;
            bflux_out[li_charged+i] = gks->bflux.fnew;
          }
          for (int i=0; i<ns_neut; ++i) {
            struct gk_neut_species *gkns = &sbapp->neut_species[i];
	    fin_neut[li_neut+i] = sbapp->neut_species[i].f1;
	    fout_neut[li_neut+i] = sbapp->neut_species[i].fnew;
            // Boundary fluxes.
            bflux_in_neut[li_neut+i] = (const struct gkyl_array **)gkns->bflux.f1;
            bflux_out_neut[li_neut+i] = gkns->bflux.fnew;
          }
        }

        gyrokinetic_multib_forward_euler(app, tcurr+dt/2, dt, fin, fout, bflux_in, bflux_out,
          fin_neut, fout_neut, bflux_in_neut, bflux_out_neut, &st);

        if (st.dt_actual < dt) {
          // Recalculate the field.
          for (int b=0; b<nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            int li_charged = b * ns_charged;
            for (int i=0; i<ns_charged; ++i)
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
          struct timespec wst = gkyl_wall_clock();
          for (int b=0; b<nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            for (int i=0; i<ns_charged; ++i) {
	      struct gk_species *gks = &sbapp->species[i]; 
              // Step f.
	      gk_species_combine(gks, gks->f1, 1.0/3.0, gks->f, 2.0/3.0, gks->fnew, &gks->local_ext);
	      gk_species_copy_range(gks, gks->f, gks->f1, &gks->local_ext);
              // Step boundary fluxes.
              gk_species_bflux_combine(sbapp, &gks->bflux, gks->bflux.f1,
                1.0/3.0, gks->bflux.f, 2.0/3.0, gks->bflux.fnew);
              gk_species_bflux_copy(sbapp, &gks->bflux, gks->bflux.f, gks->bflux.f1);
              gk_species_bflux_calc_voltime_integrated_mom(sbapp, gks, &gks->bflux, tcurr);
              gk_species_bflux_scale(sbapp, &gks->bflux, gks->bflux.f, 1.0/dt);
	    }
            for (int i=0; i<ns_neut; ++i) {
	      struct gk_neut_species *gkns = &sbapp->neut_species[i]; 
	      gk_neut_species_combine(gkns, gkns->f1, 1.0/3.0, gkns->f, 2.0/3.0, gkns->fnew, &gkns->local_ext);
	      gk_neut_species_copy_range(gkns, gkns->f, gkns->f1, &gkns->local_ext);
              // Step boundary fluxes.
              gk_neut_species_bflux_combine(sbapp, &gkns->bflux, gkns->bflux.f1,
                1.0/3.0, gkns->bflux.f, 2.0/3.0, gkns->bflux.fnew);
              gk_neut_species_bflux_copy(sbapp, &gkns->bflux, gkns->bflux.f, gkns->bflux.f1);
              gk_neut_species_bflux_calc_voltime_integrated_mom(sbapp, gkns, &gkns->bflux, tcurr);
              gk_neut_species_bflux_scale(sbapp, &gkns->bflux, gkns->bflux.f, 1.0/dt);
            }
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

//          if (app->enforce_positivity) {
//            // Apply positivity shift if requested.
//            int elc_idx = -1;
//            gkyl_array_clear(app->ps_delta_m0_ions, 0.0);
//            for (int i=0; i<ns_charged; ++i) {
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
//            for (int i=0; i<ns_charged; ++i) {
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
          for (int b=0; b<nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
            int li_charged = b * ns_charged;
            int li_neut = b * ns_neut;
            for (int i=0; i<ns_charged; ++i) {
              fout[li_charged+i] = sbapp->species[i].f;
            }
            for (int i=0; i<ns_neut; ++i) {
              fout_neut[li_neut+i] = sbapp->neut_species[i].f;
            }
          }
          gyrokinetic_multib_calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

          for (int b=0; b<nblocks_local; ++b) {
            struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
	    for (int i=0; i<ns_charged; ++i) {
              struct gk_species *gks = &sbapp->species[i];
              // Compute moment of f_new to compute moment of df/dt.
              // Need to do it after the fields are updated.
              gk_species_calc_int_mom_dt(sbapp, gks, dt, gks->fdot_mom_new);
            }

            // Compute field energy divided by dt for energy balance diagnostics.
            gk_field_calc_energy_dt(sbapp, sbapp->field, dt, sbapp->field->em_energy_red_new);
          }

          state = RK_COMPLETE;
        }
        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }

  return st;
}


