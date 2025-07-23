#include <gkyl_gyrokinetic_priv.h>

static void
gyrokinetic_forward_euler(gkyl_gyrokinetic_app* app, double tcurr, double dt,
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

  // Compute the time rate of change of the distributions, df/dt.
  gyrokinetic_rhs(app, tcurr, dt, fin, fout, bflux_out, fin_neut, fout_neut, bflux_out_neut, st);

  struct timespec wst = gkyl_wall_clock();
  // Complete update of distribution functions.
  double dta = st->dt_actual;
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_step_f(gks, fout[i], dta, fin[i]);
    gk_species_bflux_step_f(app, &gks->bflux, bflux_out[i], dta, bflux_in[i]);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_step_f(gkns, fout_neut[i], dta, fin_neut[i]);
    gk_neut_species_bflux_step_f(app, &gkns->bflux, bflux_out_neut[i], dta, bflux_in_neut[i]);
  }
  app->stat.fwd_euler_step_f_tm += gkyl_time_diff_now_sec(wst);
  app->stat.fwd_euler_tm += gkyl_time_diff_now_sec(wst_fe);

}

struct gkyl_update_status
gyrokinetic_update_ssp_rk3(gkyl_gyrokinetic_app* app, double dt0)
{
  // Take time-step using the RK3 method. Also sets the status object
  // which has the actual and suggested dts used. These can be different
  // from the actual time-step.
  const struct gkyl_array *fin[app->num_species];
  struct gkyl_array *fout[app->num_species];
  const struct gkyl_array **bflux_in[app->num_species];
  struct gkyl_array **bflux_out[app->num_species];

  const struct gkyl_array *fin_neut[app->num_neut_species];
  struct gkyl_array *fout_neut[app->num_neut_species];
  const struct gkyl_array **bflux_in_neut[app->num_neut_species];
  struct gkyl_array **bflux_out_neut[app->num_neut_species];

  struct gkyl_update_status st = { .success = true };

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        for (int i=0; i<app->num_species; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f;
          fout[i] = gks->f1;
          // Boundary fluxes.
          gk_species_bflux_clear(app, &gks->bflux, gks->bflux.f, 0.0);
          gk_species_bflux_clear(app, &gks->bflux, gks->bflux.f1, 0.0);
          bflux_in[i] = (const struct gkyl_array **)gks->bflux.f;
          bflux_out[i] = gks->bflux.f1;
        }
        for (int i=0; i<app->num_neut_species; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          fin_neut[i] = gkns->f;
	        fout_neut[i] = gkns->f1;
          // Boundary fluxes.
          gk_neut_species_bflux_clear(app, &gkns->bflux, gkns->bflux.f, 0.0);
          gk_neut_species_bflux_clear(app, &gkns->bflux, gkns->bflux.f1, 0.0);
          bflux_in_neut[i] = (const struct gkyl_array **)gkns->bflux.f;
          bflux_out_neut[i] = gkns->bflux.f1;
        }

        gyrokinetic_forward_euler(app, tcurr, dt, fin, fout, bflux_in, bflux_out,
          fin_neut, fout_neut, bflux_in_neut, bflux_out_neut, &st);
        dt = st.dt_actual;

        for (int i=0; i<app->num_species; ++i) {
          struct gk_species *gks = &app->species[i];
          // Compute moment of f_old to later compute moment of df/dt.
          // Do it before the fields are updated, but after dt is calculated.
          gk_species_calc_int_mom_dt(app, gks, dt, gks->fdot_mom_old);
        }

        // Compute field energy divided by dt for energy balance diagnostics.
        gk_field_calc_energy_dt(app, app->field, dt, app->field->em_energy_red_old);

        // Compute the fields and apply BCs.
        gyrokinetic_calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int i=0; i<app->num_species; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f1;
          fout[i] = gks->fnew;
          // Boundary fluxes.
          bflux_in[i] = (const struct gkyl_array **)gks->bflux.f1;
          bflux_out[i] = gks->bflux.fnew;
        }
        for (int i=0; i<app->num_neut_species; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          fin_neut[i] = gkns->f1;
          fout_neut[i] = gkns->fnew;
          // Boundary fluxes.
          bflux_in_neut[i] = (const struct gkyl_array **)gkns->bflux.f1;
          bflux_out_neut[i] = gkns->bflux.fnew;
        }

        gyrokinetic_forward_euler(app, tcurr+dt, dt, fin, fout, bflux_in, bflux_out,
          fin_neut, fout_neut, bflux_in_neut, bflux_out_neut, &st);

        if (st.dt_actual < dt) {

          // Recalculate the field.
          for (int i=0; i<app->num_species; ++i)
            fin[i] = app->species[i].f;
          gyrokinetic_calc_field(app, tcurr, fin);

          // Collect stats.
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_2_fail += 1;

          dt = st.dt_actual;
          state = RK_STAGE_1; // Restart from stage 1.

        } 
        else {
          struct timespec wst = gkyl_wall_clock();
          for (int i=0; i<app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            gk_species_combine(gks, gks->f1, 3.0/4.0, gks->f, 1.0/4.0, gks->fnew, &gks->local_ext);
            gk_species_bflux_combine(app, &gks->bflux, gks->bflux.f1,
              3.0/4.0, gks->bflux.f, 1.0/4.0, gks->bflux.fnew);
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            gk_neut_species_combine(gkns, gkns->f1, 3.0/4.0, gkns->f, 1.0/4.0, gkns->fnew, &gkns->local_ext);
            gk_neut_species_bflux_combine(app, &gkns->bflux, gkns->bflux.f1,
              3.0/4.0, gkns->bflux.f, 1.0/4.0, gkns->bflux.fnew);
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

          // Compute the fields and apply BCs.
          for (int i=0; i<app->num_species; ++i) {
            fout[i] = app->species[i].f1;
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            fout_neut[i] = app->neut_species[i].f1;
          }
          gyrokinetic_calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int i=0; i<app->num_species; ++i) {
          struct gk_species *gks = &app->species[i];
          fin[i] = gks->f1;
          fout[i] = gks->fnew;
          // Boundary fluxes.
          bflux_in[i] = (const struct gkyl_array **)gks->bflux.f1;
          bflux_out[i] = gks->bflux.fnew;
        }
        for (int i=0; i<app->num_neut_species; ++i) {
          struct gk_neut_species *gkns = &app->neut_species[i];
          fin_neut[i] = gkns->f1;
          fout_neut[i] = gkns->fnew;
          // Boundary fluxes.
          bflux_in_neut[i] = (const struct gkyl_array **)gkns->bflux.f1;
          bflux_out_neut[i] = gkns->bflux.fnew;
        }

        gyrokinetic_forward_euler(app, tcurr+dt/2, dt, fin, fout, bflux_in, bflux_out,
          fin_neut, fout_neut, bflux_in_neut, bflux_out_neut, &st);

        if (st.dt_actual < dt) {
          // Recalculate the field.
          for (int i=0; i<app->num_species; ++i)
            fin[i] = app->species[i].f;
          gyrokinetic_calc_field(app, tcurr, fin);

          // Collect stats.
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_3_fail += 1;

          dt = st.dt_actual;
          state = RK_STAGE_1; // Restart from stage 1.

          app->stat.nstage_2_fail += 1;
        }
        else {
          struct timespec wst = gkyl_wall_clock();
          for (int i=0; i<app->num_species; ++i) {
	          struct gk_species *gks = &app->species[i];
            // Step f.
            gk_species_combine(gks, gks->f1, 1.0/3.0, gks->f, 2.0/3.0, gks->fnew, &gks->local_ext);
            gk_species_copy_range(gks, gks->f, gks->f1, &gks->local_ext);
            // Step boundary fluxes.
            gk_species_bflux_combine(app, &gks->bflux, gks->bflux.f1,
              1.0/3.0, gks->bflux.f, 2.0/3.0, gks->bflux.fnew);
            gk_species_bflux_copy(app, &gks->bflux, gks->bflux.f, gks->bflux.f1);
            gk_species_bflux_calc_voltime_integrated_mom(app, gks, &gks->bflux, tcurr);
            gk_species_bflux_scale(app, &gks->bflux, gks->bflux.f, 1.0/dt);
          }

          for (int i=0; i<app->num_neut_species; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            gk_neut_species_combine(gkns, gkns->f1, 1.0/3.0, gkns->f, 2.0/3.0, gkns->fnew, &gkns->local_ext);
            gk_neut_species_copy_range(gkns, gkns->f, gkns->f1, &gkns->local_ext);
            // Step boundary fluxes.
            gk_neut_species_bflux_combine(app, &gkns->bflux, gkns->bflux.f1,
              1.0/3.0, gkns->bflux.f, 2.0/3.0, gkns->bflux.fnew);
            gk_neut_species_bflux_copy(app, &gkns->bflux, gkns->bflux.f, gkns->bflux.f1);
            gk_neut_species_bflux_calc_voltime_integrated_mom(app, gkns, &gkns->bflux, tcurr);
            gk_neut_species_bflux_scale(app, &gkns->bflux, gkns->bflux.f, 1.0/dt);
          }
          app->stat.time_stepper_arithmetic_tm += gkyl_time_diff_now_sec(wst);

          // Apply positivity shift if requested.
          for (int i=0; i<app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            gk_species_apply_pos_shift(app, gks);
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            struct gk_neut_species *gkns = &app->neut_species[i];
            gk_neut_species_apply_pos_shift(app, gkns);
          }

          // Enforce quasineutrality of the positivity shifts.
          gyrokinetic_pos_shift_quasineutrality(app);

          // Compute the fields and apply BCs
          for (int i=0; i<app->num_species; ++i) {
            fout[i] = app->species[i].f;
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            fout_neut[i] = app->neut_species[i].f;
          }
          gyrokinetic_calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

          for (int i=0; i<app->num_species; ++i) {
            struct gk_species *gks = &app->species[i];
            // Compute moment of f_new to compute moment of df/dt.
            // Need to do it after the fields are updated.
            gk_species_calc_int_mom_dt(app, gks, dt, gks->fdot_mom_new);

            // adapt the sources
            gk_species_source_adapt(app, gks, &gks->src, gks->lte.f_lte, tcurr);
          }

          // Compute field energy divided by dt for energy balance diagnostics.
          gk_field_calc_energy_dt(app, app->field, dt, app->field->em_energy_red_new);

          state = RK_COMPLETE;
        }
        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }

  return st;
}
