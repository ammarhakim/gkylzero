#include <gkyl_gyrokinetic_priv.h>

static void
gyrokinetic_forward_euler(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], struct gkyl_array *fout[], 
  const struct gkyl_array **bflux_in[], struct gkyl_array **bflux_out[], 
  const struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[], 
  struct gkyl_update_status *st)
{
  // Take a forward Euler step with the suggested time-step dt. This may
  // not be the actual time-step taken. However, the function will never
  // take a time-step larger than dt even if it is allowed by
  // stability. The actual time-step and dt_suggested are returned in
  // the status object.
  app->stat.nfeuler += 1;

  // Compute the time rate of change of the distributions, df/dt.
  gyrokinetic_rhs(app, tcurr, dt, fin, fout, bflux_out, fin_neut, fout_neut, st);

  // Complete update of distribution functions.
  double dta = st->dt_actual;
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_step_f(gks, fout[i], dta, fin[i]);
    if (app->species[i].info.boundary_flux_diagnostics) {
      gk_species_bflux_forward_euler(app, &gks->bflux_diag, bflux_out[i], dta, bflux_in[i]);
    }
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    gk_neut_species_step_f(gkns, fout_neut[i], dta, fin_neut[i]);
  }

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

  struct gkyl_update_status st = { .success = true };

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        if (app->fdot_diagnostics) {
          for (int i=0; i<app->num_species; ++i) {
            struct gk_species *gk_s = &app->species[i];
            // Compute moment of f_old to later compute moment of df/dt.
            // Need to do it before the fields are updated.
            gk_species_moment_calc(&gk_s->integ_moms, gk_s->local, app->local, gk_s->f); 
            gkyl_array_copy(gk_s->fdot_mom_old, gk_s->integ_moms.marr);
          }

          gkyl_array_integrate_advance(app->field->calc_em_energy, app->field->phi_smooth, 
            1.0/dt, app->field->es_energy_fac, &app->local, &app->local, app->field->em_energy_red_old);
        }
        for (int i=0; i<app->num_species; ++i) {
          struct gk_species *gk_s = &app->species[i];
          fin[i] = gk_s->f;
          fout[i] = gk_s->f1;
          if (gk_s->info.boundary_flux_diagnostics) {
            gk_species_bflux_clear(app, &gk_s->bflux_diag, gk_s->bflux_diag.f, 0.0);
            gk_species_bflux_clear(app, &gk_s->bflux_diag, gk_s->bflux_diag.f1, 0.0);
            bflux_in[i] = (const struct gkyl_array **)gk_s->bflux_diag.f;
            bflux_out[i] = gk_s->bflux_diag.f1;
          }
        }
        for (int i=0; i<app->num_neut_species; ++i) {
          fin_neut[i] = app->neut_species[i].f;
	  fout_neut[i] = app->neut_species[i].f1;
        }

        gyrokinetic_forward_euler(app, tcurr, dt, fin, fout, bflux_in, bflux_out, fin_neut, fout_neut, &st);
        // Compute the fields and apply BCs.
        gyrokinetic_calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

        dt = st.dt_actual;
        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
          if (app->species[i].info.boundary_flux_diagnostics) {
            bflux_in[i] = (const struct gkyl_array **)app->species[i].bflux_diag.f1;
            bflux_out[i] = app->species[i].bflux_diag.fnew;
          }
        }
        for (int i=0; i<app->num_neut_species; ++i) {
	  fin_neut[i] = app->neut_species[i].f1;
	  fout_neut[i] = app->neut_species[i].fnew;
        }

        gyrokinetic_forward_euler(app, tcurr+dt, dt, fin, fout, bflux_in, bflux_out, fin_neut, fout_neut, &st);

        if (st.dt_actual < dt) {

          // Recalculate the field.
          for (int i=0; i<app->num_species; ++i)
            fin[i] = app->species[i].f;
          gyrokinetic_calc_field(app, tcurr, fin);

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
          for (int i=0; i<app->num_species; ++i) {
	    struct gk_species *gks = &app->species[i];
	    gk_species_combine(gks, gks->f1, 3.0/4.0, gks->f, 1.0/4.0, gks->fnew, &gks->local_ext);
            if (gks->info.boundary_flux_diagnostics) {
              gk_species_bflux_combine_range(app, gks, &gks->bflux_diag, gks->bflux_diag.f1,
                3.0/4.0, gks->bflux_diag.f, 1.0/4.0, gks->bflux_diag.fnew, &app->local_ext);
            }
          }
          for (int i=0; i<app->num_neut_species; ++i) {
	    struct gk_neut_species *gkns = &app->neut_species[i];
	    gk_neut_species_combine(gkns, gkns->f1, 3.0/4.0, gkns->f, 1.0/4.0, gkns->fnew, &gkns->local_ext);
          }

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
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
          if (app->species[i].info.boundary_flux_diagnostics) {
            bflux_in[i] = (const struct gkyl_array **)app->species[i].bflux_diag.f1;
            bflux_out[i] = app->species[i].bflux_diag.fnew;
          }
        }
        for (int i=0; i<app->num_neut_species; ++i) {
	  fin_neut[i] = app->neut_species[i].f1;
	  fout_neut[i] = app->neut_species[i].fnew;
        }

        gyrokinetic_forward_euler(app, tcurr+dt/2, dt, fin, fout, bflux_in, bflux_out, fin_neut, fout_neut, &st);

        if (st.dt_actual < dt) {
          // Recalculate the field.
          for (int i=0; i<app->num_species; ++i)
            fin[i] = app->species[i].f;
          gyrokinetic_calc_field(app, tcurr, fin);

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
          for (int i=0; i<app->num_species; ++i) {
	    struct gk_species *gks = &app->species[i];
	    gk_species_combine(gks, gks->f1, 1.0/3.0, gks->f, 2.0/3.0, gks->fnew, &gks->local_ext);
	    gk_species_copy_range(gks, gks->f, gks->f1, &gks->local_ext);
            if (gks->info.boundary_flux_diagnostics) {
              gk_species_bflux_combine_range(app, gks, &gks->bflux_diag, gks->bflux_diag.f1,
                1.0/3.0, gks->bflux_diag.f, 2.0/3.0, gks->bflux_diag.fnew, &app->local_ext);
              gk_species_bflux_copy_range(app, gks, &gks->bflux_diag,
                gks->bflux_diag.f, gks->bflux_diag.f1, &app->local_ext);
              gk_species_bflux_scale(app, &gks->bflux_diag, gks->bflux_diag.f, 1.0/dt);
            }
          }
          for (int i=0; i<app->num_neut_species; ++i) {
	    struct gk_neut_species *gkns = &app->neut_species[i];
	    gk_neut_species_combine(gkns, gkns->f1, 1.0/3.0, gkns->f, 2.0/3.0, gkns->fnew, &gkns->local_ext);
	    gk_neut_species_copy_range(gkns, gkns->f, gkns->f1, &gkns->local_ext);
          }

          if (app->enforce_positivity) {
            // Apply positivity shift if requested.
	    gkyl_array_clear(app->ps_delta_m0_ions, 0.0);
            gkyl_array_clear(app->ps_delta_m0_elcs, 0.0);
            for (int i=0; i<app->num_species; ++i) {
              struct gk_species *gks = &app->species[i];

              // Copy f so we can calculate the moments of the change later. 
              gkyl_array_set(gks->fnew, -1.0, gks->f);

              // Shift each species.
              gkyl_positivity_shift_gyrokinetic_advance(gks->pos_shift_op, &app->local, &gks->local,
                gks->f, gks->m0.marr, gks->ps_delta_m0);

	      // Accumulate the shift density of all like-species:
              gkyl_array_accumulate(gks->ps_delta_m0s_tot, 1.0, gks->ps_delta_m0);
            }

            // Rescale each species to enforce quasineutrality.
	    for (int i=0; i<app->num_species; ++i) {
              struct gk_species *gks = &app->species[i];
              gkyl_positivity_shift_gyrokinetic_quasineutrality_scale(gks->pos_shift_op, &app->local, &gks->local,
                gks->ps_delta_m0, gks->ps_delta_m0s_tot, gks->ps_delta_m0r_tot, gks->m0.marr, gks->f);

              gkyl_array_accumulate(gks->fnew, 1.0, gks->f);
            }
          }

          // Compute the fields and apply BCs
          for (int i=0; i<app->num_species; ++i) {
            fout[i] = app->species[i].f;
          }
          for (int i=0; i<app->num_neut_species; ++i) {
            fout_neut[i] = app->neut_species[i].f;
          }
          gyrokinetic_calc_field_and_apply_bc(app, tcurr, fout, fout_neut);

          if (app->fdot_diagnostics) {
	    for (int i=0; i<app->num_species; ++i) {
              struct gk_species *gk_s = &app->species[i];
              // Compute moment of f_new to compute moment of df/dt.
              // Need to do it after the fields are updated.
              gk_species_moment_calc(&gk_s->integ_moms, gk_s->local, app->local, gk_s->f); 
              gkyl_array_set(gk_s->fdot_mom_new, 1.0/dt, gk_s->integ_moms.marr);
              gkyl_array_accumulate(gk_s->fdot_mom_new, -1.0/dt, gk_s->fdot_mom_old);
            }

            gkyl_array_integrate_advance(app->field->calc_em_energy, app->field->phi_smooth, 
              1.0/dt, app->field->es_energy_fac, &app->local, &app->local, app->field->em_energy_red_new);
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

