#include <gkyl_gyrokinetic_priv.h>
// Take a forward Euler step with the suggested time-step dt. This may
// not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by
// stability. The actual time-step and dt_suggested are returned in
// the status object.
void
gyrokinetic_forward_euler(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], struct gkyl_array *fout[], 
  const struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[], 
  struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;

  double dtmin = DBL_MAX;

  // Compute necessary moments and boundary corrections for collisions.
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      gk_species_lbo_moms(app, &app->species[i], 
        &app->species[i].lbo, fin[i]);
    }
    else if (app->species[i].collision_id == GKYL_BGK_COLLISIONS) {
      gk_species_bgk_moms(app, &app->species[i], 
        &app->species[i].bgk, fin[i]);
    }
  }

  // Compute necessary moments for cross-species collisions.
  // Needs to be done after self-collisions moments, so separate loop over species.
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) { 
      if (app->species[i].lbo.num_cross_collisions) {
        gk_species_lbo_cross_moms(app, &app->species[i], 
          &app->species[i].lbo, fin[i]);        
      }
    }
    else if (app->species[i].collision_id == GKYL_BGK_COLLISIONS) {
      if (app->species[i].bgk.num_cross_collisions) {
        gk_species_bgk_cross_moms(app, &app->species[i], 
          &app->species[i].bgk, fin[i]);        
      }
    }
    // Compute reaction rates (e.g., ionization, recombination, or charge exchange).
    if (app->species[i].has_reactions) {
      gk_species_react_cross_moms(app, &app->species[i], 
        &app->species[i].react, fin[i], fin, fin_neut);
    }
    if (app->species[i].has_neutral_reactions) {
      gk_species_react_cross_moms(app, &app->species[i], 
        &app->species[i].react_neut, fin[i], fin, fin_neut);
    }
    // Compute necessary drag coefficients for radiation operator.
    if (app->species[i].radiation_id == GKYL_GK_RADIATION) {
      gk_species_radiation_moms(app, &app->species[i], 
        &app->species[i].rad, fin, fin_neut);
    }
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    // Compute reaction cross moments (e.g., ionization, recombination, or charge exchange).
    if (app->neut_species[i].has_neutral_reactions) {
      gk_neut_species_react_cross_moms(app, &app->neut_species[i], 
        &app->neut_species[i].react_neut, fin, fin_neut);
    }
  }

  // Compute RHS of Gyrokinetic equation.
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    double dt1 = gk_species_rhs(app, s, fin[i], fout[i]);
    dtmin = fmin(dtmin, dt1);

    // Compute and store (in the ghost cell of of out) the boundary fluxes.
    // NOTE: this overwrites ghost cells that may be used for sourcing.
    if (app->update_field || app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN)
      gk_species_bflux_rhs(app, s, &s->bflux, fin[i], fout[i]);
  }

  // Compute RHS of neutrals.
  for (int i=0; i<app->num_neut_species; ++i) {
    double dt1 = gk_neut_species_rhs(app, &app->neut_species[i], fin_neut[i], fout_neut[i]);
    dtmin = fmin(dtmin, dt1);
  }

  // Compute plasma source term.
  // Done here as the RHS update for all species should be complete before
  // in case we are using boundary fluxes as a component of our source function
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].source_id) {
      gk_species_source_rhs(app, &app->species[i], 
        &app->species[i].src, fin[i], fout[i]);
    }
  }

  // Compute neutral source term.
  // Done here as the RHS update for all species should be complete before
  // in case we are using boundary fluxes as a component of our source function.
  for (int i=0; i<app->num_neut_species; ++i) {
    if (app->neut_species[i].source_id) {
      gk_neut_species_source_rhs(app, &app->neut_species[i], 
        &app->neut_species[i].src, fin_neut[i], fout_neut[i]);
    }
  }

  double dt_max_rel_diff = 0.01;
  // Check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // Compute minimum time-step across all processors.
  double dtmin_local = dtmin, dtmin_global;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dtmin_local, &dtmin_global);
  dtmin = dtmin_global;
  
  // Don't take a time-step larger that input dt.
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // Complete update of distribution functions.
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout[i], dta), 1.0, fin[i]);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    if (!app->neut_species[i].info.is_static) {
      gkyl_array_accumulate(gkyl_array_scale(fout_neut[i], dta), 1.0, fin_neut[i]);
    }
  }

}


