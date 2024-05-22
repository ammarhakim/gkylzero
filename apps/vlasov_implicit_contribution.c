#include <gkyl_vlasov_priv.h>

// Implicit contributions
void
vlasov_implicit_contribution(gkyl_vlasov_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[], const struct gkyl_array *emin,
  struct gkyl_array *fout[], struct gkyl_array *fluidout[], struct gkyl_array *emout, 
  struct gkyl_update_status *st)
{
  double dtmin = DBL_MAX;


  // compute necessary moments and boundary corrections for collisions
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_BGK_COLLISIONS && app->has_implicit_bgk_scheme) {
      vm_species_bgk_moms(app, &app->species[i], 
        &app->species[i].bgk, fin[i]);
    }
  }

  // compute RHS of Vlasov equations
  for (int i=0; i<app->num_species; ++i) {
    app->species[i].bgk.implicit_step = true;
    app->species[i].bgk.dt = dt;
    double dt1 = vm_species_rhs_implicit(app, &app->species[i], fin[i], emin, fout[i], dt);
    dtmin = fmin(dtmin, dt1);
  }

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // compute minimum time-step across all processors
  double dtmin_local = dtmin, dtmin_global;
  gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dtmin_local, &dtmin_global);
  dtmin = dtmin_global;
  
  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // complete update of distribution function
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout[i], dta), 1.0, fin[i]);
    vm_species_apply_bc(app, &app->species[i], fout[i]);
  }
}
