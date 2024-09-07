#include <gkyl_gyrokinetic_multib_priv.h>
// Take a forward Euler step with the suggested time-step dt. This may
// not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by
// stability. The actual time-step and dt_suggested are returned in
// the status object.
void
gyrokinetic_multib_forward_euler(struct gkyl_gyrokinetic_multib_app* mbapp, double tcurr, double dt,
  const struct gkyl_array *fin[], struct gkyl_array *fout[], 
  const struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[], 
  struct gkyl_update_status *st, struct gkyl_update_status *sb_st)
{
  mbapp->stat.nfeuler += 1;

  double dtmin = DBL_MAX;
  for (int lbidx = 0; lbidx < mbapp->num_local_blocks; lbidx++) {
    gkyl_gyrokinetic_app *app = mbapp->singleb_apps[lbidx];
    int lin_idx = lbidx*mbapp->num_species;
    gyrokinetic_rhs(app, tcurr, dt, &fin[lin_idx], &fout[lin_idx], &fin_neut[lin_idx], &fout_neut[lin_idx], &sb_st[lbidx]);
    dtmin = fmin(dtmin, sb_st[lbidx].dt_actual);
  }



  double dt_max_rel_diff = 0.01;
  // Check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // Compute minimum time-step across all processors.
  double dtmin_local = dtmin, dtmin_global;
  gkyl_comm_allreduce_host(mbapp->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dtmin_local, &dtmin_global);
  dtmin = dtmin_global;
  
  // Don't take a time-step larger that input dt.
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // Complete update of distribution functions.
  for (int lbidx = 0; lbidx < mbapp->num_local_blocks; lbidx++) {
    gkyl_gyrokinetic_app *app = mbapp->singleb_apps[lbidx];
    int lin_idx = lbidx*mbapp->num_species;
    for (int i=0; i<app->num_species; ++i) {
      gkyl_array_accumulate(gkyl_array_scale(fout[lin_idx+i], dta), 1.0, fin[lin_idx+i]);
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      if (!app->neut_species[i].info.is_static) {
        gkyl_array_accumulate(gkyl_array_scale(fout_neut[lin_idx+i], dta), 1.0, fin_neut[lin_idx+i]);
      }
    }
  }

}


