#include <gkyl_vlasov_priv.h>

// Implicit contributions
void
vlasov_implicit_contribution(gkyl_vlasov_app* app, double dt,
  const struct gkyl_array *fin[], struct gkyl_array *fout[])
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
    double dt1 = vm_species_rhs_implicit(app, &app->species[i], fin[i], fout[i], dt);
    dtmin = fmin(dtmin, dt1);
  }

  // complete update of distribution function
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout[i], dt), 1.0, fin[i]);
    vm_species_apply_bc(app, &app->species[i], fout[i]);
  }
}
