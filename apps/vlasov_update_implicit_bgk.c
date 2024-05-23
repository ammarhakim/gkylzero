#include <gkyl_vlasov_priv.h>

// Take time-step using the RK3 method for the explicit advective comps.
// Use the actual timestep used to update 
void
vlasov_update_implicit_bgk(gkyl_vlasov_app* app, double dt0)
{
  int ns = app->num_species;  
  const struct gkyl_array *fin[ns];
  struct gkyl_array *fout[ns];

  for (int i=0; i<ns; ++i) {
    fin[i] = app->species[i].f;
    fout[i] = app->species[i].f1;
  }

  // implicit BGK contributions
  if (app->has_implicit_bgk_scheme){
    vlasov_implicit_contribution(app, dt0, fin, fout);
  }
  
  for (int i=0; i<ns; ++i) {
    gkyl_array_copy_range(app->species[i].f, app->species[i].f1, &app->species[i].local_ext);
  };
}