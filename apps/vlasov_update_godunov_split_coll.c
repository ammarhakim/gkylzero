#include <gkyl_vlasov_priv.h>

// Take time-step using the RK3 method for the explicit advective comps.
// Use the actual timestep used to update the implicit BGK using a 
// forward euler method
struct gkyl_update_status
vlasov_update_godunov_split_coll(gkyl_vlasov_app* app, double dt0)
{

  // Take the RK3 timstep
  struct gkyl_update_status st = vlasov_update_ssp_rk3(app,dt0);

  // Take the implicit timestep
  vlasov_update_implicit_coll(app,st.dt_actual);

  return st;
}