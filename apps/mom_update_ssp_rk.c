#include <gkyl_moment_priv.h>

// internal function that takes a single time-step using a single-step
// Strang-split scheme
struct gkyl_update_status
moment_update_ssp_rk3(gkyl_moment_app* app, double dt0)
{
  double dt = dt0;
  double dt_suggested = dt0;
  
  return (struct gkyl_update_status) {
    .success = true,
    .dt_actual = dt,
    .dt_suggested = dt_suggested,
  };  
}
