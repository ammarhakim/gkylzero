#include <gkyl_pkpm_priv.h>

// Strang splitting of the locally implicit fluid-Maxwell source update and the SSP-RK3
// step of the hyperbolic terms (Wang et al. JCP 2020): a half step of the sources, the
// explicit step, and a second half step. The step the explicit update will take is only
// known afterwards (it is bounded by the CFL condition, and the suggested step of the
// very first update is the user's initial dt), so the first half step uses the step of
// the previous update and the second one covers the remainder: the sources are advanced
// by exactly the step taken, and the split is symmetric whenever the step is steady. The
// first update of a run is a plain Lie split (no previous step). The implicit midpoint
// rule is time-reversible, so a negative remainder is well defined.
struct gkyl_update_status
pkpm_update_op_split(gkyl_pkpm_app *app, double dt0)
{
  double dt_half = app->op_split_dt_prev > 0.0 ? 0.5 * fmin(dt0, app->op_split_dt_prev) : 0.0;
  if (dt_half > 0.0) {
    pkpm_fluid_em_coupling_update(app, app->pkpm_em, app->tcurr, dt_half);
  }

  struct gkyl_update_status st = pkpm_update_explicit_ssp_rk3(app, dt0);

  pkpm_fluid_em_coupling_update(
    app, app->pkpm_em, app->tcurr + st.dt_actual, st.dt_actual - dt_half
  );
  app->op_split_dt_prev = st.dt_actual;

  return st;
}
