#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_advect.h>
#include <rt_arg_parse.h>

#define sq(x) ((x) * (x))

void evalInit(double t, const double *GKYL_RESTRICT xn,
              double *GKYL_RESTRICT fout, void *ctx) {
  // Suresh and Huynh (1997) JCP 136,83-99, eq. (4.1), fig. 4.1, 4.2
  double x = xn[0], f;
  if (-0.8 <= x && x <= -0.6)
    f = exp(-log(2.0) * sq(x + 0.7) / 0.0009);
  else if (-0.4 <= x && x <= -0.2)
    f = 1.0;
  else if (0 <= x && x <= 0.2)
    f = 1.0 - fabs(10.0 * (x - 0.1));
  else if (0.4 <= x && x <= 0.6)
    f = sqrt(1.0 - 100.0 * sq(x - 0.5));
  else
    f = 0.0;

  fout[0] = f;
}

int main(int argc, char **argv) {
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 200);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  double c = 1.0;
  struct gkyl_wv_eqn *advect = gkyl_wv_advect_new(c);

  struct gkyl_moment_species fluid = {
      .name = "advect",

      .equation = advect,
      .evolve = true,
      .init = evalInit,
  };

  struct gkyl_moment app_inp = {
      .name = "advect_sine",

      .ndim = 1,
      .lower = {-1.0},
      .upper = {1.0},
      .cells = {NX},

      .num_periodic_dir = 1,
      .periodic_dirs = {0},
      .cfl_frac = 0.4,

      .num_species = 1,
      .species = {fluid},
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 20.0;

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  gkyl_moment_app_write(app, tcurr, 0);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);

    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    step += 1;
  }

  gkyl_moment_app_write(app, tcurr, 1);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  // simulation complete, free resources
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Total updates took %g secs\n", stat.total_tm);

  return 0;
}
