#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <rt_arg_parse.h>

// map (r,theta) -> (x,y)
void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double r = xc[0], th = xc[1];
  xp[0] = r*cos(th); xp[1] = r*sin(th);
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], phi = xn[1];
  int m = 2;

  double w = 1.19318673737701; // frequency of mode
  double a = 1.0, b = 0.9904672582498093;  

  double Ez_r = a*jn(m,r*w) + b*yn(m,r*w);
  double Ez = Ez_r*cos(m*phi);

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = Ez;
  fout[3] = 0.0;
  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;  
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 32);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 32*6);  

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // VM app
  struct gkyl_moment app_inp = {
    .name = "maxwell_ann_wg",

    .ndim = 2,
    // grid in computational space
    .lower = { 2.0, 0.0 },
    .upper = { 5.0, 2*GKYL_PI },
    .cells = { NX, NY },

    .mapc2p = mapc2p, // mapping of computational to physical space

    .num_periodic_dir = 1,
    .periodic_dirs = { 1 },

    .cfl_frac = 1.0,

    .field = {
      .epsilon0 = 1.0, .mu0 = 1.0,
      .evolve = 1,
      .limiter = GKYL_NO_LIMITER,
      .init = evalFieldInit,

      .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  double w = 1.19318673737701; // frequency of mode
  double tperiod = 2*M_PI/w;

  // start, end and initial time-step
  double tcurr = 0.0, tend = 2*tperiod;

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  gkyl_moment_app_write(app, tcurr, 0);
  gkyl_moment_app_calc_field_energy(app, 0.0);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);

    gkyl_moment_app_calc_field_energy(app, tcurr);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    step += 1;
  }

  gkyl_moment_app_write(app, tcurr, 1);
  gkyl_moment_app_write_field_energy(app);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  // simulation complete, free resources
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);
  
  return 0;
}
