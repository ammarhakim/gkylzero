#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

static inline double sq(double x) { return x * x; }

struct wg_ctx {
  double Lx, Ly;
  int n, m;
};

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct wg_ctx *wc = ctx;
  
  // assumes epsilon0 = mu0 = c = 1.0
  double Lx = wc->Lx, Ly = wc->Ly;
  int n = wc->n, m = wc->m;
  double x = xn[0], y = xn[1];

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = sin(n*M_PI*x/Lx)*sin(m*M_PI*y/Ly);
  fout[3] = 0.0;
  fout[4] = 0.0;
  fout[5] = cos(n*M_PI*x/Lx)*cos(m*M_PI*y/Ly);
  fout[6] = 0.0;
  fout[7] = 0.0;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 10*7);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 10*5);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct wg_ctx wc =  {
    .Lx = 7.0, .Ly = 5.0,
    .n = 7, .m = 5
  };

  struct gkyl_moment app_inp =
  {.name = "maxwell_wg_2d",

   .ndim = 2,
   .lower = {0.0, 0.0},
   .upper = {wc.Lx, wc.Ly},
   .cells = {NX, NY},

   .num_periodic_dir = 0,
   .periodic_dirs = {},
   .cfl_frac = 1.0,

   .field = {
     .epsilon0 = 1.0,
     .mu0 = 1.0,
     .evolve = 1,
     .limiter = GKYL_NO_LIMITER,

     .init = evalFieldInit,
     .ctx = &wc,
     
     .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
     .bcy = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL }
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  double omega = sqrt( sq(wc.n*M_PI/wc.Lx) + sq(wc.m*M_PI/wc.Ly) );
  double tperiod = 2*M_PI/omega; // time-period
  // start, end and initial time-step
  double tcurr = 0.0, tend = 10*tperiod;

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
    dt = fmin(status.dt_suggested, tend-tcurr);

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
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);
  
  return 0;
}
