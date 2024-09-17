#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

static inline double sq(double x) { return x*x; }

struct wg_ctx {
  double Lx, Ly;
  int n, m;
};

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
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

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 3*7);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 3*5);

  struct wg_ctx wc =  {
    .Lx = 7.0, .Ly = 5.0,
    .n = 7, .m = 5
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .init = evalFieldFunc,
    .ctx = &wc,

    .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
    .bcy = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL }
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "dg_maxwell_wg_2d",

    .cdim = 2, .vdim = 0,
    .lower = { 0.0, 0.0 },
    .upper = { wc.Lx, wc.Ly },
    .cells = { NX, NY },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 0,

    .num_species = 0,
    .species = {  },
    .field = field,
 
    .parallelism = {
      .use_gpu = app_args.use_gpu,
    },
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  double omega = sqrt( sq(wc.n*M_PI/wc.Lx) + sq(wc.m*M_PI/wc.Ly) );
  double tperiod = 2*M_PI/omega; // time-period
  // start, end and initial time-step
  double tcurr = 0.0, tend = 10*tperiod;
  double dt = tend-tcurr;

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  
  gkyl_vlasov_app_write(app, tcurr, 0);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = fmin(status.dt_suggested, tend-tcurr);
    step += 1;
  }

  gkyl_vlasov_app_write(app, tcurr, 1);
  gkyl_vlasov_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  // simulation complete, free app
  gkyl_vlasov_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of forward-Euler calls %ld\n", stat.nfeuler);
  printf("Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    printf("Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    printf("Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  printf("Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  printf("Species RHS calc took %g secs\n", stat.species_rhs_tm);
  printf("Field RHS calc took %g secs\n", stat.field_rhs_tm);
  printf("Current evaluation and accumulate took %g secs\n", stat.current_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
