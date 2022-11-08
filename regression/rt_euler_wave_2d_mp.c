#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

static inline double sq(double x) { return x * x; }

struct euler_ctx {
  double gas_gamma; // gas constant
};

void
evalEulerInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  // See https://ammar-hakim.org/sj/je/je22/je22-euler-2d.html#smooth-periodic-problem
  struct euler_ctx *app = ctx;
  double gas_gamma = app->gas_gamma;

  double x = xn[0], y = xn[1];
  double rho = 1 + 0.2*sin(M_PI*(x+y));
  double u = 1.0, v = -0.5;
  double pr = 1.0;
  
  fout[0] = rho;
  fout[1] = rho*u; fout[2] = rho*v; fout[3] = 0.0;
  fout[4] = pr/(gas_gamma-1) + 0.5*rho*(u*u+v*v);
}

struct euler_ctx
euler_ctx(void)
{
  return (struct euler_ctx) { .gas_gamma = 1.4 };
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 50);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 50);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  } 
  struct euler_ctx ctx = euler_ctx(); // context for init functions

  // equation object
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_inew( &(struct gkyl_wv_euler_inp) {
      .gas_gamma = ctx.gas_gamma,
      .rp_type = WV_EULER_RP_LAX 
    }
  );

  struct gkyl_moment_species fluid = {
    .name = "euler",

    .equation = euler,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalEulerInit,
  };

  // VM app
  struct gkyl_moment app_inp = {
    .name = "euler_wave_2d_mp",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { 2.0, 2.0 },
    .cells = { NX, NY },
    .cfl_frac = 0.9,
    .scheme_type = GKYL_MOMENT_MP,
    .mp_recon = app_args.mp_recon,    

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },

    .num_species = 1,
    .species = { fluid },
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 4.0;

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  gkyl_moment_app_write(app, tcurr, 0);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
    if (tcurr+dt > tend) dt = tend-tcurr;
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
  gkyl_wv_eqn_release(euler);
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);

  return 0;
}
