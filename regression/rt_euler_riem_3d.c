/* 3d Riemann test taken from
 * https://www.ammar-hakim.org/sj/sims/s408/s408-riemann-euler-3d.html
 *
 * This will be used to benchmark the 2d axisymmetric simulation using mapped
 * grid. See also Langseth and LeVeque, section 3.2 for the original paper.
 */

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
   // See Langseth and LeVeque, section 3.2
  struct euler_ctx *app = ctx;
  double gas_gamma = app->gas_gamma;

  double rhoi = 1.0, pri =5.0;
  double rho0 = 1.0, pr0 =1.0;
  double rho, pr;

  double x = xn[0], y = xn[1], z = xn[2];

  double r = sqrt(sq(x) + sq(y) + sq(z-0.4));
  if (r<0.2)
  {
     rho = rhoi; pr = pri;
  } else {
     rho = rho0; pr = pr0;
  }

  fout[0] = rho;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = pr/(gas_gamma-1);
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 37);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 37);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], 25);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  } 
  struct euler_ctx ctx = euler_ctx(); // context for init functions

  // equation object
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(ctx.gas_gamma);

  struct gkyl_moment_species fluid = {
    .name = "euler",

    .equation = euler,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalEulerInit,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
    .bcy = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
    .bcz = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
  };

  // VM app
  struct gkyl_moment app_inp = {
    .name = "euler_riem_3d",

    .ndim = 3,
    .lower = { 0.0, 0.0, 0.0 },
    .upper = { 1.5, 1.5, 1.0 }, 
    .cells = { NX, NY, NZ },
    .cfl_frac = 0.9,

    .num_species = 1,
    .species = { fluid },
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 0.7;

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
