#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_eqn_type.h>
#include <gkyl_vlasov.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

struct euler_ctx {
  double gas_gamma; // gas constant
};

void
evalEulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct euler_ctx *app = ctx;
  double gas_gamma = app->gas_gamma;
  pcg64_random_t rng = gkyl_pcg64_init(0);
  double x = xn[0], y = xn[1];

  double rho = 1.0, vx = 0.5, vy = 0.0, pr = 2.5;
  if (fabs(y)<0.25) {
    rho = 2.0;
    vx = -0.5;
  }

  double alpha = 1.0e-2;
  double k = 2.0*M_PI;
  for (int i=0; i<16; ++i) {
    for (int j=0; j<16; ++j) {
      vx += alpha*gkyl_pcg64_rand_double(&rng)*sin(i*k*x + j*k*y + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
      vy += alpha*gkyl_pcg64_rand_double(&rng)*sin(i*k*x + j*k*y + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
    }
  }
  
  fout[0] = rho;
  fout[1] = rho*vx; fout[2] = rho*vy; fout[3] = 0.0;
  fout[4] = 0.5*rho*(vx*vx+vy*vy) + pr/(gas_gamma-1);
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 64);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 64);  

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct euler_ctx ctx = euler_ctx(); // context for init functions

  // Euler equation object
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);

  struct gkyl_vlasov_fluid_species f = {
    .name = "euler",

    .ctx = &ctx,
    .equation = euler,
    .init = evalEulerInit,
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "dg_euler_kh_2d",

    .cdim = 2,
    .lower = { -0.5, -0.5 },
    .upper = { 0.5, 0.5 }, 
    .cells = { NX, NY },
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .cfl_frac = 0.9,

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },   

    .num_species = 0,
    .species = {},

    .num_fluid_species = 1,
    .fluid_species = { f },

    .skip_field = true, //TODO: double check this

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 10.0;
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
      fprintf(stderr, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    step += 1;
  }

  gkyl_vlasov_app_write(app, tcurr, 1);
  gkyl_vlasov_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  // simulation complete, free resources
  gkyl_wv_eqn_release(euler);
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
