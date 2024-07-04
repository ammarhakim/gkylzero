#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_mhd.h>
#include <rt_arg_parse.h>

struct mhd_ctx {
  double gas_gamma; // gas constant
};

void
evalMhdInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mhd_ctx *app = ctx;
  double x = xn[0];
  double gas_gamma = app->gas_gamma;

  double bx = 0.75;
  double rhol = 1.0, rhor = 0.125;
  double byl = 1.0, byr = -1.0;
  double pl = 1.0, pr = 0.1;

  double rho = rhor, by = byr, p = pr;
  if (x<0.5) {
    rho = rhol;
    by = byl;
    p = pl;
  }

  fout[0] = rho;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = p/(gas_gamma-1) + 0.5*(bx*bx + by*by);
  fout[5] = bx; fout[6] = by; fout[7] = 0.0;
}

struct mhd_ctx
mhd_ctx(void)
{
  return (struct mhd_ctx) { .gas_gamma = 2.0 };
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 400);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct mhd_ctx ctx = mhd_ctx(); // context for init functions

  // equation object
  const struct gkyl_wv_mhd_inp inp = {
    .gas_gamma = ctx.gas_gamma,
    .divergence_constraint = GKYL_MHD_DIVB_NONE,
  };
  struct gkyl_wv_eqn *mhd = gkyl_wv_mhd_new(&inp);

  struct gkyl_moment_species fluid = {
    .name = "mhd",

    .equation = mhd,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalMhdInit,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };
  
  // VM app
  struct gkyl_moment app_inp = {
    .name = "mhd_brio_wu",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { 1.0 }, 
    .cells = { NX },

    .num_species = 1,
    .species = { fluid },
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 0.1;

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
  gkyl_wv_eqn_release(mhd);
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);
  
  return 0;
}
