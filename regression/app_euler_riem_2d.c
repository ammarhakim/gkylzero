#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <app_arg_parse.h>

struct euler_ctx {
  double gas_gamma; // gas constant
};

void
evalEulerInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  struct euler_ctx *app = ctx;
  double gas_gamma = app->gas_gamma;

  double sloc = 0.8;
  double rho, u, v, pr;

  double x = xn[0], y = xn[1];

  double upLeft[] = {0.0, 0.3, 0.5323, 1.206, 0.0};
  double upRight[] = {0.0, 1.5, 1.5, 0.0, 0.0};
  double loLeft[] = {0.0, 0.029, 0.138, 1.206, 1.206};
  double loRight[] = {0.0, 0.3, 0.5323, 0.0, 1.206};

  if (y>sloc) {
    if (x<sloc) {
      pr = upLeft[1];
      rho = upLeft[2];
      u = upLeft[3];
      v = upLeft[4];
    }
    else {
      pr = upRight[1];
      rho = upRight[2];
      u = upRight[3];
      v = upRight[4];
    }
  }
  else {
    if (x<sloc) {
      pr = loLeft[1];
      rho = loLeft[2];
      u = loLeft[3];
      v = loLeft[4];
    }
    else {
      pr = loRight[1];
      rho = loRight[2];
      u = loRight[3];
      v = loRight[4];
    }
  }
  fout[0] = rho;
  fout[1] = rho*u; fout[2] = rho*v; fout[3] = 0.0;
  fout[4] = 0.5*rho*(u*u+v*v) + pr/(gas_gamma-1);
}

struct euler_ctx
euler_ctx(void)
{
  return (struct euler_ctx) { .gas_gamma = 1.4 };
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = get_parse_app_args(argc, argv); 
  struct euler_ctx ctx = euler_ctx(); // context for init functions

  // equation object
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(ctx.gas_gamma);

  struct gkyl_moment_species fluid = {
    .name = "euler",

    .equation = euler,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalEulerInit,
  };

  // VM app
  struct gkyl_moment app_inp = {
    .name = "euler_riem_2d",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { 1.0, 1.0 }, 
    .cells = { 200, 200 },

    .num_species = 1,
    .species = { fluid },
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 0.8;

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

  // write stats to file
  FILE *fp = fopen("euler_riem_2d_stat.json", "a");
  if (fp) {
    gkyl_moment_stat_write_json(fp, stat);
    fclose(fp);
  }
  
  return 0;
}
