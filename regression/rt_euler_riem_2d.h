#pragma once

#include <gkyl_wv_euler.h>

// Everything is supposed to be private

#include <math.h>
#include <stdio.h>

#include <rt_euler_riem_2d.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

struct euler_ctx {
  double gas_gamma; // gas constant
};

static void
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

static const char*
get_sim_name(enum gkyl_wv_euler_rp rp_type)
{
  switch (rp_type) {
    case WV_EULER_RP_ROE:
      return "euler_riem_2d_roe";
      break;
    case WV_EULER_RP_HLLC:
      return "euler_riem_2d_hllc";
      break;
    case WV_EULER_RP_LAX:
      return "euler_riem_2d_lax";
      break;
  }
}

// Run regression test using specified RP type  
static int
rt_euler_riem_2d_run(int argc, char **argv, enum gkyl_wv_euler_rp rp_type)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 200);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 200);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  } 
  struct euler_ctx ctx = euler_ctx(); // context for init functions

  // equation object
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_inew(
    &(struct gkyl_wv_euler_inp) {
      .gas_gamma = ctx.gas_gamma,
      .rp_type = rp_type
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
    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { 1.0, 1.0 }, 
    .cells = { NX, NY },

    .num_species = 1,
    .species = { fluid },
  };
  strcpy(app_inp.name, get_sim_name(rp_type));

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 0.8;

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  gkyl_moment_app_write(app, tcurr, 0);
  gkyl_moment_app_calc_integrated_mom(app, tcurr);  

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);

    gkyl_moment_app_calc_integrated_mom(app, tcurr);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    step += 1;
  }

  gkyl_moment_app_write(app, tcurr, 1);
  gkyl_moment_app_write_integrated_mom(app);
    
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
