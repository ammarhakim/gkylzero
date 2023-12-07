#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

#define SQ(x) ((x)*(x))

struct euler_ctx {
  double gas_gamma; // gas constant
  double inlet_xlen; // lenght of inlet 
  double bump_xlen; // x extent of bump
  double height; // height of channel
  double R; // radius of bump
};

void
evalEulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct euler_ctx *app = ctx;
  double halfx = app->inlet_xlen+0.5*app->bump_xlen;
  double gas_gamma = app->gas_gamma;

  double rho = 1.0, pr = 1.0;
  double cs = sqrt(gas_gamma*pr/rho);
  double u = 2.0*cs;

  if (xn[0] < -halfx*(1-0.25)) {
    fout[0] = rho;
    fout[1] = rho*u;  fout[2] = 0.0; fout[3] = 0.0;
    fout[4] = pr/(gas_gamma-1) + 0.5*rho*u*u;
  }
  else {
    fout[0] = rho*1e-5;
    fout[1] = 0.0;  fout[2] = 0.0; fout[3] = 0.0;
    fout[4] = pr*1e-5/(gas_gamma-1);
  }
}

// map (r,theta) -> (x,y)
void
mapc2p(double t, const double *xc_in, double* GKYL_RESTRICT xp, void *ctx)
{
  struct euler_ctx *app = ctx;
  
  double bump_xlen = app->bump_xlen;
  double halfx = app->inlet_xlen+0.5*bump_xlen;
  double R = app->R, H = app->height;
  
  double xc = xc_in[0], yc = xc_in[1];
  double zeta_min = sqrt(SQ(R) - SQ(0.5*bump_xlen));
  
  xp[0] = xc; xp[1] = yc;

  if (fabs(xc) < 0.5*bump_xlen) {
    double eta = xc;
    double zeta = sqrt(SQ(R) - SQ(eta));
    double yb = zeta-zeta_min;
    xp[1] = (H-yb)*yc/H + yb;
  }
}

struct euler_ctx
euler_ctx(void)
{
  return (struct euler_ctx) {
    .gas_gamma = 1.4,
    .inlet_xlen = 1.0,
    .bump_xlen = 2.0,
    .height = 2.0,
    .R = 5.0
  };
}

void
write_data(struct gkyl_tm_trigger *iot, const gkyl_moment_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr))
    gkyl_moment_app_write(app, tcurr, iot->curr-1);
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 150);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 75);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct euler_ctx ctx = euler_ctx(); // context for init functions

  double halfx = ctx.inlet_xlen+0.5*ctx.bump_xlen;

  // equation object
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(ctx.gas_gamma, false);

  struct gkyl_moment_species fluid = {
    .name = "euler",

    .equation = euler,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalEulerInit,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_COPY },
  };

  // VM app
  struct gkyl_moment app_inp = {
    .name = "euler_bump_in_channel",

    .ndim = 2,
    // grid in computational space
    .lower = { -halfx, 0.0 },
    .upper = { halfx, ctx.height },
    .cells = { NX, NY },

    .mapc2p = mapc2p, // mapping of computational to physical space
    .c2p_ctx = &ctx,
    
    .cfl_frac = 0.9,

    .num_species = 1,
    .species = { fluid },
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 10.0;
  int nframe = 1;
  
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

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

    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_moment_app_write(app, tcurr, 1000);
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
