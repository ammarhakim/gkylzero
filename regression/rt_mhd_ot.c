#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_moment_prim_mhd.h>
#include <gkyl_util.h>
#include <gkyl_wv_mhd.h>
#include <rt_arg_parse.h>

struct mhd_ctx {
  enum gkyl_wv_mhd_div_constraint divergence_constraint;
  double gas_gamma; // gas constant
};

void
calcq(double gas_gamma, const double pv[8], double q[8])
{
  double rho = pv[0], u = pv[1], v = pv[2], w = pv[3], pr = pv[4];
  q[0] = rho;
  q[1] = rho*u; q[2] = rho*v; q[3] = rho*w;
  q[5] = pv[5]; q[6] = pv[6]; q[7] = pv[7];  // B field
  double pb = 0.5*(pv[5]*pv[5]+pv[6]*pv[6]+pv[7]*pv[7]); // magnetic pressure
  q[4] = pr/(gas_gamma-1) + 0.5*rho*(u*u+v*v+w*w) + pb;
}

void
evalMhdInit(double t, const double* GKYL_RESTRICT xn,
  double* GKYL_RESTRICT fout, void *ctx)
{
  struct mhd_ctx *app = ctx;
  double x = xn[0], y =xn[1];
  double gas_gamma = app->gas_gamma;

  double rho = 25/(36*M_PI);
  double p = 5/(12*M_PI);
  double vx = sin(2*M_PI*y);
  double vy = -sin(2*M_PI*x);
  double vz = 0;
  double B0 = 1/sqrt(4*M_PI);
  double Bx = B0*sin(2*M_PI*y);
  double By = B0*sin(4*M_PI*x);
  double Bz = 0;
  double v[8] = {rho, vx, vy, vz, p, Bx, By, Bz};

  gkyl_mhd_cons_vars(gas_gamma, v, fout);

  if (app->divergence_constraint==GKYL_MHD_DIVB_GLM)
    fout[8] = 0; // divB correction potential
}

struct mhd_ctx
mhd_ctx(void)
{
  return (struct mhd_ctx) { .gas_gamma = 5./3. };
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 128);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 128);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct mhd_ctx ctx = {
    .gas_gamma = 5./3.,
    .divergence_constraint = GKYL_MHD_DIVB_EIGHT_WAVES,
  };

  // equation object
  const struct gkyl_wv_mhd_inp inp = {
    .gas_gamma = ctx.gas_gamma,
    .divergence_constraint = ctx.divergence_constraint,
    .glm_ch = 1.0, // initial value; will be updated with max speed in each step
    .glm_alpha = 0.4, // passed to source
  };
  struct gkyl_wv_eqn *mhd = gkyl_wv_mhd_new(&inp);

  struct gkyl_moment_species fluid = {
    .name = "mhd",

    .equation = mhd,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalMhdInit,
  };

  // VM app
  struct gkyl_moment app_inp = {
    .name = "mhd_ot_glm",

    .ndim = 2,
    .lower = { -0.5, -0.5 },
    .upper = { 0.5, 0.5 }, 
    .cells = { NX, NY },

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },
    .cfl_frac = 0.6,

    .num_species = 1,
    .species = { fluid },
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 0.5;
  int nframe = 5;

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
