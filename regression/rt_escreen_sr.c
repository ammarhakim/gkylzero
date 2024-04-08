#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct escreen_inp {
  double tend;
  double charge, mass; // for electrons
  int conf_cells, vel_cells;
  double vel_extents[2];
};

struct escreen_ctx {
  double T; // electron thermal velocity
  double vdrift; // drift velocity
};

static inline double sq(double x) { return x*x; }

void
evalDistFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct escreen_ctx *app = ctx;
  double x = xn[0], p = xn[1];
  double T = app->T, vdrift = app->vdrift;

  // modified Bessel function of the second kind evaluated for T = mc^2 (K_2(1))
  double K_2 = 1.6248388986351774828107073822838437146593935281628733843345054697;
  // modified Bessel function of the second kind evaluated for T = 0.1 mc^2 (K_2(10))
  //double K_2 = 0.0000215098170069327687306645644239671272492068461808732468335569;
  // modified Bessel function of the second kind evaluated for T = 0.04 mc^2 (K_2(25))
  //double K_2 = 3.7467838080691090570137658745889511812329380156362352887017e-12;
  // Lorentz factor for drift velocity
  double gamma = 1.0/sqrt(1 - vdrift*vdrift);

  double n = 1.0;
  double mc2_T = 1.0/T;

  double fv = n/K_2*exp(-mc2_T*gamma*(sqrt(1 + p*p) - vdrift*p));
  fout[0] = fv;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct escreen_ctx *app = ctx;
  double x = xn[0];

  double E_x = 100.0;

  fout[0] = E_x; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

struct escreen_ctx
create_default_ctx(void)
{
  return (struct escreen_ctx) {
    .T = 1.0,
    .vdrift = 0.0,
  };
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_vlasov_app_write(app, tcurr, iot->curr-1);
    gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, iot->curr-1);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 16);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 1024);  

  struct escreen_ctx ctx;
  
  ctx = create_default_ctx(); // context for init functions

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_SR,
    .charge = -1.0,
    .mass = 1.0,
    .lower = { -768.0 },
    .upper = { 768.0 }, 
    .cells = { VX },

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFunc,
      .ctx_func = &ctx,
    },

    // source is the same as initial condition
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .projection = {
        .proj_id = GKYL_PROJ_FUNC,
        .func = evalDistFunc,
        .ctx_func = &ctx,
      },
    },

    .num_diag_moments = 2,
    .diag_moments = { "M0", "M1i" },
  };

  // positrons
  struct gkyl_vlasov_species pos = {
    .name = "pos",
    .model_id = GKYL_MODEL_SR,
    .charge = 1.0,
    .mass = 1.0,
    .lower = { -768.0 },
    .upper = { 768.0 }, 
    .cells = { VX },

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFunc,
      .ctx_func = &ctx,
    },

    // source is the same as initial condition
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .projection = {
        .proj_id = GKYL_PROJ_FUNC,
        .func = evalDistFunc,
        .ctx_func = &ctx,
      },
    },

    .num_diag_moments = 2,
    .diag_moments = { "M0", "M1i" },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "escreen_sr",
    
    .cdim = 1, .vdim = 1,
    .lower = { 0.0 },
    .upper = { 1.0 },
    .cells = { NX },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, pos },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };
  
  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 25.0;
  double dt = tend-tcurr;
  int nframe = 250;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

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
    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_vlasov_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  // simulation complete, free objects
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
