#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gkyl_vlasov.h>
#include <app_arg_parse.h>
#include <rxi_ini.h>

struct twostream_inp {
  double tend;
  double charge, mass; // for electrons
  int conf_cells, vel_cells;
  double vel_extents[2];
};

struct twostream_ctx {
  double knumber; // wave-number
  double vth; // electron thermal velocity
  double vdrift; // drift velocity
  double perturbation;
};

static inline double sq(double x) { return x*x; }

void
evalDistFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct twostream_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double alpha = app->perturbation, k = app->knumber, vt = app->vth, vdrift = app->vdrift;

  double fv = 1/sqrt(8*M_PI*sq(vt))*(exp(-sq(v-vdrift)/(2*sq(vt)))+exp(-sq(v+vdrift)/(2*sq(vt))));
  fout[0] = (1+alpha*cos(k*x))*fv;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct twostream_ctx *app = ctx;
  double x = xn[0];
  double alpha = app->perturbation, k = app->knumber;

  double E_x = -alpha*sin(k*x)/k;
  
  fout[0] = E_x; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

struct twostream_inp
create_twostream_inp(rxi_ini_t *inp)
{
  struct twostream_inp tsinp = {
    .charge = -1.0,
    .mass = 1.0,
  };

  int read_failed = 0;
  // read from input file
  if (!rxi_ini_sget(inp, "conf-grid", "tend", "%lg", &tsinp.tend)) {
    fprintf(stderr, "Must provide 'tend' in section '[conf-grid]'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "conf-grid", "cells", "%d", &tsinp.conf_cells)) {
    fprintf(stderr, "Must provide 'cells' in section '[conf-grid]'!\n");
    read_failed = 1;
  }

  rxi_ini_sget(inp, "electrons", "charge", "%lg", &tsinp.charge);
  rxi_ini_sget(inp, "electrons", "mass", "%lg", &tsinp.mass);

  if (!rxi_ini_sget(inp, "electrons", "cells", "%d", &tsinp.vel_cells)) {
    fprintf(stderr, "Must provide 'cells' in section '[electrons]'!\n");
    read_failed = 1;
  }

  if (!rxi_ini_sget(inp, "electrons", "lower", "%lg", &tsinp.vel_extents[0])) {
    fprintf(stderr, "Must provide 'lower' in section '[electrons]'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "electrons", "upper", "%lg", &tsinp.vel_extents[1])) {
    fprintf(stderr, "Must provide 'upper' in section '[electrons]'!\n");
    read_failed = 1;
  }

  if (read_failed) {
    fprintf(stderr, "... aborting!\n");
    exit(1);
  }  

  return tsinp;
}

struct twostream_ctx
create_ctx(rxi_ini_t *inp)
{
  struct twostream_ctx ctx = {
    .perturbation = 1.0e-6
  };

  int read_failed = 0;
  if (!rxi_ini_sget(inp, "electrons", "knumber", "%lg", &ctx.knumber)) {
    fprintf(stderr, "Must provide 'knumber' in section '[electrons]'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "electrons", "vth", "%lg", &ctx.vth)) {
    fprintf(stderr, "Must provide 'vth' in section '[electrons]'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "electrons", "vdrift", "%lg", &ctx.vdrift)) {
    fprintf(stderr, "Must provide 'vdrift' in section '[electrons]'!\n");
    read_failed = 1;
  }
  rxi_ini_sget(inp, "electrons", "perturbation", "%lg", &ctx.perturbation);

  if (read_failed) {
    fprintf(stderr, "... aborting!\n");
    exit(1);
  }
  
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = get_parse_app_args(argc, argv);

  if (strcmp(app_args.file_name, APP_ARGS_DEFAULT_FILE_NAME) == 0)
    strcpy(app_args.file_name, "twostream.ini");

  printf("Name is %s\n", app_args.file_name);
  
  rxi_ini_t *inp = rxi_ini_load(app_args.file_name);

  struct twostream_ctx ctx = create_ctx(inp); // context for init functions
  struct twostream_inp tsinp = create_twostream_inp(inp); // input parameters
  
  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = tsinp.charge,
    .mass = tsinp.mass,
    .lower = { tsinp.vel_extents[0]},
    .upper = { tsinp.vel_extents[1] }, 
    .cells = { tsinp.vel_cells },

    .evolve = 1,
    .ctx = &ctx,
    .init = evalDistFunc,

    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2" },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalFieldFunc,
  };

  // VM app
  struct gkyl_vm vm = {
    .cdim = 1, .vdim = 1,
    .lower = { -M_PI/ctx.knumber },
    .upper = { M_PI/ctx.knumber },
    .cells = { tsinp.conf_cells },
    .poly_order = 2,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 1,
    .species = { elc },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };
  // construct sim name based on input file name
  const char *inp_last_slash = strrchr(app_args.file_name, '/');
  const char *inp_no_slash = inp_last_slash ? inp_last_slash+1 : app_args.file_name;
  strncpy(vm.name, inp_no_slash, strcspn(inp_no_slash, ".ini"));
  
  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = tsinp.tend;
  double dt = tend-tcurr;

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  
  gkyl_vlasov_app_write(app, tcurr, 0);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 0);

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
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 1);
  gkyl_vlasov_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  // simulation complete, free objects
  rxi_ini_free(inp);  
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
