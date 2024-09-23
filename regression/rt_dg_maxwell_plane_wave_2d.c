#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  // assumes epsilon0 = mu0 = c = 1.0
  double L = 1.0;
  double kwave = 2.0;
  double lwave = 2.0;
  double pi = 3.141592653589793238462643383279502884;
  double phi = 2*pi/L*(kwave*xn[0] + lwave*xn[1]);
  double knorm = sqrt(kwave*kwave + lwave*lwave);
  double kxn = kwave/knorm;
  double kyn = lwave/knorm;
  //n = (-1, 1, 1), n_hat = 1/math.sqrt(3)
  double E0 = 1.0/sqrt(3.0);
  fout[0] = -E0*cos(phi);
  fout[1] = E0*cos(phi);
  fout[2] = E0*cos(phi);
  fout[3] = E0*cos(phi)*2*pi/L*kyn;
  fout[4] = -E0*cos(phi)*2*pi/L*kxn;
  fout[5] = E0*cos(phi)*2*pi/L*(-kxn - kyn);
  fout[6] = 0.0;
  fout[7] = 0.0;  
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
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 16);

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,
    .ctx = 0,
    .init = evalFieldFunc
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "dg_maxwell_plane_wave_2d",

    .cdim = 2, .vdim = 0,
    .lower = { 0.0, 0.0 },
    .upper = { 1.0, 1.0 },
    .cells = { NX, NY },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },

    .num_species = 0,
    .species = {  },
    .field = field,
 
    .parallelism = {
      .use_gpu = app_args.use_gpu,
    },
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 2.0;
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
      printf("** Update method failed! Aborting simulation ....\n");
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

  // simulation complete, free app
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
