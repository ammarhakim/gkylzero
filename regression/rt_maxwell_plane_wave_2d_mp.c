#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 128);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 128);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  
  struct gkyl_moment app_inp = {
    .name = "maxwell_plane_wave_2d_mp",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { 1.0, 1.0 }, 
    .cells = { NX, NY },

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },
    .cfl_frac = 1.0,

    .scheme_type = GKYL_MOMENT_MP,
    .mp_recon = app_args.mp_recon,
    .skip_mp_limiter = app_args.skip_limiters,

    .field = {
      .epsilon0 = 1.0, .mu0 = 1.0,
      .evolve = 1,
      .limiter = GKYL_NO_LIMITER,
      .init = evalFieldInit,
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 2.0;

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  gkyl_moment_app_write(app, tcurr, 0);
  gkyl_moment_app_calc_field_energy(app, 0.0);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);

    gkyl_moment_app_calc_field_energy(app, tcurr);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    step += 1;
  }

  gkyl_moment_app_write(app, tcurr, 1);
  gkyl_moment_app_write_field_energy(app);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  // simulation complete, free resources
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);
  
  return 0;
}
