#include <math.h>
#include <stdio.h>

#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>

void
evalFieldInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
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
  // VM app
  struct gkyl_moment app_inp = {
    .name = "maxwell_plane_wave",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { 1.0, 1.0 }, 
    .cells = { 128, 128 },

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },
    .cfl_frac = 1.0,

    .field = {
      .epsilon0 = 1.0, .mu0 = 1.0,
      .elcErrorSpeedFactor = 0.0, .mgnErrorSpeedFactor = 0.0,
      .evolve = 1,
      .ctx = NULL,
      .limiter = GKYL_NO_LIMITER,
      .init = evalFieldInit,
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 2.0;

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  gkyl_moment_app_write(app, tcurr, 0);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1;
  while (tcurr < tend) {
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
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);
  
  return 0;
}
