#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_iso_euler.h>
#include <rt_arg_parse.h>

void
evalElcInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double elcMass = 1.0/1836.2;

  fout[0] = elcMass*(1 - cos(x));
  fout[1] = elcMass*(1 - cos(x))*cos(x); fout[2] = 0.0; fout[3] = 0.0; 
}

void
evalIonInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double ionMass = 1.0;

  fout[0] = ionMass*(1 - sin(x));
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0; 
}

void
evalFieldInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  // electric field
  fout[0] = sin(x), fout[1] = 0.0; fout[2] = 0.0;
  // magnetic field
  fout[3] = 0.0, fout[4] = 0.0; fout[5] = 0.0;

  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalElcAppAccel(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  double elcCharge = -1.0;
  double elcMass = 1.0/1836.2;
  double qbym = elcCharge/elcMass;
  // first two terms are momentum sources and need to divide by qbym
  fout[0] = 4.0/3.0*cos(t)*cos(t)*cos(x)*sin(x)/qbym - sin(t)*sin(x)/qbym + (cos(t)*sin(x));
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalIonAppAccel(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  double ionCharge = 1.0;
  double ionMass = 1.0;
  double qbym = ionCharge/ionMass;
  // first two terms are momentum sources and need to divide by qbym
  fout[0] = -4.0/3.0*sin(t)*sin(t)*cos(x)*sin(x)/qbym + cos(t)*cos(x)/qbym - (cos(t)*sin(x));
  fout[1] = 0.0;
  fout[2] = 0.0;
}


void
evalAppCurrent(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = -cos(x)*sin(t) + cos(t)*sin(x) + sin(t)*sin(x);
  fout[1] = 0.0;
  fout[2] = 0.0;
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
  // electron/ion equations
  // equal thermal velocities (electrons hotter than ions)
  // thermal velocities sub-relativistic (c = 1)
  struct gkyl_wv_eqn *elc_iso_euler = gkyl_wv_iso_euler_new(0.01);
  struct gkyl_wv_eqn *ion_iso_euler = gkyl_wv_iso_euler_new(0.01);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = -1.0, .mass = 1.0/1836.2,

    .equation = elc_iso_euler,
    .evolve = 1,
    .init = evalElcInit,
    .app_accel_func = evalElcAppAccel,
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = 1.0, .mass = 1.0,

    .equation = ion_iso_euler,
    .evolve = 1,
    .init = evalIonInit,
    .app_accel_func = evalIonAppAccel,
  };  

  // VM app
  struct gkyl_moment app_inp = {
    .name = "5m_iso_moms",

    .ndim = 1,
    .lower = { -2.0*M_PI },
    .upper = { 2.0*M_PI }, 
    .cells = { 128 },

    .num_species = 2,
    .species = { elc, ion },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .field = {
      .epsilon0 = 1.0, .mu0 = 1.0,
      
      .evolve = 1,
      .init = evalFieldInit,
      .app_current_func = evalAppCurrent,
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 1.0;
  int nframe = 1;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);
  long step = 1;
  while ((tcurr < tend) && (step <= app_args.num_steps)) {
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

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Source updates took %g secs\n", stat.sources_tm);
  printf("Total updates took %g secs\n", stat.total_tm);  

  // simulation complete, free resources
  gkyl_wv_eqn_release(elc_iso_euler);
  gkyl_wv_eqn_release(ion_iso_euler);
  gkyl_moment_app_release(app);
  
  return 0;
}
