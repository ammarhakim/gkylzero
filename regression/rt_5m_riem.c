#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

void
evalElcInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  double gasGamma = 5./3.;
  double elcMass = 1/1836.2;
  double ionMass = 1.0;

  double rho, pr;
  if (xn[0] < 0.5) {
    rho = 1.0*elcMass/ionMass;
    pr = 5.0e-5;
  }
  else {
    rho = 0.125*elcMass/ionMass;
    pr = 5.0e-6;
  }
  fout[0] = rho;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = pr/(gasGamma-1);  
}

void
evalIonInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  double gasGamma = 5./3.;

  double rho, pr;
  if (xn[0] < 0.5) {
    rho = 1.0;
    pr = 5.0e-5;
  }
  else {
    rho = 0.125;
    pr = 5.0e-6;
  }
  fout[0] = rho;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = pr/(gasGamma-1);    
}

void
evalFieldInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  double Bx = 0.75e-2;
  double Bz = xn[0] < 0.5 ? 1.0e-2 : -1.0e-2;

  //double Bx = 0.0;
  //double Bz = 0.0;

  // electric field
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // magnetic field
  fout[3] = Bx, fout[4] = 0.0; fout[5] = Bz;

  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);
  // electron/ion equations
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(5.0/3.0);
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(5.0/3.0);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = -1.0, .mass = 1.0/1836.2,

    .equation = elc_euler,
    .evolve = 1,
    .init = evalElcInit,
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = 1.0, .mass = 1.0,

    .equation = ion_euler,
    .evolve = 1,
    .init = evalIonInit,
  };  

  // VM app
  struct gkyl_moment app_inp = {
    .name = "5m_riem",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { 1.0 }, 
    .cells = { 1024  },

    .num_species = 2,
    .species = { elc, ion },

    .field = {
      .epsilon0 = 1.0, .mu0 = 1.0,
      
      .evolve = 1,
      .init = evalFieldInit,
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 10.0;

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  gkyl_moment_app_write(app, tcurr, 0);

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

    step += 1;
  }

  gkyl_moment_app_write(app, tcurr, 1);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);
  

  // simulation complete, free resources
  gkyl_wv_eqn_release(elc_euler);
  gkyl_wv_eqn_release(ion_euler);
  gkyl_moment_app_release(app);
  
  return 0;
}
