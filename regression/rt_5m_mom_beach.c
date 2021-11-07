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
  // electron mass in kg
  double elcMass = 9.10938215e-31;
  double elcCharge = -1.602176487e-19;
  double epsilon0 = 8.854187817620389850536563031710750260608e-12;

  // plasma frequency ~ (1 - x)^5
  double wpdt = 25*(1-xn[0])*(1-xn[0])*(1-xn[0])*(1-xn[0])*(1-xn[0]);
  double xupper = 1.0;
  double xlower = 0.0;
  double dx100 = (xupper-xlower)/100;
  // dx100/speed of light
  double deltaT = dx100/299792458.0;
  double factor = deltaT*deltaT*elcCharge*elcCharge/(elcMass*epsilon0);
  double ne = wpdt*wpdt/factor;
  fout[0] = elcMass*ne;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  // 1 ev electrons (10^4 Kelvin)
  fout[4] = ne*1.602176487e-19/(gasGamma-1);  
}

void
evalFieldInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  // electric field
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // magnetic field
  fout[3] = 0.0, fout[4] = 0.0; fout[5] = 0.0;

  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalAppCurrent(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  double J0 = 1.0e-12;   // Amps/m^3.
  double xupper = 1.0;
  double xlower = 0.0;
  double nx = 400.0;
  double xLastEdge = xupper-(xupper-xlower)/nx;
  double dx100 = (xupper-xlower)/100;
  // dx100/speed of light
  double deltaT = dx100/299792458.0;
  double driveOmega = M_PI/10.0/deltaT;
  if (xn[0]>xLastEdge)
    fout[1] = -J0*sin(driveOmega*t);
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
  // electron equation
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(5.0/3.0);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = -1.602176487e-19, .mass = 9.10938215e-31,

    .equation = elc_euler,
    .evolve = 1,
    .init = evalElcInit,
  };

  // VM app
  struct gkyl_moment app_inp = {
    .name = "5m_mom_beach",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { 1.0 }, 
    .cells = { 400  },

    .num_species = 1,
    .species = { elc },

    .field = {
      .epsilon0 = 8.854187817620389850536563031710750260608e-12, 
      .mu0 = 12.56637061435917295385057353311801153679e-7,
      
      .evolve = 1,
      .init = evalFieldInit,
      .app_current_func = evalAppCurrent,
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 5.0e-9;

  int nframe = 100;
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
  

  // simulation complete, free resources
  gkyl_wv_eqn_release(elc_euler);
  gkyl_moment_app_release(app);
  
  return 0;
}
