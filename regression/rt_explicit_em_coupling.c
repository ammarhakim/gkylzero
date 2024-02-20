#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_coldfluid.h>
#include <rt_arg_parse.h>

/**
 * This test is only intended as an example of how the explicit em coupling can 
 * work when the energies are non-relativistic. However, the test will break down
 * if the velocities become relativistic. (Such as if the field energies are too high)
 */ 

void
evalColdInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double xcell = xn[0];
  double rho;

  // create a unfform density
  rho = 20.e23;

  // no initial momentum density (fout = rho*u)
  fout[0] = rho;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
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

  // laser profile:
  const double position    = -1.0e-6;       // This point is on the laser plane
  const double E_max        = 1e10; //10.0e12;      // Maximum amplitude of the laser field (in V/m)
  const double profile_duration = 15.e-15;  // The duration of the laser (in s)
  const double profile_t_peak = 30.e-15;    // Time at which the laser reaches its peak (in s)
  const double wavelength = 0.8e-6;         // The wavelength of the laser (in m)

  // compute last edge
  double xupper = 35.0e-6;
  double xlaser = -0.0e-6;
  double xlower = -35.0e-6;
  double nx = 2800; 
  double dx = (xupper - xlower)/nx;
  double xLastEdge = xlaser+dx;

  // constants
  const double c = 299792458.0;
  const double pi = 3.14159265358979311599796346854;
  const double mu_0 = 12.56637061435917295385057353311801153679e-7;

  //Laser amplitude (t + dt/2.0)
  const double factor = 2.0*E_max/(mu_0*c*dx);
  if (xn[0]<=xLastEdge && xn[0]>(xLastEdge-dx)){
    fout[0] = 0.0;
    fout[1] = factor*sin(2.0*pi*c*t/wavelength)* exp(- pow((t - profile_t_peak),2) / pow(profile_duration,2)); // 1.01
    fout[2] = 0.0;
  } else {
    fout[0] = 0.0;
    fout[1] = 0.0;
    fout[2] = 0.0;
  }
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 2800); 

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // equation object
  struct gkyl_wv_eqn *coldf = gkyl_wv_coldfluid_new();

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = -1.602176487e-19, .mass = 9.10938215e-31,

    .equation = coldf,
    .split_type = GKYL_WAVE_FWAVE,
    .evolve = 1,
    .init = evalColdInit,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  // Fluid app
  struct gkyl_moment app_inp = {
    .name = "rt_fluid_explicit_coupling_test",

    .ndim = 1,
    .lower = { -35.0e-6 },
    .upper = { 35.0e-6 }, 
    .cells = { NX },

    .cfl_frac = 0.90, 

    .num_species = 1,
    .species = { elc }, 
    .field = {
      .epsilon0 = 8.854187817620389850536563031710750260608e-12, 
      .mu0 = 12.56637061435917295385057353311801153679e-7,
      
      .evolve = 1,
      .init = evalFieldInit,
      .app_current_func = evalAppCurrent,
      .use_explicit_em_coupling = 1,

      .bcx = { GKYL_FIELD_COPY, GKYL_FIELD_COPY },
    }
  };



  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend =  1e-13; //5.8e-13; //7.0e-12;

  int nframe = 5; // 22;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = (tend-tcurr)/nframe };

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1;
  while ((tcurr < tend) && (step <= app_args.num_steps)) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status;
    status = gkyl_moment_update(app, dt);
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
  gkyl_wv_eqn_release(coldf);
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);  
  
  return 0;
}