#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_cold_sr_fluid.h>
#include <rt_arg_parse.h>

void
evalColdInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double xcell = xn[0];
  double rho;

  // create a density ramp
  if (xn[0]<=10.e-6) {
    rho = 1.0e10;
  } else if ((xn[0]>10.e-6) && (xn[0]<30.e-6)) {
    rho = 1.0e10 + 20.e23*((xn[0]*5.e4 + -0.5));
  } else {
    rho = 1.0e10 + 20.e23;
  }  

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
  const double position    = -1.0e-6;         // This point is on the laser plane
  const double E_max        = 10.e12;       // Maximum amplitude of the laser field (in V/m)
  const double profile_duration = 15.e-15;  // The duration of the laser (in s)
  const double profile_t_peak = 30.e-15;    // Time at which the laser reaches its peak (in s)
  const double wavelength = 0.8e-6;         // The wavelength of the laser (in m)

  // compute last edge
  double xupper = 170.0e-6;
  double xlaser = -11.0e-6;
  double xlower = -50.0e-6;
  double nx = 37548.0; // 18774.0; //10240.0;
  double dx = (xupper - xlower)/nx;
  double xLastEdge = xlaser+dx;

  // constants
  const double c = 299792458.0;
  const double pi = 3.14159265358979311599796346854;
  const double mu_0 = 12.56637061435917295385057353311801153679e-7;

  //Laser amplitude (t + dt/2.0)
  //const double ad_hoc_factor = 16.0*10*2.4e5*(1.1795);
  const double factor = 2.0*E_max/(mu_0*c*dx);
  if (xn[0]<=xLastEdge && xn[0]>(xLastEdge-dx)){
    fout[0] = 0.0;
    //fout[1] = ad_hoc_factor*sin(2.0*pi*c*t/wavelength)*(2.0/(mu_0*c))*E_max * exp(- pow((t - profile_t_peak),2) / pow(profile_duration,2));
    fout[1] = 1.00*factor*sin(2.0*pi*c*t/wavelength)* exp(- pow((t - profile_t_peak),2) / pow(profile_duration,2)); // 1.01
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 37548); // 18774

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // equation object
  struct gkyl_wv_eqn *coldf_sr = gkyl_wv_cold_sr_fluid_new();

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = -1.602176487e-19, .mass = 9.10938215e-31,

    .equation = coldf_sr,
    .split_type = GKYL_WAVE_FWAVE,
    .evolve = 1,
    .init = evalColdInit,
    //.limiter = GKYL_ZERO, // First order test

    //.bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = 1.602176487e-19, .mass = 1.67262192e-27,

    .equation = coldf_sr,
    .split_type = GKYL_WAVE_FWAVE,
    .evolve = 1,
    .init = evalColdInit,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  // Fluid app
  struct gkyl_moment app_inp = {
    .name = "cold_sr_fluid_wfa",

    .ndim = 1,
    .lower = { -50.0e-6 },
    .upper = { 170.0e-6 }, 
    .cells = { NX },

    .cfl_frac = 0.90, // 0.9,

    .num_species = 2,
    .species = { elc, ion }, //, ion 
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
  double tcurr = 0.0, tend =  5.8e-13; //7.0e-12;

  int nframe = 22;
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
  gkyl_wv_eqn_release(coldf_sr);
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);  
  
  return 0;
}

