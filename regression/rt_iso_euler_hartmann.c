#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_iso_euler.h>
#include <rt_arg_parse.h>

void
evalElcInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  double pi = 3.141592653589793238462643383279502884;
  double gasGamma = 5.0/3.0;
  double elcMass = 9.10938215e-31;
  double mu0 = 12.56637061435917295385057353311801153679e-7;

  double vt_elc = 3.0e6;
  double T_elc = elcMass*vt_elc*vt_elc/2.0;
  // plasma beta = 2*mu0*n*T/B^2 = 0.01
  double beta = 0.01;
  // 1 Tesla magnetic field
  double B0 = 10.0;
  // Compute density from specified parameters
  double n = beta*B0*B0/(2*mu0*T_elc);

  double rhoe = n*elcMass;
  double ere = n*T_elc/(gasGamma-1);

  fout[0] = rhoe;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
}

void
evalIonInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  double pi = 3.141592653589793238462643383279502884;
  double gasGamma = 5.0/3.0;
  double ionMass = 1.672621637e-27;
  double mu0 = 12.56637061435917295385057353311801153679e-7;
  
  double vt_ion = 6.0e4;
  double T_ion = ionMass*vt_ion*vt_ion/2.0;
  // plasma beta = 2*mu0*n*T/B^2 = 0.01
  double beta = 0.01;
  // 1 Tesla magnetic field
  double B0 = 10.0;
  // Compute density from specified parameters
  double n = beta*B0*B0/(2*mu0*T_ion);

  double rhoi = n*ionMass;
  double eri = n*T_ion/(gasGamma-1);

  fout[0] = rhoi;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
}

void
evalFieldInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  // 1 Tesla magnetic field in x direction, no other fields
  // electric field
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // magnetic field
  fout[3] = 1.0, fout[4] = 0.0; fout[5] = 0.0;

  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalAppAccel(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  // Constant force with value of acceleration on earth (9.8 m/s^2) to start
  fout[0] = 0.0;
  fout[1] = 9.8;
  fout[2] = 0.0;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);
  // electron/ion equations
  struct gkyl_wv_eqn *elc_iso_euler = gkyl_wv_iso_euler_new(3.0e6);
  struct gkyl_wv_eqn *ion_iso_euler = gkyl_wv_iso_euler_new(6.0e4);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = -1.602176487e-19, .mass = 9.10938215e-31,

    .equation = elc_iso_euler,
    .evolve = 1,
    .init = evalElcInit,
    .app_accel_func = evalAppAccel,
    .type_brag = GKYL_BRAG_MAG_FULL,
    .coll_fac = 1.0,  

    .bcx = { GKYL_SPECIES_NO_SLIP, GKYL_SPECIES_NO_SLIP },
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = 1.602176487e-19, .mass = 1.672621637e-27,

    .equation = ion_iso_euler,
    .evolve = 1,
    .init = evalIonInit,
    .app_accel_func = evalAppAccel,
    .type_brag = GKYL_BRAG_MAG_FULL,
    .coll_fac = 1.0,

    .bcx = { GKYL_SPECIES_NO_SLIP, GKYL_SPECIES_NO_SLIP },    
  };  

  // VM app
  struct gkyl_moment app_inp = {
    .name = "iso_euler_hartmann",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { 1.0 }, 
    .cells = { 64 },

    .num_periodic_dir = 0,
    .periodic_dirs = { },
    .cfl_frac = 1.0,

    .num_species = 2,
    .species = { elc, ion },
    .field = {
      .epsilon0 = 8.854187817620389850536563031710750260608e-12,
      .mu0 = 12.56637061435917295385057353311801153679e-7,
      
      .evolve = 1,
      .init = evalFieldInit,
      
      .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 1.0e-6;

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
