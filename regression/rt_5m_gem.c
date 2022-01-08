#include <gkyl_alloc.h>
#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

void
evalElcInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double gasGamma = 5.0/3.0;
  double elcMass = 1/25.0;
  double elcCharge = -1.0;
  double TiOverTe = 5.0;
  double lambda = 0.5;
  double n0 = 1.0;
  double nbOverN0 = 0.2;
  double B0 = 0.1;
  double plasmaBeta = 1.0;
  double TeFrac = 1.0/(1.0 + TiOverTe);
  double sech2 = (1.0/cosh(y/lambda))*(1.0/cosh(y/lambda));

  double n = n0*(sech2 + nbOverN0);
  double Jz = -(B0/lambda)*sech2;
  double Ttotal = plasmaBeta*(B0*B0)/2.0/n0;

  double rhoe = n*elcMass;
  double ezmom = (elcMass/elcCharge)*Jz*TeFrac;
  double ere = n*Ttotal*TeFrac/(gasGamma-1) + 0.5*ezmom*ezmom/rhoe;

  fout[0] = rhoe;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = ezmom;
  fout[4] = ere;  
}

void
evalIonInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double gasGamma = 5.0/3.0;
  double ionMass = 1.0;
  double ionCharge = 1.0;
  double TiOverTe = 5.0;
  double lambda = 0.5;
  double n0 = 1.0;
  double nbOverN0 = 0.2;
  double B0 = 0.1;
  double plasmaBeta = 1.0;
  double TiFrac = TiOverTe/(1.0 + TiOverTe);
  double sech2 = (1.0/cosh(y/lambda))*(1.0/cosh(y/lambda));

  double n = n0*(sech2 + nbOverN0);
  double Jz = -(B0/lambda)*sech2;
  double Ttotal = plasmaBeta*(B0*B0)/2.0/n0;

  double rhoi = n*ionMass;
  double izmom = (ionMass/ionCharge)*Jz*TiFrac;
  double eri = n*Ttotal*TiFrac/(gasGamma-1) + 0.5*izmom*izmom/rhoi;

  fout[0] = rhoi;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = izmom;
  fout[4] = eri;    
}

void
evalFieldInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double pi = 3.141592653589793238462643383279502884;
  double lambda = 0.5;
  double B0 = 0.1;
  double psi0 = 0.1*B0;

  double Lx = 25.6;
  double Ly = 12.8;

  double Bxb = B0 * tanh(y / lambda);
  double Bx = Bxb - psi0 *(pi / Ly) * cos(2 * pi * x / Lx) * sin(pi * y / Ly);
  double By = psi0 * (2 * pi / Lx) * sin(2 * pi * x / Lx) * cos(pi * y / Ly);
  double Bz = 0.0;

  // electric field
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // magnetic field
  fout[3] = Bx, fout[4] = By; fout[5] = Bz;

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
    .charge = -1.0, .mass = 1.0/25.0,

    .equation = elc_euler,
    .evolve = 1,
    .init = evalElcInit,

    .bcy = { GKYL_MOMENT_SPECIES_WALL, GKYL_MOMENT_SPECIES_WALL },
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = 1.0, .mass = 1.0,

    .equation = ion_euler,
    .evolve = 1,
    .init = evalIonInit,

    .bcy = { GKYL_MOMENT_SPECIES_WALL, GKYL_MOMENT_SPECIES_WALL },    
  };  

  // VM app
  struct gkyl_moment app_inp = {
    .name = "5m_gem",

    .ndim = 2,
    .lower = { -12.8, -6.4 },
    .upper = { 12.8, 6.4 }, 
    .cells = { 64, 32 },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },
    .cfl_frac = 1.0,

    .num_species = 2,
    .species = { elc, ion },

    .field = {
      .epsilon0 = 1.0, .mu0 = 1.0,
      .mag_error_speed_fact = 1.0,
      
      .evolve = 1,
      .init = evalFieldInit,
      
      .bcy = { GKYL_MOMENT_FIELD_COND, GKYL_MOMENT_FIELD_COND },
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 250.0;

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
  gkyl_wv_eqn_release(elc_euler);
  gkyl_wv_eqn_release(ion_euler);
  gkyl_moment_app_release(app);  
  
  return 0;
}
