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
  double pi = 3.141592653589793238462643383279502884;
  double gasGamma = 5.0/3.0;
  double Atwood = 0.17;
  double gHat = 0.09;
  double B0 = 0.1;
  double ionCharge = 1.0;
  double ionMass = 25.0;
  double elcMass = 1.0;
  double OmegaCi = ionCharge*B0/ionMass;
  double mu0 = 1.0;
  double n1 = 1.0;
  double n2 = n1*(1-Atwood)/(1+Atwood);
  double betaIon = 0.071;
  double Te_Ti = 0.1;
  double T_ion = (betaIon/n1)*B0*B0/(2*mu0);
  double T_elc = Te_Ti*T_ion;
  double T = (T_elc + T_ion);
  double Valf = B0/sqrt(mu0*n1*ionMass); 
  double grav = gHat*OmegaCi*Valf; // gravitational acceleration
  double maxPert = 0.01;

  double Lx = 3.0;
  double Ly = 3.75;
  double xloc = 0.5*Lx;
  double n = 0.0;
  if (x < xloc) 
    n = n1;
  else
    n = n2;
  pcg64_random_t rng = gkyl_pcg64_init(0);

  double rhoe = n*elcMass;
  double exmom = 0.0;
  double kx = 2.0*pi/Lx;
  double ky = 2.0*pi/Ly;
  for (int i=0; i<32; ++i)
    for (int j=0; j<32; ++j)
      exmom += rhoe*maxPert*Valf*gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*pi*gkyl_pcg64_rand_double(&rng))*sin(j*ky*y + 2.0*pi*gkyl_pcg64_rand_double(&rng));
  double ere = n*T_elc/(gasGamma-1) + 0.5*exmom*exmom/rhoe;

  fout[0] = rhoe;
  fout[1] = exmom; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = ere;  
}

void
evalIonInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double pi = 3.141592653589793238462643383279502884;
  double gasGamma = 5.0/3.0;
  double Atwood = 0.17;
  double gHat = 0.09;
  double B0 = 0.1;
  double ionCharge = 1.0;
  double ionMass = 25.0;
  double elcMass = 1.0;
  double OmegaCi = ionCharge*B0/ionMass;
  double mu0 = 1.0;
  double n1 = 1.0;
  double n2 = n1*(1-Atwood)/(1+Atwood);
  double betaIon = 0.071;
  double Te_Ti = 0.1;
  double T_ion = (betaIon/n1)*B0*B0/(2*mu0);
  double T_elc = Te_Ti*T_ion;
  double T = (T_elc + T_ion);
  double Valf = B0/sqrt(mu0*n1*ionMass); 
  double grav = gHat*OmegaCi*Valf; // gravitational acceleration
  double maxPert = 0.01;

  double Lx = 3.0;
  double Ly = 3.75;
  double xloc = 0.5*Lx;
  double n = 0.0;
  if (x < xloc) 
    n = n1;
  else
    n = n2;
  pcg64_random_t rng = gkyl_pcg64_init(0);

  double rhoi = n*ionMass;
  double ixmom = 0.0;
  double kx = 2.0*pi/Lx;
  double ky = 2.0*pi/Ly;
  for (int i=0; i<32; ++i)
    for (int j=0; j<32; ++j)
      ixmom += rhoi*maxPert*Valf*gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*pi*gkyl_pcg64_rand_double(&rng))*sin(j*ky*y + 2.0*pi*gkyl_pcg64_rand_double(&rng));
  double eri = n*T_ion/(gasGamma-1) + 0.5*ixmom*ixmom/rhoi;

  fout[0] = rhoi;
  fout[1] = ixmom; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = eri;    
}

void
evalFieldInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Atwood = 0.17;
  double gHat = 0.09;
  double B0 = 0.1;
  double ionCharge = 1.0;
  double ionMass = 25.0;
  double OmegaCi = ionCharge*B0/ionMass;
  double mu0 = 1.0;
  double n1 = 1.0;
  double n2 = n1*(1-Atwood)/(1+Atwood);
  double betaIon = 0.071;
  double Te_Ti = 0.1;
  double T_ion = (betaIon/n1)*B0*B0/(2*mu0);
  double T_elc = Te_Ti*T_ion;
  double T = (T_elc + T_ion);
  double Valf = B0/sqrt(mu0*n1*ionMass); 
  double grav = gHat*OmegaCi*Valf; // gravitational acceleration

  double Lx = 3.0;
  double xloc = 0.5*Lx;
  double Bz = 0.0;
  if (x < xloc) 
    Bz = sqrt(B0*B0 + 2*mu0*(ionMass*grav*n1*x));
  else
    Bz = sqrt(B0*B0 + 2*mu0*((n1- n2)*T + ionMass*grav*n1*xloc + ionMass*grav*n2*(x-xloc)));

  // electric field
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // magnetic field
  fout[3] = 0.0, fout[4] = 0.0; fout[5] = Bz;

  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalAppAccel(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double Atwood = 0.17;
  double gHat = 0.09;
  double B0 = 0.1;
  double ionCharge = 1.0;
  double ionMass = 25.0;
  double OmegaCi = ionCharge*B0/ionMass;
  double mu0 = 1.0;
  double n1 = 1.0;
  double Valf = B0/sqrt(mu0*n1*ionMass); 
  double grav = gHat*OmegaCi*Valf; // gravitational acceleration
  fout[0] = grav;
  fout[1] = 0.0;
  fout[2] = 0.0;
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
    .charge = -1.0, .mass = 1.0,

    .equation = elc_euler,
    .evolve = 1,
    .init = evalElcInit,
    .app_accel_func = evalAppAccel,

    .bcx = { GKYL_SPECIES_WALL, GKYL_SPECIES_WALL },
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = 1.0, .mass = 25.0,

    .equation = ion_euler,
    .evolve = 1,
    .init = evalIonInit,
    .app_accel_func = evalAppAccel,

    .bcx = { GKYL_SPECIES_WALL, GKYL_SPECIES_WALL },    
  };  

  // VM app
  struct gkyl_moment app_inp = {
    .name = "5m_rt",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { 3.0, 3.75 }, 
    .cells = { 64, 64 },

    .num_periodic_dir = 1,
    .periodic_dirs = { 1 },
    .cfl_frac = 1.0,

    .num_species = 2,
    .species = { elc, ion },

    .field = {
      .epsilon0 = 1.0, .mu0 = 1.0,
      .mag_error_speed_fact = 1.0,
      
      .evolve = 1,
      .init = evalFieldInit,
      
      .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(app_inp);

  // start, end and initial time-step
  // Omega_ci = 0.004, tend = 1/Omega_ci
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
