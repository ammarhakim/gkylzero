#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_ten_moment.h>
#include <rt_arg_parse.h>

void
evalElcInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double elcMass = 1.0;
  double elcCharge = -1.0;
  double ionMass = 36.0;
  double ionCharge = 1.0;
  double Te_Ti = 0.1;

  // Initial conditions
  double n0 = 1.0;
  double nbOverN0 = 0.001;
  double vtElc = 0.06;
  double elcTemp = vtElc*vtElc*elcMass/2.0;
  double ionTemp = elcTemp/Te_Ti;
  double vtIon = sqrt(2.0*ionTemp/ionMass);

  // electron beta
  double beta = 1.0/11.0; // total beta = 1.0 = ion beta + electron beta
  double vAe = vtElc/sqrt(beta);
  double B0 = vAe; // derived from normalization of elcMass and n0
  double vAi = vAe/sqrt(ionMass);

  double omegaCi = ionCharge*B0/ionMass;
  double omegaCe = ionCharge*B0/elcMass;

  double larmor_elc = vtElc/omegaCe;
  double larmor_ion = vtIon/omegaCi;

  // perturbation
  double noiseAmp = 0.0001;
  double mode = 8.0;

  double l = larmor_ion; // Current sheet width
  double Lx = 6.4*l;
  double Ly = 12.8*l;
  double sech2 = (1.0/cosh(y/l))*(1.0/cosh(y/l));

  double n = n0*sech2;
  double nBackground = n0*nbOverN0;
  double TeFrac = elcTemp / (elcTemp + ionTemp);
  double TiFrac = 1.0 - TeFrac;
  double i_y = 1.0;
  double i_x = mode;
  double JxNoise = -noiseAmp*(i_y*M_PI/Ly)*sin(i_y*M_PI*y/Ly)*sin(i_x*2.0*M_PI*x/Lx)/mode;
  double JyNoise = -noiseAmp*(i_x*2.0*M_PI/Lx)*cos(i_y*M_PI*y/Ly)*cos(i_x*2.0*M_PI*x/Lx)/mode;

  double Jx  = (B0/l)*(-sech2) + JxNoise;
  double Jy  = JyNoise;

  double rhoe = n*elcMass;
  double exmom = (elcMass/elcCharge)*Jx*TeFrac;
  double eymom = (elcMass/elcCharge)*Jy*TeFrac;
  double pre = n*elcTemp;

  fout[0] = rhoe;
  fout[1] = exmom; fout[2] = eymom; fout[3] = 0.0;
  fout[4] = pre + exmom*exmom/rhoe; fout[5] = exmom*eymom/rhoe; fout[6] = 0.0;  
  fout[7] = pre + eymom*eymom/rhoe; fout[8] = 0.0; fout[9] = pre;  
}

void
evalIonInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double elcMass = 1.0;
  double elcCharge = -1.0;
  double ionMass = 36.0;
  double ionCharge = 1.0;
  double Te_Ti = 0.1;

  // Initial conditions
  double n0 = 1.0;
  double nbOverN0 = 0.001;
  double vtElc = 0.06;
  double elcTemp = vtElc*vtElc*elcMass/2.0;
  double ionTemp = elcTemp/Te_Ti;
  double vtIon = sqrt(2.0*ionTemp/ionMass);

  // electron beta
  double beta = 1.0/11.0; // total beta = 1.0 = ion beta + electron beta
  double vAe = vtElc/sqrt(beta);
  double B0 = vAe; // derived from normalization of elcMass and n0
  double vAi = vAe/sqrt(ionMass);

  double omegaCi = ionCharge*B0/ionMass;
  double omegaCe = ionCharge*B0/elcMass;

  double larmor_elc = vtElc/omegaCe;
  double larmor_ion = vtIon/omegaCi;

  // perturbation
  double noiseAmp = 0.0001;
  double mode = 8.0;

  double l = larmor_ion; // Current sheet width
  double Lx = 6.4*l;
  double Ly = 12.8*l;
  double sech2 = (1.0/cosh(y/l))*(1.0/cosh(y/l));

  double n = n0*sech2;
  double nBackground = n0*nbOverN0;
  double TeFrac = elcTemp / (elcTemp + ionTemp);
  double TiFrac = 1.0 - TeFrac;
  double i_y = 1.0;
  double i_x = mode;
  double JxNoise = -noiseAmp*(i_y*M_PI/Ly)*sin(i_y*M_PI*y/Ly)*sin(i_x*2.0*M_PI*x/Lx)/mode;
  double JyNoise = -noiseAmp*(i_x*2.0*M_PI/Lx)*cos(i_y*M_PI*y/Ly)*cos(i_x*2.0*M_PI*x/Lx)/mode;

  double Jx  = (B0/l)*(-sech2) + JxNoise;
  double Jy  = JyNoise;

  double rhoi = n*ionMass;
  double ixmom = (ionMass/ionCharge)*Jx*TiFrac;
  double iymom = (ionMass/ionCharge)*Jy*TiFrac;
  double pri = n*ionTemp;

  fout[0] = rhoi;
  fout[1] = ixmom; fout[2] = iymom; fout[3] = 0.0;
  fout[4] = pri + ixmom*ixmom/rhoi; fout[5] = ixmom*iymom/rhoi; fout[6] = 0.0;  
  fout[7] = pri + iymom*iymom/rhoi; fout[8] = 0.0; fout[9] = pri; 
}

void
evalFieldInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double elcMass = 1.0;
  double elcCharge = -1.0;
  double ionMass = 36.0;
  double ionCharge = 1.0;
  double Te_Ti = 0.1;

  // Initial conditions
  double n0 = 1.0;
  double nbOverN0 = 0.001;
  double vtElc = 0.06;
  double elcTemp = vtElc*vtElc*elcMass/2.0;
  double ionTemp = elcTemp/Te_Ti;
  double vtIon = sqrt(2.0*ionTemp/ionMass);

  // electron beta
  double beta = 1.0/11.0; // total beta = 1.0 = ion beta + electron beta
  double vAe = vtElc/sqrt(beta);
  double B0 = vAe; // derived from normalization of elcMass and n0
  double vAi = vAe/sqrt(ionMass);

  double omegaCi = ionCharge*B0/ionMass;
  double omegaCe = ionCharge*B0/elcMass;

  double larmor_elc = vtElc/omegaCe;
  double larmor_ion = vtIon/omegaCi;

  // perturbation
  double noiseAmp = 0.0001;
  double mode = 8.0;

  double l = larmor_ion; // Current sheet width
  double Lx = 6.4*l;
  double Ly = 12.8*l;

  double i_y = 1.0;
  double i_x = mode;
  double BzNoise = noiseAmp*cos(i_y*M_PI*y/Ly)*sin(i_x*2.0*M_PI*x/Lx)/mode;
  double Bzb = -B0*tanh(y/l);
  double Bx = 0.0;
  double By = 0.0;
  double Bz = Bzb + BzNoise;

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
  struct gkyl_wv_eqn *elc_ten_moment = gkyl_wv_ten_moment_new();
  struct gkyl_wv_eqn *ion_ten_moment = gkyl_wv_ten_moment_new();

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = -1.0, .mass = 1.0,
    // d_e = 1.0, k0 = 1/d_e
    .k0 = 1.0,
    .equation = elc_ten_moment,
    .evolve = 1,
    .init = evalElcInit,

    .bcy = { GKYL_MOMENT_SPECIES_WALL, GKYL_MOMENT_SPECIES_WALL },
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = 1.0, .mass = 36.0,
    // d_i = 6.0, k0 = 1/d_i
    .k0 = 1.0/6.0,
    .equation = ion_ten_moment,
    .evolve = 1,
    .init = evalIonInit,

    .bcy = { GKYL_MOMENT_SPECIES_WALL, GKYL_MOMENT_SPECIES_WALL },    
  };  

  // VM app
  struct gkyl_moment app_inp = {
    .name = "10m_lhdi",

    .ndim = 2,
    .lower = { -3.2*5.720775535473554, -6.4*5.720775535473554 },
    .upper = { 3.2*5.720775535473554, 6.4*5.720775535473554 }, 
    .cells = { 64, 128  },

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
  // OmegaCi^{-1} ~ 181 -> tend ~ 6 OmegaCi^{-1}
  double tcurr = 0.0, tend = 1100.0;

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
  gkyl_wv_eqn_release(elc_ten_moment);
  gkyl_wv_eqn_release(ion_ten_moment);
  gkyl_moment_app_release(app);  
  
  return 0;
}
