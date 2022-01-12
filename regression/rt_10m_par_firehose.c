#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_ten_moment.h>
#include <rt_arg_parse.h>

void
evalElcInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  double pi = 3.141592653589793238462643383279502884;
  // Physical constants (using normalized code units).
  double gasGamma = 5.0/3.0;
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permiability of free space.
  double lightSpeed = 1/sqrt(mu0*epsilon0); // Speed of light.
  double ionMass = 1836.0; // Proton mass.
  double ionCharge = 1.0; // Proton charge.
  double elcMass = 1.0; // Electron mass.
  double elcCharge = -1.0; // Electron charge.
  double Te_Ti = 1.0; // Ratio of electron to proton temperature
  double vAe = 0.0125; // Electron Alfven speed (how non-relativistic the system is).
  double vAi = vAe/sqrt(ionMass/elcMass); // Proton Alfven speed.
  double n0 = 1.0; // Initial reference number density.
  double B0 = vAe*sqrt(mu0*n0*elcMass);
  double wpi = sqrt(n0*ionCharge*ionCharge/(epsilon0*ionMass));
  double di = lightSpeed/wpi;

  // Parallel firehose unstable for betaPerp - betaPar + 2 < 0.
  double beta = 300.0/pi; // Trace proton plasma beta = 2 mu0 ni Ti/B0^2 = (betaPar + 2 betaPerp)/3
  double dbeta = 100.0; // dbeta = beta_parallel - beta_perp
  double betaPar = beta + 2.*dbeta/3.; // parallel proton plasma beta = 2 mu0 ni T_paralleli/B0^2
  double betaPerp = beta - dbeta/3.; // perp proton plasma beta = 2 mu0 ni T_perpi/B0^2  

  double vtElc = vAe*sqrt(beta);
  double vtIon = vtElc/sqrt(ionMass*Te_Ti/elcMass);
  double elcTemp = vtElc*vtElc/2.0;

  double ionTempPar = vAe*vAe*(betaPar*elcMass/2);
  double ionTempPerp = vAe*vAe*(betaPerp*elcMass/2);

  double rhoe = elcMass*n0;
  double momxe = 0.0;
  double momye = 0.0;
  double momze = 0.0;
  double pxxe = n0*elcTemp + momxe*momxe/rhoe;
  double pxye = momxe*momye/rhoe;
  double pxze = momxe*momze/rhoe;
  double pyye = n0*elcTemp + momye*momye/rhoe;
  double pyze = momye*momye/rhoe;
  double pzze = n0*elcTemp + momze*momze/rhoe;

  fout[0] = rhoe;
  fout[1] = momxe; fout[2] = momye; fout[3] = momze;
  fout[4] = pxxe; fout[5] = pxye; fout[6] = pxze;  
  fout[7] = pyye; fout[8] = pyze; fout[9] = pzze;
}

void
evalIonInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  double pi = 3.141592653589793238462643383279502884;
  // Physical constants (using normalized code units).
  double gasGamma = 5.0/3.0;
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permiability of free space.
  double lightSpeed = 1/sqrt(mu0*epsilon0); // Speed of light.
  double ionMass = 1836.0; // Proton mass.
  double ionCharge = 1.0; // Proton charge.
  double elcMass = 1.0; // Electron mass.
  double elcCharge = -1.0; // Electron charge.
  double Te_Ti = 1.0; // Ratio of electron to proton temperature
  double vAe = 0.0125; // Electron Alfven speed (how non-relativistic the system is).
  double vAi = vAe/sqrt(ionMass/elcMass); // Proton Alfven speed.
  double n0 = 1.0; // Initial reference number density.
  double B0 = vAe*sqrt(mu0*n0*elcMass);
  double wpi = sqrt(n0*ionCharge*ionCharge/(epsilon0*ionMass));
  double di = lightSpeed/wpi;

  // Parallel firehose unstable for betaPerp - betaPar + 2 < 0.
  double beta = 300.0/pi; // Trace proton plasma beta = 2 mu0 ni Ti/B0^2 = (betaPar + 2 betaPerp)/3
  double dbeta = 100.0; // dbeta = beta_parallel - beta_perp
  double betaPar = beta + 2.*dbeta/3.; // parallel proton plasma beta = 2 mu0 ni T_paralleli/B0^2
  double betaPerp = beta - dbeta/3.; // perp proton plasma beta = 2 mu0 ni T_perpi/B0^2  

  double vtElc = vAe*sqrt(beta);
  double vtIon = vtElc/sqrt(ionMass*Te_Ti/elcMass);
  double elcTemp = vtElc*vtElc/2.0;

  double ionTempPar = vAe*vAe*(betaPar*elcMass/2);
  double ionTempPerp = vAe*vAe*(betaPerp*elcMass/2);

  double rhoi = ionMass*n0;
  double momxi = 0.0;
  double momyi = 0.0;
  double momzi = 0.0;
  double pxxi = n0*ionTempPar + momxi*momxi/rhoi;
  double pxyi = momxi*momyi/rhoi;
  double pxzi = momxi*momzi/rhoi;
  double pyyi = n0*ionTempPerp + momyi*momyi/rhoi;
  double pyzi = momyi*momyi/rhoi;
  double pzzi = n0*ionTempPerp + momzi*momzi/rhoi;

  fout[0] = rhoi;
  fout[1] = momxi; fout[2] = momyi; fout[3] = momzi;
  fout[4] = pxxi; fout[5] = pxyi; fout[6] = pxzi;  
  fout[7] = pyyi; fout[8] = pyzi; fout[9] = pzzi;    
}

void
evalFieldInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  double pi = 3.141592653589793238462643383279502884;
  // Physical constants (using normalized code units).
  double gasGamma = 5.0/3.0;
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permiability of free space.
  double lightSpeed = 1/sqrt(mu0*epsilon0); // Speed of light.
  double ionMass = 1836.0; // Proton mass.
  double ionCharge = 1.0; // Proton charge.
  double elcMass = 1.0; // Electron mass.
  double elcCharge = -1.0; // Electron charge.
  double Te_Ti = 1.0; // Ratio of electron to proton temperature
  double vAe = 0.0125; // Electron Alfven speed (how non-relativistic the system is).
  double vAi = vAe/sqrt(ionMass/elcMass); // Proton Alfven speed.
  double n0 = 1.0; // Initial reference number density.
  double B0 = vAe*sqrt(mu0*n0*elcMass);
  double wpi = sqrt(n0*ionCharge*ionCharge/(epsilon0*ionMass));
  double di = lightSpeed/wpi;

  pcg64_random_t rng = gkyl_pcg64_init(0);

  double Bx = B0;
  double By = 0.0, Bz = 0.0;

  double alpha = 1.0e-6*B0;
  double Lx = 12845.57;
  double kx = 2.0*pi/Lx;
  for (int i=0; i<48; ++i) {
    By -= alpha*gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*pi*gkyl_pcg64_rand_double(&rng));
    Bz -= alpha*gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*pi*gkyl_pcg64_rand_double(&rng));
  }

  // electric field
  fout[0] = 0.0; fout[1] = 0.0; fout[2] = 0.0;
  // magnetic field
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;

  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
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

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  double elc_k0 = 0.001; // rho_e ~ 10, k0e = 0.01/rho_e ~ 0.001
  double ion_k0 = 0.00002; // rho_i ~ 420, k0e = 0.01/rho_e ~ 0.00002
  
  // electron/ion equations
  struct gkyl_wv_eqn *elc_ten_moment = gkyl_wv_ten_moment_new(elc_k0);
  struct gkyl_wv_eqn *ion_ten_moment = gkyl_wv_ten_moment_new(ion_k0);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = -1.0, .mass = 1.0,
    .equation = elc_ten_moment,
    .evolve = 1,
    .init = evalElcInit,
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = 1.0, .mass = 1836.0,
    .equation = ion_ten_moment,
    .evolve = 1,
    .init = evalIonInit,
  };  

  // VM app
  struct gkyl_moment app_inp = {
    .name = "10m_par_firehose",

    .ndim = 1,
    .lower = { 0.0 },
    // Lx = 300 di ~ 12845.7
    .upper = { 12845.57 }, 
    .cells = { 560 },

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
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(app_inp);

  // start, end and initial time-step
  // Omega_ci^{-1} ~ 1e5
  double tcurr = 0.0, tend = 1.0e7;
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
