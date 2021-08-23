#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_ten_moment.h>
#include <app_arg_parse.h>

void
evalElcInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double pi = 3.141592653589793238462643383279502884;
  // Physical constants (using normalized code units).
  double gasGamma = 5.0/3.0;
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permiability of free space.
  double lightSpeed = 1/sqrt(mu0*epsilon0); // Speed of light.
  double ionMass = 1.0; // Proton mass.
  double ionCharge = 1.0; // Proton charge.
  double elcMass = ionMass/100.0; // Electron mass.
  double elcCharge = -1.0; // Electron charge.
  double vAe = 0.2;
  double n0 = 1.0;
  double B0 = vAe*sqrt(n0*elcMass);
  double wpi = sqrt(n0*ionCharge*ionCharge/(epsilon0*ionMass));
  double di = lightSpeed/wpi;
  double w0 = 1.0*di;
  double psi0 = 0.1*B0*di;
  double guide1 = 0.099*B0;
  double guide2 = 0.099*B0;

  double b1 = 1.696*B0;
  double b2 = 1.00*B0;
  double n1 = 0.062*n0;
  double n2 = 1.0*n0;
  double beta2 = 2.748;
  double T_i2 = beta2*(b2*b2) / (2.0*n2*mu0);
  double T_i1_T_i2 = 7.73 / 1.374;
  double T_e1_T_i2 = 1.288 / 1.374;

  double T_e1 = T_i2*T_e1_T_i2;
  double T_i1  = T_i2*T_i1_T_i2;

  // derived quantities (asymptotic values)
  double T_e2 = (0.5*(b1*b1-b2*b2)+0.5*(guide1*guide1-guide2*guide2)+n1*(T_i1+T_e1)-n2*T_i2)/n2; //so the system is in force balance

  double Lx = 40.96*di;
  double Ly = 20.48*di;

  double b1x = 0.5*(b2+b1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(b2-b1);
  double b1y = 0.0;
  double b1z = 0.5*(guide2-guide1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(guide2+guide1);

  double Tvali = 0.5*(T_i2-T_i1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(T_i2+T_i1);
  double Tvale = 0.5*(T_e2-T_e1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(T_e2+T_e1);
  double n = (0.5*(b1*b1 - b1x*b1x) +0.5*(guide1*guide1 - b1z*b1z) + n1*(T_i1+T_e1))/(Tvali+Tvale);

  double Bx = b1x - psi0*4.0*pi/Ly*sin(2.0*pi*x/Lx)*sin(4.0*pi*y/Ly);
  double By = b1y + psi0*2.0*pi/Lx*cos(2.0*pi*x/Lx)*(1-cos(4.0*pi*y/Ly));
  double Bz = b1z;

  double TeFrac = Tvale/(Tvale + Tvali);
  double TiFrac = Tvali/(Tvale + Tvali);

  double Jx = 0.5*(guide2-guide1)/w0*((1.0/cosh((y-Ly*.25)/w0))*(1.0/cosh((y-Ly*.25)/w0)) 
    - (1.0/cosh((y-Ly*.75)/w0))*(1.0/cosh((y-Ly*.75)/w0)) 
    + (1.0/cosh((y-Ly*1.25)/w0))*(1.0/cosh((y-Ly*1.25)/w0)) 
    - (1.0/cosh((y+Ly*.25)/w0))*(1.0/cosh((y+Ly*.25)/w0)));
  double Jy = 0.0;
  double Jz  = -0.5*(b2+b1)/w0*((1.0/cosh((y-Ly*.25)/w0))*(1.0/cosh((y-Ly*.25)/w0)) 
    - (1.0/cosh((y-Ly*.75)/w0))*(1.0/cosh((y-Ly*.75)/w0)) 
    + (1.0/cosh((y-Ly*1.25)/w0))*(1.0/cosh((y-Ly*1.25)/w0)) 
    - (1.0/cosh((y+Ly*.25)/w0))*(1.0/cosh((y+Ly*.25)/w0))) 
    - psi0*sin(2.0*pi*x/Lx)*((2.0*pi/Lx)*(2.0*pi/Lx)*(1 - cos(4.0*pi*y/Ly)) + (4.0*pi/Ly)*(4.0*pi/Ly)*cos(4.0*pi*y/Ly));
   
  double Jxe = Jx*TeFrac;
  double Jye = Jy*TeFrac;
  double Jze = Jz*TeFrac;

  double rhoe = elcMass*n;
  double momxe = (elcMass/elcCharge)*Jxe;
  double momye = (elcMass/elcCharge)*Jye;
  double momze = (elcMass/elcCharge)*Jze;
  double pxxe = n*Tvale + momxe*momxe/rhoe;
  double pxye = momxe*momye/rhoe;
  double pxze = momxe*momze/rhoe;
  double pyye = n*Tvale + momye*momye/rhoe;
  double pyze = momye*momye/rhoe;
  double pzze = n*Tvale + momze*momze/rhoe;

  fout[0] = rhoe;
  fout[1] = momxe; fout[2] = momye; fout[3] = momze;
  fout[4] = pxxe; fout[5] = pxye; fout[6] = pxze;  
  fout[7] = pyye; fout[8] = pyze; fout[9] = pzze;
}

void
evalIonInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double pi = 3.141592653589793238462643383279502884;
  // Physical constants (using normalized code units).
  double gasGamma = 5.0/3.0;
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permiability of free space.
  double lightSpeed = 1/sqrt(mu0*epsilon0); // Speed of light.
  double ionMass = 1.0; // Proton mass.
  double ionCharge = 1.0; // Proton charge.
  double elcMass = ionMass/100.0; // Electron mass.
  double elcCharge = -1.0; // Electron charge.
  double vAe = 0.2;
  double n0 = 1.0;
  double B0 = vAe*sqrt(n0*elcMass);
  double wpi = sqrt(n0*ionCharge*ionCharge/(epsilon0*ionMass));
  double di = lightSpeed/wpi;
  double w0 = 1.0*di;
  double psi0 = 0.1*B0*di;
  double guide1 = 0.099*B0;
  double guide2 = 0.099*B0;

  double b1 = 1.696*B0;
  double b2 = 1.00*B0;
  double n1 = 0.062*n0;
  double n2 = 1.0*n0;
  double beta2 = 2.748;
  double T_i2 = beta2*(b2*b2) / (2.0*n2*mu0);
  double T_i1_T_i2 = 7.73 / 1.374;
  double T_e1_T_i2 = 1.288 / 1.374;

  double T_e1 = T_i2*T_e1_T_i2;
  double T_i1  = T_i2*T_i1_T_i2;

  // derived quantities (asymptotic values)
  double T_e2 = (0.5*(b1*b1-b2*b2)+0.5*(guide1*guide1-guide2*guide2)+n1*(T_i1+T_e1)-n2*T_i2)/n2; //so the system is in force balance

  double Lx = 40.96*di;
  double Ly = 20.48*di;

  double b1x = 0.5*(b2+b1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(b2-b1);
  double b1y = 0.0;
  double b1z = 0.5*(guide2-guide1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(guide2+guide1);
  
  double Tvali = 0.5*(T_i2-T_i1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(T_i2+T_i1);
  double Tvale = 0.5*(T_e2-T_e1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(T_e2+T_e1);
  double n = (0.5*(b1*b1 - b1x*b1x) +0.5*(guide1*guide1 - b1z*b1z) + n1*(T_i1+T_e1))/(Tvali+Tvale);

  double Bx = b1x - psi0*4.0*pi/Ly*sin(2.0*pi*x/Lx)*sin(4.0*pi*y/Ly);
  double By = b1y + psi0*2.0*pi/Lx*cos(2.0*pi*x/Lx)*(1-cos(4.0*pi*y/Ly));
  double Bz = b1z;

  double TeFrac = Tvale/(Tvale + Tvali);
  double TiFrac = Tvali/(Tvale + Tvali);

  double Jx = 0.5*(guide2-guide1)/w0*((1.0/cosh((y-Ly*.25)/w0))*(1.0/cosh((y-Ly*.25)/w0)) 
    - (1.0/cosh((y-Ly*.75)/w0))*(1.0/cosh((y-Ly*.75)/w0)) 
    + (1.0/cosh((y-Ly*1.25)/w0))*(1.0/cosh((y-Ly*1.25)/w0)) 
    - (1.0/cosh((y+Ly*.25)/w0))*(1.0/cosh((y+Ly*.25)/w0)));
  double Jy = 0.0;
  double Jz  = -0.5*(b2+b1)/w0*((1.0/cosh((y-Ly*.25)/w0))*(1.0/cosh((y-Ly*.25)/w0)) 
    - (1.0/cosh((y-Ly*.75)/w0))*(1.0/cosh((y-Ly*.75)/w0)) 
    + (1.0/cosh((y-Ly*1.25)/w0))*(1.0/cosh((y-Ly*1.25)/w0)) 
    - (1.0/cosh((y+Ly*.25)/w0))*(1.0/cosh((y+Ly*.25)/w0))) 
    - psi0*sin(2.0*pi*x/Lx)*((2.0*pi/Lx)*(2.0*pi/Lx)*(1 - cos(4.0*pi*y/Ly)) + (4.0*pi/Ly)*(4.0*pi/Ly)*cos(4.0*pi*y/Ly));
   
  double Jxi = Jx*TiFrac;
  double Jyi = Jy*TiFrac;
  double Jzi = Jz*TiFrac;

  double rhoi = ionMass*n;
  double momxi = (ionMass/ionCharge)*Jxi;
  double momyi = (ionMass/ionCharge)*Jyi;
  double momzi = (ionMass/ionCharge)*Jzi;
  double pxxi = n*Tvali + momxi*momxi/rhoi;
  double pxyi = momxi*momyi/rhoi;
  double pxzi = momxi*momzi/rhoi;
  double pyyi = n*Tvali + momyi*momyi/rhoi;
  double pyzi = momyi*momyi/rhoi;
  double pzzi = n*Tvali + momzi*momzi/rhoi;

  fout[0] = rhoi;
  fout[1] = momxi; fout[2] = momyi; fout[3] = momzi;
  fout[4] = pxxi; fout[5] = pxyi; fout[6] = pxzi;  
  fout[7] = pyyi; fout[8] = pyzi; fout[9] = pzzi;   
}

void
evalFieldInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double pi = 3.141592653589793238462643383279502884;
  // Physical constants (using normalized code units).
  double gasGamma = 5.0/3.0;
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permiability of free space.
  double lightSpeed = 1/sqrt(mu0*epsilon0); // Speed of light.
  double ionMass = 1.0; // Proton mass.
  double ionCharge = 1.0; // Proton charge.
  double elcMass = ionMass/100.0; // Electron mass.
  double elcCharge = -1.0; // Electron charge.
  double vAe = 0.2;
  double n0 = 1.0;
  double B0 = vAe*sqrt(n0*elcMass);
  double wpi = sqrt(n0*ionCharge*ionCharge/(epsilon0*ionMass));
  double di = lightSpeed/wpi;
  double w0 = 1.0*di;
  double psi0 = 0.1*B0*di;
  double guide1 = 0.099*B0;
  double guide2 = 0.099*B0;

  double b1 = 1.696*B0;
  double b2 = 1.00*B0;
  double n1 = 0.062*n0;
  double n2 = 1.0*n0;
  double beta2 = 2.748;
  double T_i2 = beta2*(b2*b2) / (2.0*n2*mu0);
  double T_i1_T_i2 = 7.73 / 1.374;
  double T_e1_T_i2 = 1.288 / 1.374;

  double T_e1 = T_i2*T_e1_T_i2;
  double T_i1  = T_i2*T_i1_T_i2;

  // derived quantities (asymptotic values)
  double T_e2 = (0.5*(b1*b1-b2*b2)+0.5*(guide1*guide1-guide2*guide2)+n1*(T_i1+T_e1)-n2*T_i2)/n2; //so the system is in force balance

  double Lx = 40.96*di;
  double Ly = 20.48*di;

  double b1x = 0.5*(b2+b1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(b2-b1);
  double b1y = 0.0;
  double b1z = 0.5*(guide2-guide1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(guide2+guide1);
  
  double Tvali = 0.5*(T_i2-T_i1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(T_i2+T_i1);
  double Tvale = 0.5*(T_e2-T_e1)*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1)+0.5*(T_e2+T_e1);
  double n = (0.5*(b1*b1 - b1x*b1x) +0.5*(guide1*guide1 - b1z*b1z) + n1*(T_i1+T_e1))/(Tvali+Tvale);

  double Bx = b1x - psi0*4.0*pi/Ly*sin(2.0*pi*x/Lx)*sin(4.0*pi*y/Ly);
  double By = b1y + psi0*2.0*pi/Lx*cos(2.0*pi*x/Lx)*(1-cos(4.0*pi*y/Ly));
  double Bz = b1z;

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
    .charge = -1.0, .mass = 1.0/100.0,
    .k0 = 1.0,
    .equation = elc_ten_moment,
    .evolve = 1,
    .init = evalElcInit,
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = 1.0, .mass = 1.0,
    .k0 = 0.1,
    .equation = ion_ten_moment,
    .evolve = 1,
    .init = evalIonInit, 
  };  

  // VM app
  struct gkyl_moment app_inp = {
    .name = "10m_burch",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { 40.96, 20.48 }, 
    .cells = { 256, 128  },

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },
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
  gkyl_wv_eqn_release(elc_ten_moment);
  gkyl_wv_eqn_release(ion_ten_moment);
  gkyl_moment_app_release(app);  
  
  return 0;
}
