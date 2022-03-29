#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

struct burch_ctx {
  // parameters for plasma
  double gasGamma;
  double charge;
  double elcMass, ionMass;
  double T_elc, T_ion;
  double beta;
  double B0;
  double n0;

  double epsilon0, mu0;

  double coll_fac;
  double tau;
  double eta_par, eta_perp;

  double w0, Lx, Ly;
};

struct burch_ctx
create_ctx(void)
{
  double gasGamma = 5.0/3.0;
  double charge = 1.0;
  double ionMass = 1.0;
  double elcMass = ionMass/100.0;
  double epsilon0 = 1.0;
  double mu0 = 1.0;
  // double charge = 1.602176487e-19;
  // double elcMass = 9.10938215e-31;
  // double ionMass = 1.672621637e-27;
  // double epsilon0 = 8.854187817620389850536563031710750260608e-12;
  // double mu0 = 12.56637061435917295385057353311801153679e-7;
  double c = 1.0/sqrt(epsilon0*mu0);
  double vAe = c;
  // plasma density in particles per m^3
  double n0 = 1.0;
  // Compute magnetic field from specified parameters
  double B0 = vAe*sqrt(mu0*n0*elcMass);

  // plasma beta = 2*mu0*n*T/B^2
  double beta = 2.748;
  double b1 = 1.696*B0;
  double b2 = 1.00*B0;
  double n1 = 0.062*n0;
  double n2 = 1.0*n0;
  double T_ion = beta*(b2*b2)/(2.0*n2*mu0);
  double T_e1_T_i2 = 1.288/1.374;
  double T_elc = T_ion*T_e1_T_i2;
  
  double coll_fac = 1e3;
  double tau = coll_fac*6.0*sqrt(2.0*M_PI*elcMass*T_elc*M_PI*T_elc*M_PI*T_elc)*epsilon0*epsilon0/(charge*charge*charge*charge*n0);

  double eta_par = 0.96*n0*T_ion*sqrt(2.0)*tau;
  double eta_perp = 0.3*n0*T_ion/(tau*charge*charge*B0*B0/(ionMass*ionMass));

  double wpi = sqrt(charge*charge*n0/(epsilon0*ionMass));
  double di = c/wpi;
  double w0 = 1.0*di;
  double Lx = 40.96*di;
  double Ly = 20.48*di;
  struct burch_ctx ctx = {
    .gasGamma = gasGamma,
    .charge = charge,
    .elcMass = elcMass,
    .ionMass = ionMass,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .T_elc = T_elc,
    .T_ion = T_ion,
    .beta = beta,
    .B0 = B0,
    .n0 = n0,
    .coll_fac = coll_fac,
    .tau = tau,
    .eta_par = eta_par,
    .eta_perp = eta_perp,
    .w0 = w0,
    .Lx = Lx,
    .Ly = Ly,
  };
  return ctx;
}

void
evalElcInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  struct burch_ctx *app = ctx;
  double pi = 3.141592653589793238462643383279502884;
  // Physical constants (using normalized code units).
  double gasGamma = app->gasGamma;
  double epsilon0 = app->epsilon0; // Permittivity of free space.
  double mu0 = app->mu0; // Permiability of free space.
  double lightSpeed = 1/sqrt(mu0*epsilon0); // Speed of light.
  double ionMass = app->ionMass; // Proton mass.
  double ionCharge = app->charge; // Proton charge.
  double elcMass = app->elcMass; // Electron mass.
  double elcCharge = -1.0*app->charge; // Electron charge.
  double n0 = app->n0;
  double B0 = app->B0;
  double vAe = B0/sqrt(mu0*n0*elcMass);

  double wpi = sqrt(ionCharge*ionCharge*n0/(epsilon0*ionMass));
  double di = lightSpeed/wpi;
  double w0 = 1.0*di;
  double Lx = 40.96*di;
  double Ly = 20.48*di;

  //double w0 = app->w0;
  double psi0 = 0.1*B0*w0;
  double guide1 = 0.099*B0;
  double guide2 = 0.099*B0;

  double b1 = 1.696*B0;
  double b2 = 1.00*B0;
  double n1 = 0.062*n0;
  double n2 = 1.0*n0;
  double beta2 = app->beta;
  double T_i2 = beta2*(b2*b2) / (2.0*n2*mu0);
  double T_i1_T_i2 = 7.73 / 1.374;
  double T_e1_T_i2 = 1.288 / 1.374;

  double T_e1 = T_i2*T_e1_T_i2;
  double T_i1  = T_i2*T_i1_T_i2;

  // derived quantities (asymptotic values)
  double T_e2 = (0.5*(b1*b1-b2*b2)+0.5*(guide1*guide1-guide2*guide2)+n1*(T_i1+T_e1)-n2*T_i2)/n2; //so the system is in force balance

  //double Lx = app->Lx;
  //double Ly = app->Ly;

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

  double Jx = 0.5/mu0*(guide2-guide1)/w0*((1.0/cosh((y-Ly*.25)/w0))*(1.0/cosh((y-Ly*.25)/w0)) 
    - (1.0/cosh((y-Ly*.75)/w0))*(1.0/cosh((y-Ly*.75)/w0)) 
    + (1.0/cosh((y-Ly*1.25)/w0))*(1.0/cosh((y-Ly*1.25)/w0)) 
    - (1.0/cosh((y+Ly*.25)/w0))*(1.0/cosh((y+Ly*.25)/w0)));
  double Jy = 0.0;
  double Jz = -0.5/mu0*(b2+b1)/w0*((1.0/cosh((y-Ly*.25)/w0))*(1.0/cosh((y-Ly*.25)/w0))
    - (1.0/cosh((y-Ly*.75)/w0))*(1.0/cosh((y-Ly*.75)/w0)) 
    + (1.0/cosh((y-Ly*1.25)/w0))*(1.0/cosh((y-Ly*1.25)/w0)) 
    - (1.0/cosh((y+Ly*.25)/w0))*(1.0/cosh((y+Ly*.25)/w0)))
    - 1.0/mu0*psi0*sin(2.0*pi*x/Lx)*((2.0*pi/Lx)*(2.0*pi/Lx)*(1 - cos(4.0*pi*y/Ly)) + (4.0*pi/Ly)*(4.0*pi/Ly)*cos(4.0*pi*y/Ly));
   
  double Jxe = Jx*TeFrac;
  double Jye = Jy*TeFrac;
  double Jze = Jz*TeFrac;

  double rhoe = elcMass*n;
  double momxe = (elcMass/elcCharge)*Jxe;
  double momye = (elcMass/elcCharge)*Jye;
  double momze = (elcMass/elcCharge)*Jze;
  double ere = n*Tvale/(gasGamma-1) + 0.5*(momxe*momxe + momye*momye + momze*momze)/rhoe;

  fout[0] = rhoe;
  fout[1] = momxe; fout[2] = momye; fout[3] = momze;
  fout[4] = ere;  
}

void
evalIonInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  struct burch_ctx *app = ctx;
  double pi = 3.141592653589793238462643383279502884;
  // Physical constants (using normalized code units).
  double gasGamma = app->gasGamma;
  double epsilon0 = app->epsilon0; // Permittivity of free space.
  double mu0 = app->mu0; // Permiability of free space.
  double lightSpeed = 1.0/sqrt(mu0*epsilon0); // Speed of light.
  double ionMass = app->ionMass; // Proton mass.
  double ionCharge = app->charge; // Proton charge.
  double elcMass = app->elcMass; // Electron mass.
  double elcCharge = -1.0*app->charge; // Electron charge.
  double n0 = app->n0;
  double B0 = app->B0;
  double vAe = B0/sqrt(mu0*n0*elcMass);

  double w0 = app->w0;
  double psi0 = 0.1*B0*w0;
  double guide1 = 0.099*B0;
  double guide2 = 0.099*B0;

  double b1 = 1.696*B0;
  double b2 = 1.00*B0;
  double n1 = 0.062*n0;
  double n2 = 1.0*n0;
  double beta2 = app->beta;
  double T_i2 = beta2*(b2*b2) / (2.0*n2*mu0);
  double T_i1_T_i2 = 7.73 / 1.374;
  double T_e1_T_i2 = 1.288 / 1.374;

  double T_e1 = T_i2*T_e1_T_i2;
  double T_i1  = T_i2*T_i1_T_i2;

  // derived quantities (asymptotic values)
  double T_e2 = (0.5*(b1*b1-b2*b2)+0.5*(guide1*guide1-guide2*guide2)+n1*(T_i1+T_e1)-n2*T_i2)/n2; //so the system is in force balance

  double Lx = app->Lx;
  double Ly = app->Ly;

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

  double Jx = 0.5/mu0*(guide2-guide1)/w0*((1.0/cosh((y-Ly*.25)/w0))*(1.0/cosh((y-Ly*.25)/w0)) 
    - (1.0/cosh((y-Ly*.75)/w0))*(1.0/cosh((y-Ly*.75)/w0)) 
    + (1.0/cosh((y-Ly*1.25)/w0))*(1.0/cosh((y-Ly*1.25)/w0)) 
    - (1.0/cosh((y+Ly*.25)/w0))*(1.0/cosh((y+Ly*.25)/w0)));
  double Jy = 0.0;
  double Jz  = -0.5/mu0*(b2+b1)/w0*((1.0/cosh((y-Ly*.25)/w0))*(1.0/cosh((y-Ly*.25)/w0)) 
    - (1.0/cosh((y-Ly*.75)/w0))*(1.0/cosh((y-Ly*.75)/w0)) 
    + (1.0/cosh((y-Ly*1.25)/w0))*(1.0/cosh((y-Ly*1.25)/w0)) 
    - (1.0/cosh((y+Ly*.25)/w0))*(1.0/cosh((y+Ly*.25)/w0))) 
    - 1.0/mu0*psi0*sin(2.0*pi*x/Lx)*((2.0*pi/Lx)*(2.0*pi/Lx)*(1 - cos(4.0*pi*y/Ly)) + (4.0*pi/Ly)*(4.0*pi/Ly)*cos(4.0*pi*y/Ly));
   
  double Jxi = Jx*TiFrac;
  double Jyi = Jy*TiFrac;
  double Jzi = Jz*TiFrac;

  double rhoi = ionMass*n;
  double momxi = (ionMass/ionCharge)*Jxi;
  double momyi = (ionMass/ionCharge)*Jyi;
  double momzi = (ionMass/ionCharge)*Jzi;
  double eri = n*Tvali/(gasGamma-1) + 0.5*(momxi*momxi + momyi*momyi + momzi*momzi)/rhoi;

  fout[0] = rhoi;
  fout[1] = momxi; fout[2] = momyi; fout[3] = momzi;
  fout[4] = eri;     
}

void
evalFieldInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  struct burch_ctx *app = ctx;
  double pi = 3.141592653589793238462643383279502884;
  // Physical constants (using normalized code units).
  double gasGamma = app->gasGamma;
  double epsilon0 = app->epsilon0; // Permittivity of free space.
  double mu0 = app->mu0; // Permiability of free space.
  double lightSpeed = 1.0/sqrt(mu0*epsilon0); // Speed of light.
  double ionMass = app->ionMass; // Proton mass.
  double ionCharge = app->charge; // Proton charge.
  double elcMass = app->elcMass; // Electron mass.
  double elcCharge = -1.0*app->charge; // Electron charge.
  double n0 = app->n0;
  double B0 = app->B0;
  double vAe = B0/sqrt(mu0*n0*elcMass);

  double w0 = app->w0;
  double psi0 = 0.1*B0*w0;
  double guide1 = 0.099*B0;
  double guide2 = 0.099*B0;

  double b1 = 1.696*B0;
  double b2 = 1.00*B0;
  double n1 = 0.062*n0;
  double n2 = 1.0*n0;
  double beta2 = app->beta;
  double T_i2 = beta2*(b2*b2) / (2.0*n2*mu0);
  double T_i1_T_i2 = 7.73 / 1.374;
  double T_e1_T_i2 = 1.288 / 1.374;

  double T_e1 = T_i2*T_e1_T_i2;
  double T_i1  = T_i2*T_i1_T_i2;

  // derived quantities (asymptotic values)
  double T_e2 = (0.5*(b1*b1-b2*b2)+0.5*(guide1*guide1-guide2*guide2)+n1*(T_i1+T_e1)-n2*T_i2)/n2; //so the system is in force balance

  double Lx = app->Lx;
  double Ly = app->Ly;

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
  struct burch_ctx ctx = create_ctx(); // context for init functions
  // electron/ion equations
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(ctx.gasGamma);
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(ctx.gasGamma);

  printf("Parallel viscosity = %lg\n", ctx.eta_par);
  printf("Perpendicular viscosity = %lg\n", ctx.eta_perp);
  printf("Ratio eta_par/eta_perp = %lg\n", ctx.eta_par/ctx.eta_perp);
  printf("Collision time = %lg\n", ctx.tau);
  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = -1.0*ctx.charge, .mass = ctx.elcMass,

    .equation = elc_euler,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalElcInit,
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.charge, .mass = ctx.ionMass,

    .equation = ion_euler,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalIonInit,
  };  

  // VM app
  struct gkyl_moment app_inp = {
    .name = "brag_burch",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { ctx.Lx, ctx.Ly }, 
    .cells = { 256, 128  },

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },
    .cfl_frac = 0.001,

    .num_species = 2,
    .species = { elc, ion },
    .type_brag = GKYL_BRAG_MAG_FULL,
    .coll_fac = ctx.coll_fac,
    .field = {
      .epsilon0 = ctx.epsilon0,
      .mu0 = ctx.mu0,
      
      .evolve = 1,
      .ctx = &ctx,
      .init = evalFieldInit,
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 1.0;

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
