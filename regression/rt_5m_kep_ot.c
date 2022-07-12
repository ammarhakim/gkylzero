#include <gkyl_alloc.h>
#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

struct ot_5m_ctx {
  double epsilon0;
  double mu0;
  double n; // initial number density

  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass

  double vte; // electron thermal velocity
  double vti; // ion thermal velocity
  double beta; // plasma beta
  double B0; // reference magnetic field
  double vA; // Alfven speed
  double omegaCe; // electron cyclotron frequency
  double omegaCi; // ion cyclotron frequency
  double larmor_e; // electron gyroradius
  double larmor_i; // ion gyroradius
  double L; // size of the box

  double lambdaD;
  double wpe;
  double Lambda;

  double taue; // electron-electron collision time
  double taui; // ion-ion collision time
};

void
evalElcInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  struct ot_5m_ctx *app = ctx;

  double pi = M_PI;
  double eps0 = app->epsilon0;
  double mu0 = app->mu0;
  double gasGamma = 5.0/3.0;
  double elcMass = app->massElc;
  double ionMass = app->massIon;
  double elcCharge = app->chargeElc;
  double ionCharge = app->chargeIon;

  double n0 = app->n;

  double beta = app->beta;
  double B0 = app->B0;
  double vtElc = app->vte;
  double vAi = app->vA;
  double vtIon = app->vti;
  double omegaCi = app->omegaCi;
  double larmor_i = app->larmor_i;

  double u0x = 0.2*vAi;
  double u0y = 0.2*vAi;
  double B0x = 0.2*B0;
  double B0y = 0.2*B0;

  double Lx = app->L;
  double Ly = app->L;

  double _2pi = 2.0*pi;
  double _4pi = 2.0*_2pi;

  double Jz = (B0y*(_4pi/Lx)*cos(_4pi*x/Lx) + B0x*(_2pi/Ly)*cos(_2pi*y/Ly)) / mu0;
  
  double vdrift_x = -u0x*sin(_2pi*y/Ly);
  double vdrift_y = u0y*sin(_2pi*x/Lx);
  double vdrift_z = -Jz / (n0*ionCharge);

  double rhoe = n0*elcMass;
  double exmom = vdrift_x*rhoe;
  double eymom = vdrift_y*rhoe;
  double ezmom = vdrift_z*rhoe;
  double ere = n0*vtElc*vtElc*elcMass/(gasGamma-1) + 0.5*exmom*exmom/rhoe + 0.5*eymom*eymom/rhoe + 0.5*ezmom*ezmom/rhoe;

  fout[0] = rhoe;
  fout[1] = exmom; fout[2] = eymom; fout[3] = ezmom;
  fout[4] = ere;  
}

void
evalIonInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  struct ot_5m_ctx *app = ctx;

  double pi = M_PI;
  double eps0 = app->epsilon0;
  double mu0 = app->mu0;
  double gasGamma = 5.0/3.0;
  double elcMass = app->massElc;
  double ionMass = app->massIon;
  double elcCharge = app->chargeElc;
  double ionCharge = app->chargeIon;

  double n0 = app->n;

  double beta = app->beta;
  double B0 = app->B0;
  double vtElc = app->vte;
  double vAi = app->vA;
  double vtIon = app->vti;
  double omegaCi = app->omegaCi;
  double larmor_i = app->larmor_i;

  double u0x = 0.2*vAi;
  double u0y = 0.2*vAi;
  double B0x = 0.2*B0;
  double B0y = 0.2*B0;

  double Lx = app->L;
  double Ly = app->L;

  double _2pi = 2.0*pi;
  double _4pi = 2.0*_2pi;

  double Jz = (B0y*(_4pi/Lx)*cos(_4pi*x/Lx) + B0x*(_2pi/Ly)*cos(_2pi*y/Ly)) / mu0;
  
  double vdrift_x = -u0x*sin(_2pi*y/Ly);
  double vdrift_y = u0y*sin(_2pi*x/Lx);
  double vdrift_z = 0.;

  double rhoi = n0*ionMass;
  double ixmom = vdrift_x*rhoi;
  double iymom = vdrift_y*rhoi;
  double izmom = vdrift_z*rhoi;
  double eri = n0*vtIon*vtIon*ionMass/(gasGamma-1) + 0.5*ixmom*ixmom/rhoi + 0.5*iymom*iymom/rhoi + 0.5*izmom*izmom/rhoi;

  fout[0] = rhoi;
  fout[1] = ixmom; fout[2] = iymom; fout[3] = izmom;
  fout[4] = eri;    
}

void
evalFieldInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  struct ot_5m_ctx *app = ctx;

  double pi = M_PI;
  double eps0 = app->epsilon0;
  double mu0 = app->mu0;
  double gasGamma = 5.0/3.0;
  double elcMass = app->massElc;
  double ionMass = app->massIon;
  double elcCharge = app->chargeElc;
  double ionCharge = app->chargeIon;

  double n0 = app->n;

  double beta = app->beta;
  double B0 = app->B0;
  double vtElc = app->vte;
  double vAi = app->vA;
  double vtIon = app->vti;
  double omegaCi = app->omegaCi;
  double larmor_i = app->larmor_i;

  double u0x = 0.2*vAi;
  double u0y = 0.2*vAi;
  double B0x = 0.2*B0;
  double B0y = 0.2*B0;

  double Lx = app->L;
  double Ly = app->L;

  double _2pi = 2.0*pi;
  double _4pi = 2.0*_2pi;

  double Jz = (B0y*(_4pi/Lx)*cos(_4pi*x/Lx) + B0x*(_2pi/Ly)*cos(_2pi*y/Ly)) / mu0;

  double Bx = -B0x*sin(_2pi*y/Ly);
  double By = B0y*sin(_4pi*x/Lx);
  double Bz = B0;

  // Assumes ni = ne = n0
  double u_xe = -u0x*sin(_2pi*y/Ly);
  double u_ye = u0y*sin(_2pi*x/Lx);
  double u_ze = -Jz / (ionCharge*n0);

  // E = - v_e x B ~  (J - u) x B
  double Ex = - (u_ye*Bz - u_ze*By);
  double Ey = - (u_ze*Bx - u_xe*Bz);
  double Ez = - (u_xe*By - u_ye*Bx);


  // electric field
  fout[0] = Ex, fout[1] = Ey; fout[2] = Ez;
  // magnetic field
  fout[3] = Bx, fout[4] = By; fout[5] = Bz;

  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
}

struct ot_5m_ctx
create_ctx(void)
{
  double epsilon0 = 1.0;
  double mu0 = 1.0;
  double massElc = 1.0;
  double massIon = 25.0;
  double charge = 1.0;
  double n0 = 1.0;
  double vte = 1.0/4.0;
  double vti = vte/sqrt(massIon/massElc); // equal proton and electron temperatures

  double beta = 1.0;
  double B0 = sqrt(mu0*n0*vte*vte*massElc/beta);
  double vA = B0/sqrt(mu0*n0*massIon);

  double omegaCe = charge*B0/massElc;
  double omegaCi = charge*B0/massIon;

  double larmor_e = vte/omegaCe;
  double larmor_i = vti/omegaCi;

  double L = 16.0*M_PI*larmor_i;

  double taue = 6.0*sqrt(2.0*M_PI*M_PI*M_PI)*massElc*massElc*vte*vte*vte*epsilon0*epsilon0/(charge*charge*charge*charge*n0);
  double taui = 6.0*sqrt(2.0*M_PI*M_PI*M_PI)*massIon*massIon*vti*vti*vti*epsilon0*epsilon0/(charge*charge*charge*charge*n0);

  struct ot_5m_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .n = n0,
    .chargeElc = -charge,
    .massElc = massElc,
    .chargeIon = charge,
    .massIon = massIon,

    .vte = vte, 
    .vti = vti, 

    .beta = beta,
    .B0 = B0,
    .vA = vA,

    .omegaCe = omegaCe,
    .omegaCi = omegaCi,
    .larmor_e = larmor_e,
    .larmor_i = larmor_i,

    .taue = taue,
    .taui = taui,

    .L = L,    
  };
  return ctx;
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 256);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 256);

  struct ot_5m_ctx ctx = create_ctx(); // context for init functions
  printf("Ion cyclotron frequency = %lg\n", ctx.omegaCi);

  printf("Electron-electron collision time = %lg\n", ctx.taue);
  printf("Ion-ion collision time = %lg\n", ctx.taui);

  printf("Braginskii expansion parameter tau*omegaC, electrons = %lg\n", ctx.taue*ctx.omegaCe);
  printf("Braginskii expansion parameter tau*omegaC, ions = %lg\n", ctx.taui*ctx.omegaCi);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  // electron/ion equations
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(5.0/3.0);
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(5.0/3.0);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,

    .equation = elc_euler,
    .evolve = 1,
    .init = evalElcInit,
    .ctx = &ctx,

    .type_brag = GKYL_BRAG_UNMAG_FULL,
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,

    .equation = ion_euler,
    .evolve = 1,
    .init = evalIonInit,
    .ctx = &ctx,

    .type_brag = GKYL_BRAG_UNMAG_FULL,
  };  

  // VM app
  struct gkyl_moment app_inp = {
    .name = "5m_kep_ot",

    .ndim = 2,
    .lower = { 0., 0. },
    .upper = { ctx.L, ctx.L }, 
    .cells = { NX, NY },

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },
    .cfl_frac = 1.0,

    .fluid_scheme = GKYL_MOMENT_FLUID_KEP,

    .num_species = 2,
    .species = { elc, ion },
    .coll_fac = 1.0,
    .field = {
      .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
      .mag_error_speed_fact = 1.0,
      
      .evolve = 1,
      .init = evalFieldInit,
      .ctx = &ctx,
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 100.0/ctx.omegaCi;
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
  gkyl_wv_eqn_release(elc_euler);
  gkyl_wv_eqn_release(ion_euler);
  gkyl_moment_app_release(app);  
  
  return 0;
}
