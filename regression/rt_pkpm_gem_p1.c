#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct pkpm_gem_ctx {
  double epsilon0;
  double mu0;
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te_Ti; // electron to ion temperature ratio
  double n0;
  double vAe;
  double B0;
  double guide;
  double w0;
  double psi0;
  double beta;
  double T_e;
  double T_i;
  double vtElc;
  double vtIon;
  double nuElc;
  double nuIon;
  double Lx;
  double Ly;
  double tend;
  double min_dt;
  bool use_gpu;
};

static inline double
maxwellian(double n, double v, double vth)
{
  double v2 = v*v;
  return n/sqrt(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

static inline double
sech2(double x)
{
  return 1.0/(cosh(x)*cosh(x));
}

void
evalDistFuncElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], vx = xn[2];

  double vt = app->vtElc;
  double T_e = app->T_e;
  double T_i = app->T_i;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double n0 = app->n0;
  double B0 = app->B0;
  double w0 = app->w0;

  double b1x = B0*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1);

  double n = (0.5*(B0*B0 - b1x*b1x) + n0*(T_i+T_e))/(T_i+T_e);
  
  double fv = maxwellian(n, vx, vt);
    
  fout[0] = fv;
  fout[1] = vt*vt*fv;
}
void
evalDistFuncIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], vx = xn[2];

  double vt = app->vtIon;
  double T_e = app->T_e;
  double T_i = app->T_i;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double n0 = app->n0;
  double B0 = app->B0;
  double w0 = app->w0;

  double b1x = B0*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1);

  double n = (0.5*(B0*B0 - b1x*b1x) + n0*(T_i+T_e))/(T_i+T_e);
  
  double fv = maxwellian(n, vx, vt);
    
  fout[0] = fv;
  fout[1] = vt*vt*fv;
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;
  
  double x = xn[0], y = xn[1];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double B0 = app->B0;
  double w0 = app->w0;
  double psi0 = app->psi0;  
  double T_e = app->T_e;
  double T_i = app->T_i;

  double TeFrac = T_e/(T_i+T_e);

  double pi_2 = 2.0*M_PI;
  double pi_4 = 4.0*M_PI;

  double Jx = 0.0;
  double Jy = 0.0;
  double Jz  = -B0/w0*(sech2((y-Ly*.25)/w0) - sech2((y-Ly*.75)/w0)+ sech2((y-Ly*1.25)/w0) - sech2((y+Ly*.25)/w0)) 
    - psi0*sin(pi_2*x/Lx)*((pi_2/Lx)*(pi_2/Lx)*(1 - cos(pi_4*y/Ly)) + (pi_4/Ly)*(pi_4/Ly)*cos(pi_4*y/Ly));

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = me*Jz/qe*TeFrac;
}

void
evalFluidIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;
  
  double x = xn[0], y = xn[1];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double B0 = app->B0;
  double w0 = app->w0;
  double psi0 = app->psi0;  
  double T_e = app->T_e;
  double T_i = app->T_i;

  double TiFrac = T_i/(T_i+T_e);

  double pi_2 = 2.0*M_PI;
  double pi_4 = 4.0*M_PI;

  double Jx = 0.0;
  double Jy = 0.0;
  double Jz  = -B0/w0*(sech2((y-Ly*.25)/w0) - sech2((y-Ly*.75)/w0)+ sech2((y-Ly*1.25)/w0) - sech2((y+Ly*.25)/w0)) 
    - psi0*sin(pi_2*x/Lx)*((pi_2/Lx)*(pi_2/Lx)*(1 - cos(pi_4*y/Ly)) + (pi_4/Ly)*(pi_4/Ly)*cos(pi_4*y/Ly));

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = mi*Jz/qi*TiFrac;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;

  double x = xn[0], y = xn[1];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double B0 = app->B0;
  double guide = app->guide;
  double w0 = app->w0;
  double psi0 = app->psi0; 

  double pi_2 = 2.0*M_PI;
  double pi_4 = 4.0*M_PI;

  double b1x = B0*(tanh((y-Ly*.25)/w0)-tanh((y-Ly*.75)/w0) + tanh((y-Ly*1.25)/w0)-tanh((y+Ly*.25)/w0)+1);
  double b1y = 0.0;
  double b1z = guide;

  double E_x = 0.0;
  double E_y = 0.0;
  double E_z = 0.0;
  double B_x = b1x - psi0*pi_4/Ly*sin(pi_2*x/Lx)*sin(pi_4*y/Ly);
  double B_y = b1y + psi0*pi_2/Lx*cos(pi_2*x/Lx)*(1-cos(pi_4*y/Ly));
  double B_z = b1z;
  
  fout[0] = E_x; fout[1] = E_y, fout[2] = E_z;
  fout[3] = B_x; fout[4] = B_y; fout[5] = B_z;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;
  fout[0] = app->nuIon;
}

struct pkpm_gem_ctx
create_ctx(void)
{
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // pemiability of free space

  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge
  double massIon = 25.0; // ion mass
  double chargeIon = 1.0; // ion charge

  double Te_Ti = 0.2; // ratio of electron to ion temperature
  double n0 = 1.0; // initial number density

  double vAe = 1.0/4.0;
  double B0 = vAe*sqrt(mu0*n0*massElc);
  double beta = 5.0/6.0;
  double vtElc = vAe*sqrt(beta*Te_Ti);

  double guide = 0.1*B0;

  double T_e = vtElc*vtElc/2.0;
  double T_i = T_e/Te_Ti;

  // ion velocities
  double vAi = vAe/sqrt(massIon);
  double vtIon = vtElc/sqrt(massIon); //Ti/Te = 1.0

  // ion cyclotron frequency and gyroradius
  double omegaCi = chargeIon*B0/massIon;
  double di = vAi/omegaCi;

  // Layer width and perturbation
  double w0 = 0.5*di;
  double psi0 = 0.1*B0*di;

  // collision frequencies
  double nuElc = 0.01*omegaCi;
  double nuIon = 0.01*omegaCi/sqrt(massIon);

  // domain size and simulation time
  double Lx = 8.0*M_PI*di;
  double Ly = 8.0*M_PI*di;
  double tend = 10.0/omegaCi;
  
  struct pkpm_gem_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .chargeElc = chargeElc,
    .massElc = massElc,
    .chargeIon = chargeIon,
    .massIon = massIon,
    .Te_Ti = Te_Ti,
    .T_e = T_e,
    .T_i = T_i,
    .n0 = n0,
    .vAe = vAe,
    .B0 = B0,
    .guide = guide,
    .beta = beta,
    .w0 = w0,
    .psi0 = psi0,
    .vtElc = vtElc,
    .vtIon = vtIon,
    .nuElc = nuElc,
    .nuIon = nuIon,
    .Lx = Lx,
    .Ly = Ly,
    .tend = tend,
    .min_dt = 1.0e-2, 
  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_vlasov_app_write(app, tcurr, iot->curr-1);
    gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, iot->curr-1);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 128);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 128);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 16);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
     
  struct pkpm_gem_ctx ctx = create_ctx(); // context for init functions

  // electron momentum                                                                                              
  struct gkyl_vlasov_fluid_species fluid_elc = {
    .name = "fluid_elc",
    .num_eqn = 3,
    .pkpm_species = "elc",
    .ctx = &ctx,
    .init = evalFluidElc,
    .nuHyp = 1.0e-4,
  };  
  
  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -8.0 * ctx.vtElc},
    .upper = { 8.0 * ctx.vtElc}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
    },    

    .num_diag_moments = 0,
  };

  // ion momentum                                                                                              
  struct gkyl_vlasov_fluid_species fluid_ion = {
    .name = "fluid_ion",
    .num_eqn = 3,
    .pkpm_species = "ion",
    .ctx = &ctx,
    .init = evalFluidIon,
    .nuHyp = 1.0e-4, 
  };  
  
  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -8.0 * ctx.vtIon},
    .upper = { 8.0 * ctx.vtIon}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncIon,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
    },    

    .num_diag_moments = 0,
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "pkpm_gem_p1",

    .cdim = 2, .vdim = 1,
    .lower = { 0.0, 0.0 },
    .upper = { ctx.Lx, ctx.Ly },
    .cells = { NX, NY },
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    //.cfl_frac = 0.8,
    
    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },

    .num_species = 2,
    .species = { elc, ion },
    .num_fluid_species = 2,
    .fluid_species = { fluid_elc, fluid_ion },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.tend;
  double dt = tend-tcurr;
  int nframe = 100;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    if (status.dt_actual < ctx.min_dt) {
      printf("** Time step crashing! Aborting simulation and writing out last output ....\n");
      gkyl_vlasov_app_write(app, tcurr, 1000);
      gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 1000);
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_vlasov_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  // simulation complete, free app
  gkyl_vlasov_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of forward-Euler calls %ld\n", stat.nfeuler);
  printf("Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    printf("Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    printf("Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  printf("Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  printf("Species RHS calc took %g secs\n", stat.species_rhs_tm);
  printf("Species collisions took %g secs\n", stat.species_coll_mom_tm);
  printf("Species collisions took %g secs\n", stat.species_coll_tm);
  printf("Field RHS calc took %g secs\n", stat.field_rhs_tm);
  printf("Current evaluation and accumulate took %g secs\n", stat.current_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
