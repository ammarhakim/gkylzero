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
  double nbOverN0;
  double vAe;
  double B0;
  double guide;
  double w0;
  double psi0;
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
maxwellian(double n, double v, double temp, double mass)
{
  double v2 = v*v;
  return n/sqrt(2.0*M_PI*temp/mass)*exp(-v2/(2.0*temp/mass));
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

  double me = app->massElc;
  double mi = app->massIon;
  double T_e = app->T_e;
  double T_i = app->T_i;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double n0 = app->n0;
  double nbOverN0 = app->nbOverN0;
  double B0 = app->B0;
  double w0 = app->w0;

  double n = n0*(sech2(y/w0) + nbOverN0);
  
  double fv = maxwellian(n, vx, T_e, me);
    
  fout[0] = fv;
  fout[1] = T_e/me*fv;
}
void
evalDistFuncIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_gem_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], vx = xn[2];

  double me = app->massElc;
  double mi = app->massIon;
  double T_e = app->T_e;
  double T_i = app->T_i;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double n0 = app->n0;
  double nbOverN0 = app->nbOverN0;
  double B0 = app->B0;
  double w0 = app->w0;

  double n = n0*(sech2(y/w0) + nbOverN0);
  
  double fv = maxwellian(n, vx, T_i, mi);
    
  fout[0] = fv;
  fout[1] = T_i/mi*fv;
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
  double Jz = -B0/w0*sech2(y/w0) - psi0*cos(pi_2*x/Lx)*cos(M_PI*y/Ly)*((pi_2/Lx)*(pi_2/Lx) + (M_PI/Ly)*(M_PI/Ly));

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
  double Jz  = -B0/w0*sech2(y/w0) - psi0*cos(pi_2*x/Lx)*cos(M_PI*y/Ly)*((pi_2/Lx)*(pi_2/Lx) + (M_PI/Ly)*(M_PI/Ly));

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

  double b1x = B0*tanh(y/w0);
  double b1y = 0.0;
  double b1z = guide;

  double E_x = 0.0;
  double E_y = 0.0;
  double E_z = 0.0;
  double B_x = b1x - psi0 * (M_PI / Ly) * cos(2 * M_PI * x / Lx) * sin(M_PI * y / Ly);
  double B_y = b1y + psi0 * (2 * M_PI / Lx) * sin(2 * M_PI * x / Lx) * cos(M_PI * y / Ly);
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
  double massIon = 100.0; // ion mass
  double chargeIon = 1.0; // ion charge

  double Te_Ti = 0.2; // ratio of electron to ion temperature
  double n0 = 1.0; // initial number density
  double nbOverN0 = 0.2; // background number density

  double vAe = 0.04;
  double vAi = vAe/sqrt(massIon);
  double beta_elc = 1.0/6.0;
  // B0 has 1/sqrt(2) factor because B is equally subdivided between
  // guide field and in-plane field
  double B0 = vAe*sqrt(mu0*n0*massElc)/sqrt(2.0);  
  double guide = B0;
  double tot_B = sqrt(B0*B0 + guide*guide);

  double T_e = beta_elc*tot_B*tot_B/(2.0*mu0*n0);
  double T_i = T_e/Te_Ti;
  double vtElc = sqrt(2.0*T_e/massElc);
  double vtIon = sqrt(2.0*T_i/massIon);

  // ion cyclotron frequency and gyroradius
  double omegaCi = chargeIon*tot_B/massIon;
  double di = vAi/omegaCi;

  // Layer width and perturbation
  double w0 = 0.5*di;
  double psi0 = 0.1*tot_B*di;

  // collision frequencies
  double nuElc = 0.01*omegaCi;
  double nuIon = 0.01*omegaCi/sqrt(massIon);

  // domain size and simulation time
  double Lx = 8.0*M_PI*di;
  double Ly = 4.0*M_PI*di;
  double tend = 35.0/omegaCi;
  
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
    .nbOverN0 = nbOverN0,
    .vAe = vAe,
    .B0 = B0,
    .guide = guide,
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 64);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 32);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 32);

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
    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
  };  
  
  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -4.0 * ctx.vtElc},
    .upper = { 4.0 * ctx.vtElc}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
    },    

    .num_diag_moments = 0,
    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
  };

  // ion momentum                                                                                              
  struct gkyl_vlasov_fluid_species fluid_ion = {
    .name = "fluid_ion",
    .num_eqn = 3,
    .pkpm_species = "ion",
    .ctx = &ctx,
    .init = evalFluidIon,
    .nuHyp = 1.0e-4, 
    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
  };  
  
  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -4.0 * ctx.vtIon},
    .upper = { 4.0 * ctx.vtIon}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncIon,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
    },    

    .num_diag_moments = 0,
    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,
    .bcy = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "pkpm_gem_p1",

    .cdim = 2, .vdim = 1,
    .lower = { -ctx.Lx/2.0, -ctx.Ly/2.0 },
    .upper = { ctx.Lx/2.0, ctx.Ly/2.0 },
    .cells = { NX, NY },
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    //.cfl_frac = 0.8,
    
    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

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
  int nframe = 35;
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
