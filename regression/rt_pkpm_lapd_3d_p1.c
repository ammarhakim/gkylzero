#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct pkpm_lapd_ctx {
  double epsilon0;
  double mu0;
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te_Ti; // electron to ion temperature ratio
  double n0;
  double di; // ion inertial length for setting density gradient length scale
  double vAe;
  double B0;
  double beta;
  double vtElc;
  double vtIon;
  double nuElc;
  double nuIon;
  double L;
  double Lz;
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

void
evalDistFuncElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_lapd_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], z = xn[2], vx = xn[3];
  // Radial density profile in x and y
  double sigma2 = app->di*app->di;
  double n = app->n0*exp(-(x*x + y*y)/(2.0*sigma2));
  double fv = maxwellian(n, vx, app->vtElc);
    
  fout[0] = fv;
  fout[1] = app->vtElc*app->vtElc*fv;
}
void
evalDistFuncIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_lapd_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], z = xn[2], vx = xn[3];
  // Radial density profile in x and y
  double sigma2 = app->di*app->di;
  double n = app->n0*exp(-(x*x + y*y)/(2.0*sigma2));
  double fv = maxwellian(n, vx, app->vtIon);
    
  fout[0] = fv;
  fout[1] = app->vtIon*app->vtIon*fv;
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_lapd_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], z = xn[2];
  // No initial flows
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalFluidIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_lapd_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], z = xn[2];
  // No initial flows
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_lapd_ctx *app = ctx;

  double x = xn[0], y = xn[1], z = xn[2];
  // No initial fields
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalExtEmFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_lapd_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double B_z = app->B0;
  // External B field is in the z direction
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = B_z;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_lapd_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_lapd_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = app->nuIon;
}

struct pkpm_lapd_ctx
create_ctx(void)
{
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // pemiability of free space

  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge
  double massIon = 1836.153; // ion mass
  double chargeIon = 1.0; // ion charge

  double Te_Ti = 1.0; // ratio of electron to ion temperature
  double n0 = 1.0; // initial number density

  double vAe = 0.5;
  double B0 = vAe*sqrt(mu0*n0*massElc);
  double beta = 0.01;
  double vtElc = vAe*sqrt(beta/2.0);

  // ion velocities
  double vAi = vAe/sqrt(massIon);
  double vtIon = vtElc/sqrt(massIon); //Ti/Te = 1.0

  // Cyclotron frequencies and inertial lengths
  double omegaCi = chargeIon*B0/massIon;
  double omegaCe = chargeIon*B0/massElc;
  double di = vAi/omegaCi;
  double de = vAe/omegaCe;

  // Plasma frequencies and Debye length
  double wpe = sqrt(chargeIon*chargeIon*n0/(epsilon0*massElc));
  double wpi = sqrt(chargeIon*chargeIon*n0/(epsilon0*massIon));
  double lambdaD = vtElc/wpe;

  // collision frequencies
  double nuElc = 0.001*omegaCi;
  double nuIon = 0.001*omegaCi/sqrt(massIon)*(Te_Ti*sqrt(Te_Ti));

  // domain size and simulation time
  double L = 2*di;
  double Lz = 4*di;
  double tend = 10.0/omegaCi;

  // Denormalized values for reference
  // Compute electron thermal velocity from speed of light normalization
  double vtElc_denorm = vtElc*299792458.0;
  // Compute electron temperature in ev from thermal velocity
  double Telc_denorm = vtElc_denorm*vtElc_denorm*9.109383632e-31;
  // Compute density in 1/m^3 from collisionality (and plasma parameter)
  // Note plasma parameter is computed as Lambda = wpe/nu_ee, which is missing a factor of 2pi/ln(Lambda) (1/2-1/3)
  double nElc_denorm = pow(8.85418781e-12*Telc_denorm, 3)/(wpe/nuElc*wpe/nuElc*pow(1.60217663e-19, 6));
  // Compute Debye length in meters from density and temperature
  double lambdaD_denorm = sqrt(8.85418781e-12*Telc_denorm/(1.60217663e-19*1.60217663e-19*nElc_denorm));
  // Compute background magnetic field in Tesla from beta, density, and temperature
  double B0_denorm = sqrt(2.0*1.256637062e-6*nElc_denorm*Telc_denorm/beta);
  // Compute axial extent in meters
  double axial_denorm = Lz/lambdaD*lambdaD_denorm;
  // Compute radial extent in meters
  double radial_denorm = L/lambdaD*lambdaD_denorm;

  printf("Axial direction is %g Debye lengths\n", Lz/lambdaD);
  printf("Electron Plasma parameter (wpe/nu_ee) = %g\n", wpe/nuElc);
  printf("Electron thermal velocity in m/s = %g\n", vtElc_denorm);
  printf("Electron temperature in ev = %g\n", Telc_denorm/1.60217663e-19);
  printf("Electron density in 1/m^3 = %g\n", nElc_denorm);
  printf("Electron Debye length in meters = %g\n", lambdaD_denorm);
  printf("Background magnetic field in Tesla = %g\n", B0_denorm);
  printf("Axial extent in meters = %g\n", axial_denorm);
  printf("Radial extent in meters = %g\n", radial_denorm);
  
  struct pkpm_lapd_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .chargeElc = chargeElc,
    .massElc = massElc,
    .chargeIon = chargeIon,
    .massIon = massIon,
    .Te_Ti = Te_Ti,
    .n0 = n0,
    .di = di,
    .vAe = vAe,
    .B0 = B0,
    .beta = beta,
    .vtElc = vtElc,
    .vtIon = vtIon,
    .nuElc = nuElc,
    .nuIon = nuIon,
    .L = L,
    .Lz = Lz,
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 16);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[0], 32);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 16);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
     
  struct pkpm_lapd_ctx ctx = create_ctx(); // context for init functions

  // electron momentum                                   
  struct gkyl_vlasov_fluid_species fluid_elc = {
    .name = "fluid_elc",
    .num_eqn = 3,
    .pkpm_species = "elc",
    .ctx = &ctx,
    .init = evalFluidElc,
    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
    .bcy = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
    .bcz = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
  };  
  
  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -6.0 * ctx.vtElc},
    .upper = { 6.0 * ctx.vtElc}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
      //.normNu = true,
    },    

    .num_diag_moments = 0,
    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
    .bcy = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
    .bcz = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
  };

  // ion momentum
  struct gkyl_vlasov_fluid_species fluid_ion = {
    .name = "fluid_ion",
    .num_eqn = 3,
    .pkpm_species = "ion",
    .ctx = &ctx,
    .init = evalFluidIon,
    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
    .bcy = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
    .bcz = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
  };  
  
  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -6.0 * ctx.vtIon},
    .upper = { 6.0 * ctx.vtIon}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncIon,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
      //.normNu = true,
    },    

    .num_diag_moments = 0,
    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
    .bcy = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
    .bcz = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,
    // Plasma EM field BCs are PEC, external field goes into conducting wall
    .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL }, 
    .bcy = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL }, 
    .bcz = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL }, 
    .ext_em = evalExtEmFunc,
    .ext_em_ctx = &ctx,
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "pkpm_3d_lapd_p1",

    .cdim = 3, .vdim = 1,
    .lower = { -ctx.L, -ctx.L, 0.0 },
    .upper = { ctx.L, ctx.L, ctx.Lz },
    .cells = { NX, NX, NZ },
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    //.cfl_frac = 0.8,
    
    .num_periodic_dir = 0,
    .periodic_dirs = { },

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

  gkyl_vlasov_app_cout(app, stdout, "\n");
  gkyl_vlasov_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_vlasov_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_vlasov_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_vlasov_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_vlasov_app_cout(app, stdout, "Fluid Species RHS calc took %g secs\n", stat.fluid_species_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species PKPM Vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_vlasov_app_cout(app, stdout, "EM Variables (bvar) calculation took %g secs\n", stat.field_em_vars_tm);
  gkyl_vlasov_app_cout(app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm);
  gkyl_vlasov_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  // simulation complete, free app
  gkyl_vlasov_app_release(app);
  
  return 0;
}
