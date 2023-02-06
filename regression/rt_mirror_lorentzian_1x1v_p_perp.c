#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct mirror_ctx {
  double epsilon0;
  double mu0;

  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass

  double vte; // electron thermal velocity
  double vti; // ion thermal velocity
  double Lx; // size of the box
  double n0; // initial number density

  double lambdaD;
  double wpe;
  double Lambda;

  double gamma; // FWHM of Lorentzian
  double loc; // location of Lorentzian
  double mag; // magnitude of Lorentzian

  double b_z0; // magnetic field at z = 0, for parameter checking
  double b_zL; // magnetic field at z = Lx, for parameter checking
};

static inline double sq(double x) { return x*x; }
static inline double cu(double x) { return x*x*x; }

double calcMagB(const double * GKYL_RESTRICT xn, void *ctx)
{
  struct mirror_ctx *app = ctx;
  double x = xn[0];
  double gamma = app->gamma;
  double loc = app->loc;
  double mag = app->mag;
  double magB = mag/(gamma*(1.0 + sq((x-loc)/gamma))) + mag/(gamma*(1.0 + sq((x+loc)/gamma)));
  return magB;
}



void
evalDistFuncElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mirror_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vte, n0 = app->n0;
  double fv = n0*exp(-sq(x))/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  double magB = calcMagB(xn,ctx);
  fout[0] = fv/magB;
}

void
evalDistFuncIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mirror_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vti, n0 = app->n0;
  double fv = n0*exp(-sq(x))/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  double magB = calcMagB(xn,ctx);
  fout[0] = fv/magB;
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mirror_ctx *app = ctx;
  double x = xn[0];
  double vt = app->vte, n = app->n0;
  double magB = calcMagB(xn,ctx);
  // No initial flow (u = 0)
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = n*exp(-sq(x))*vt*vt/magB;
}

void
evalFluidIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mirror_ctx *app = ctx;
  double x = xn[0];
  double vt = app->vti, n = app->n0;
  double magB = calcMagB(xn,ctx);
  // No initial flow (u = 0)
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = n*exp(-sq(x))*vt*vt/magB;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mirror_ctx *app = ctx;
  double x = xn[0];
  double B_x = 1.0;
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = B_x; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mirror_ctx *app = ctx;
  fout[0] = app->wpe/app->Lambda*log(app->Lambda);
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mirror_ctx *app = ctx;
  double Te_Ti = app->vte*app->vte*app->massElc/(app->vti*app->vti*app->massIon);
  double nu_ee = app->wpe/app->Lambda*log(app->Lambda);
  fout[0] = nu_ee/sqrt(app->massIon/app->massElc)*(Te_Ti*sqrt(Te_Ti));
}

void
evalMagB(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mirror_ctx *app = ctx;
  double x = xn[0];
  double magB = calcMagB(xn,ctx);
  fout[0] = magB;
}

struct mirror_ctx
create_ctx(void)
{
  double epsilon0 = 8.8541878128e-12;
  double massElc = 9.109384e-31;
  double charge = 1.602177e-19;
  double n0 = 3e19;
  double vte = 1.0e7;
  double lambdaD = sqrt(epsilon0*vte*vte*massElc/(charge*charge*n0));

  double loc = 1.0;
  double gamma = 0.1;
  double mag = 1.0;
  double b_z0 = mag/(gamma*(1.0 + sq((loc)/gamma))) + mag/(gamma*(1.0 + sq((loc)/gamma)));
  double Lx = 2.0;
  double b_zL = mag/(gamma*(1.0 + sq((Lx + loc)/gamma))) + mag/(gamma*(1.0 + sq((Lx - loc)/gamma)));
  struct mirror_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = 1.256637062e-6,
    .n0 = n0,
    .chargeElc = -charge,
    .massElc = massElc,
    .chargeIon = charge,
    .massIon = 3.34449469e-27,

    .vte = vte, //~1 kev electrons
    .vti = 7.0e5, //~10 kev deuterium
    // total length is 2 m
    // Debye length is 0.04 mm ~ 25k Debye lengths
    // Plasma parameter ~ 4e6
    // Electron plasma frequency ~ 2e11/s
    .lambdaD = lambdaD,
    .wpe = vte/lambdaD,
    .Lambda = n0*lambdaD*lambdaD*lambdaD,
    .Lx = Lx,

    .gamma = gamma,
    .loc = loc,
    .mag = mag,
    .b_z0 = b_z0,
    .b_zL = b_zL
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

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  
  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 256);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 64);  
  
  struct mirror_ctx ctx = create_ctx(); // context for init functions
  printf("Debye length = %lg\n", ctx.lambdaD);
  printf("Plasma frequency = %lg\n", ctx.wpe);
  printf("Plasma Parameter = %lg\n", ctx.Lambda);
  printf("Electron-Electron collision frequency = %lg\n", ctx.wpe/ctx.Lambda*log(ctx.Lambda));
  // electron beta at z = 0
  printf("Electron beta = %lg\n", 2.0*ctx.mu0*ctx.n0*ctx.vte*ctx.vte*ctx.massElc/(ctx.b_z0*ctx.b_z0));
  // electron gyroradius at z = 0
  printf("Electron gyroradius = %lg\n", ctx.vte/(ctx.chargeIon*ctx.b_z0/ctx.massElc));

  // magnetic field at z = Lx
  printf("Magnetic field at z = Lx, %lg\n", ctx.b_zL);
  // Debye length at z = Lx
  printf("Debye length at z = Lx, %lg\n", ctx.lambdaD/sqrt(exp(-sq(ctx.Lx))));
  
  // electron Pperp
  struct gkyl_vlasov_fluid_species fluid_elc = {
    .name = "fluid_elc",
    .num_eqn = 4,
    .pkpm_species = "elc",
    .ctx = &ctx,
    .init = evalFluidElc,

    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
  };  
  
  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -6.0 * ctx.vte},
    .upper = { 6.0 * ctx.vte}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncElc,

    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
    
    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
    },    

    .magB_ctx = &ctx,
    .magB = evalMagB,

    .num_diag_moments = 0,
  };

  // ion Pperp                                                                                              
  struct gkyl_vlasov_fluid_species fluid_ion = {
    .name = "fluid_ion",
    .num_eqn = 4,
    .pkpm_species = "ion",
    .ctx = &ctx,
    .init = evalFluidIon,

    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
  };  
  
  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -6.0 * ctx.vti},
    .upper = { 6.0 * ctx.vti}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncIon,

    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
    
    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
    },    

    .magB_ctx = &ctx,
    .magB = evalMagB,

    .num_diag_moments = 0,
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,

    .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL }
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "mirror_lorentzian_1x1v_p_perp",

    .cdim = 1, .vdim = 1,
    .lower = { -ctx.Lx },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

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
  double tcurr = 0.0, tend = 1.0e-7;
  double dt = tend-tcurr;
  int nframe = 1;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_vlasov_app_calc_integrated_mom(app, tcurr);
  gkyl_vlasov_app_calc_field_energy(app, tcurr);
  
  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);

    gkyl_vlasov_app_calc_integrated_mom(app, tcurr);
    gkyl_vlasov_app_calc_field_energy(app, tcurr);    
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_vlasov_app_write_integrated_mom(app);
  gkyl_vlasov_app_write_field_energy(app);  
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
  printf("Species Collisions calc took %g secs\n", stat.species_coll_tm);
  printf("Field RHS calc took %g secs\n", stat.field_rhs_tm);
  printf("Current evaluation and accumulate took %g secs\n", stat.current_tm);
  printf("Collisional moment calc took %g secs\n", stat.species_coll_mom_tm);
  printf("Fluid Species RHS calc took %g secs\n", stat.fluid_species_rhs_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
