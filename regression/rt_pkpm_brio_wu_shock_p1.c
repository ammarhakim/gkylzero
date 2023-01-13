#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct pkpm_brio_wu_ctx {
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double vte; // electron thermal velocity
  double vti; // ion thermal velocity
  double beta; // plasma beta (for strength of magnetic field)
  double Lx; // size of the box
  double n0; // initial number density
};

static inline double sq(double x) { return x*x; }

static inline double
maxwellian(double n, double v, double u, double vth)
{
  double v2 = (v - u)*(v - u);
  return n/sqrt(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

void
evalDistFuncElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_brio_wu_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vte, n0 = app->n0;
  double mass = app->massElc;
  double Lx = app->Lx;
  if (x<0.5*Lx) {
    fout[0] = maxwellian(n0, v, 0.0, vt);
    fout[1] = vt*vt*maxwellian(n0, v, 0.0, vt);
  }
  else {
    fout[0] = maxwellian(0.125*n0, v, 0.0, sqrt(0.1*vt*vt/0.125));
    fout[1] = 0.1*vt*vt/0.125*maxwellian(0.125*n0, v, 0.0, sqrt(0.1*vt*vt/0.125));
  }
}

void
evalDistFuncIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_brio_wu_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vti, n0 = app->n0;
  double mass = app->massIon;
  double Lx = app->Lx;
  if (x<0.5*Lx) {
    fout[0] = maxwellian(n0, v, 0.0, vt);
    fout[1] = vt*vt*maxwellian(n0, v, 0.0, vt);
  }
  else {
    fout[0] = maxwellian(0.125*n0, v, 0.0, sqrt(0.1*vt*vt/0.125));
    fout[1] = 0.1*vt*vt/0.125*maxwellian(0.125*n0, v, 0.0, sqrt(0.1*vt*vt/0.125));
  }
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_brio_wu_ctx *app = ctx;
  double x = xn[0];
  double vt = app->vte, n0 = app->n0;
  double mass = app->massElc;
  double Lx = app->Lx;
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalFluidIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_brio_wu_ctx *app = ctx;
  double x = xn[0];
  double vt = app->vti, n0 = app->n0;
  double mass = app->massIon;
  double Lx = app->Lx;
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_brio_wu_ctx *app = ctx;
  double x = xn[0];
  double beta = app->beta;
  double vt = app->vte;
  double mass = app->massElc;
  double n0 = app->n0;
  double Lx = app->Lx;
  double Telc = vt*vt*mass;
  double B0 = sqrt(2.0*n0*Telc/beta);
  double B_x = B0/sqrt(2.0);
  double B_z = 0.0;
  if (x<0.5*Lx)
    B_z = B0/sqrt(2.0);
  else
    B_z = -B0/sqrt(2.0);
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = B_x; fout[4] = 0.0; fout[5] = B_z;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_brio_wu_ctx *app = ctx;
  fout[0] = 4.0e-5;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_brio_wu_ctx *app = ctx;
  fout[0] = 4.0e-5/sqrt(app->massIon);
}

struct pkpm_brio_wu_ctx
create_ctx(void)
{
  struct pkpm_brio_wu_ctx ctx = {
    .chargeElc = -1.0,
    .massElc = 1.0,
    .chargeIon = 1.0,
    .massIon = 1836.153,
    .vte = 0.04,
    .vti = ctx.vte/sqrt(ctx.massIon),
    .beta = 0.1, 
    .Lx = 409.6,
    .n0 = 1.0
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 448);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 64);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct pkpm_brio_wu_ctx ctx = create_ctx(); // context for init functions

  // electron Tperp                                                                                              
  struct gkyl_vlasov_fluid_species fluid_elc = {
    .name = "fluid_elc",
    .num_eqn = 3,
    .pkpm_species = "elc",
    .ctx = &ctx,
    .init = evalFluidElc,
    .nuHyp = 1.0e-5,
  };  
  
  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -16.0 * ctx.vte},
    .upper = { 16.0 * ctx.vte}, 
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

  // ion Tperp                                                                                              
  struct gkyl_vlasov_fluid_species fluid_ion = {
    .name = "fluid_ion",
    .num_eqn = 3,
    .pkpm_species = "ion",
    .ctx = &ctx,
    .init = evalFluidIon,
    .nuHyp = 1.0e-5,
  };  
  
  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -16.0 * ctx.vti},
    .upper = { 16.0 * ctx.vti}, 
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
    .name = "pkpm_brio_wu_p1",

    .cdim = 1, .vdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 1,
    .basis_type = app_args.basis_type, 
    .cfl_frac = 0.8,

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
  double tcurr = 0.0, tend = 10000.0;
  double dt = tend-tcurr;

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  
  gkyl_vlasov_app_write(app, tcurr, 0);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 0);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    step += 1;
  }

  gkyl_vlasov_app_write(app, tcurr, 1);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 1);
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
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
