#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm.h>
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
  double tend; // end time
  double min_dt; // minimum size time step for breaking out of simulation
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
    .vte = 0.05,
    .vti = ctx.vte/sqrt(ctx.massIon),
    .beta = 0.1, 
    .Lx = 2000.0,
    .n0 = 1.0,
    .tend = 25000.0,
    .min_dt = 0.001,
  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_pkpm_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) 
    gkyl_pkpm_app_write(app, tcurr, iot->curr-1);
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 1792);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 96);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct pkpm_brio_wu_ctx ctx = create_ctx(); // context for init functions
  
  // electrons
  struct gkyl_pkpm_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -12.0 * ctx.vte},
    .upper = { 12.0 * ctx.vte}, 
    .cells = { VX },

    .ctx_dist = &ctx,
    .ctx_fluid = &ctx,
    .init_dist = evalDistFuncElc,
    .init_fluid = evalFluidElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
    },    

    .diffusion = {.D = 1.0e-5, .order=4},
  }; 
  
  // ions
  struct gkyl_pkpm_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -24.0 * ctx.vti},
    .upper = { 24.0 * ctx.vti}, 
    .cells = { VX },

    .ctx_dist = &ctx,
    .ctx_fluid = &ctx,
    .init_dist = evalDistFuncIon,
    .init_fluid = evalFluidIon,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
    },    

    .diffusion = {.D = 1.0e-5, .order=4},
  };

  // field
  struct gkyl_pkpm_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc
  };

  // pkpm app
  struct gkyl_pkpm pkpm = {
    .name = "pkpm_brio_wu_p2",

    .cdim = 1, .vdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 2,
    .basis_type = app_args.basis_type, 

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&pkpm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.tend;
  double dt = tend-tcurr;
  int nframe = 250;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_pkpm_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_pkpm_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    if (status.dt_actual < ctx.min_dt) {
      printf("** Time step crashing! Aborting simulation and writing out last output ....\n");
      gkyl_pkpm_app_write(app, tcurr, 1000);
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  // fetch simulation statistics
  struct gkyl_pkpm_stat stat = gkyl_pkpm_app_stat(app);

  gkyl_pkpm_app_cout(app, stdout, "\n");
  gkyl_pkpm_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_pkpm_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_pkpm_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_pkpm_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_pkpm_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_pkpm_app_cout(app, stdout, "Fluid Species RHS calc took %g secs\n", stat.fluid_species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species PKPM Vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_pkpm_app_cout(app, stdout, "EM Variables (bvar) calculation took %g secs\n", stat.field_em_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm);
  gkyl_pkpm_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);

  gkyl_pkpm_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_pkpm_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  // simulation complete, free app
  gkyl_pkpm_app_release(app);
  
  return 0;
}
