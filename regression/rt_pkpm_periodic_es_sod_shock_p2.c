#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm.h>
#include <rt_arg_parse.h>

struct pkpm_es_sod_shock_ctx {
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te_Ti; // electron to ion temperature ratio
  double vte; // electron thermal velocity
  double vti; // ion thermal velocity
  double n0; // initial number density
  double B0; // reference magnetic field 
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double Lx; // Domain size (x-direction).
  int Nx; // Cell count (x-direction).
  double t_end; // Final simulation time.
  double init_dt; // Initial time step guess so first step does not generate NaN
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
  struct pkpm_es_sod_shock_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double n_l = app->n0;
  double n_r = 0.125*n_l;
  double vt_l = app->vte;
  double vt_r = sqrt(0.1/0.125)*vt_l;
  if (fabs(x)<0.5*app->Lx) {
    fout[0] = maxwellian(n_l, v, 0.0, vt_l);
    fout[1] = vt_l*vt_l*maxwellian(n_l, v, 0.0, vt_l);
  }
  else {
    fout[0] = maxwellian(n_r, v, 0.0, vt_r);
    fout[1] = vt_r*vt_r*maxwellian(n_r, v, 0.0, vt_r);
  }
}

void 
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_es_sod_shock_ctx *app = ctx;
  double x = xn[0];
  // No initial flow (u = 0)
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalDistFuncIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_es_sod_shock_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double n_l = app->n0;
  double n_r = 0.125*n_l;
  double vt_l = app->vti;
  double vt_r = sqrt(0.1/0.125)*vt_l;
  if (fabs(x)<0.5*app->Lx) {
    fout[0] = maxwellian(n_l, v, 0.0, vt_l);
    fout[1] = vt_l*vt_l*maxwellian(n_l, v, 0.0, vt_l);
  }
  else {
    fout[0] = maxwellian(n_r, v, 0.0, vt_r);
    fout[1] = vt_r*vt_r*maxwellian(n_r, v, 0.0, vt_r);
  }
}

void 
evalFluidIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_es_sod_shock_ctx *app = ctx;
  double x = xn[0];
  // No initial flow (u = 0)
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_es_sod_shock_ctx *app = ctx;
  double x = xn[0];
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalExtEmFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_es_sod_shock_ctx *app = ctx;
  double x = xn[0];
  double B_x = app->B0;
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = B_x; fout[4] = 0.0; fout[5] = 0.0;
}


void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_es_sod_shock_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_es_sod_shock_ctx *app = ctx;
  fout[0] = app->nuIon;
}

struct pkpm_es_sod_shock_ctx
create_ctx(void)
{
  double Lx = 128.0; 
  int Nx = 256; 
  double dx = Lx/Nx;
  double t_end = 100.0;
  // initial dt guess so first step does not generate NaN
  double init_dt = ((Lx/Nx)/6.0)/(5.0);

  struct pkpm_es_sod_shock_ctx ctx = {
    .chargeElc = -1.0,
    .massElc = 1.0,
    .chargeIon = 1.0,
    .massIon = 1836.153,
    .Te_Ti = 4.0,
    .vte = 1.0,
    .vti = ctx.vte/sqrt(ctx.Te_Ti*ctx.massIon),
    .n0 = 1.0, 
    .B0 = 1.0, 
    .nuElc = 1.0e-4, 
    .nuIon = 1.0e-4/sqrt(ctx.massIon)*(ctx.Te_Ti*sqrt(ctx.Te_Ti)), 
    .Lx = Lx,
    .Nx = Nx, 
    .t_end = t_end, 
    .init_dt = init_dt, 
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  struct pkpm_es_sod_shock_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 48);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // electrons
  struct gkyl_pkpm_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -6.0 * ctx.vte},
    .upper = { 6.0 * ctx.vte}, 
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
  };

  // field
  struct gkyl_pkpm_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,

    .ext_em = evalExtEmFunc,
    .ext_em_ctx = &ctx,
  };

  // pkpm app
  struct gkyl_pkpm pkpm = {
    .name = "pkpm_periodic_es_sod_shock_p2",

    .cdim = 1, .vdim = 1,
    .lower = { -ctx.Lx },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .use_explicit_source = true, 

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&pkpm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.t_end;
  double dt = ctx.init_dt;

  // initialize simulation
  gkyl_pkpm_app_apply_ic(app, tcurr);
  
  gkyl_pkpm_app_write(app, tcurr, 0);
  gkyl_pkpm_app_calc_field_energy(app, tcurr);
  gkyl_pkpm_app_calc_integrated_L2_f(app, tcurr);
  gkyl_pkpm_app_calc_integrated_mom(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_pkpm_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);

    gkyl_pkpm_app_calc_field_energy(app, tcurr);
    gkyl_pkpm_app_calc_integrated_L2_f(app, tcurr);
    gkyl_pkpm_app_calc_integrated_mom(app, tcurr);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    step += 1;
  }

  gkyl_pkpm_app_write(app, tcurr, 1);
  gkyl_pkpm_app_write_field_energy(app);
  gkyl_pkpm_app_write_integrated_L2_f(app);
  gkyl_pkpm_app_write_integrated_mom(app);
  gkyl_pkpm_app_stat_write(app);

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
