#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm.h>
#include <rt_arg_parse.h>

struct langmuir_ctx {
  double charge; // charge
  double mass; // mass
  double vt; // thermal velocity
  double perturb; // perturbation amplitude
  double k0; // wave number
  double Lx; // Domain size (x-direction).
  int Nx; // Cell count (x-direction).
  double t_end; // Final simulation time.
  double init_dt; // Initial time step guess so first step does not generate NaN
};

static inline double sq(double x) { return x*x; }

static inline double
maxwellian(double n, double v, double vth)
{
  double v2 = v*v;
  return n/sqrt(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

void
evalDistFuncElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct langmuir_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double k = app->k0, alpha = app->perturb;
  fout[0] = (1 + alpha*cos(k*x))*maxwellian(1.0, v, 1.0);
  // T_perp/m = 1.0
  fout[1] = (1 + alpha*cos(k*x))*maxwellian(1.0, v, 1.0);
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  // no initial flow
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct langmuir_ctx *app = ctx;

  double k = app->k0, alpha = app->perturb;
  double x = xn[0];
  double E_x = -alpha*sin(k*x)/k;
  
  fout[0] = E_x; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalExtEmFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_sheath_ctx *app = ctx;
  double x = xn[0];
  double B_x = 1.0;
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = B_x; fout[4] = 0.0; fout[5] = 0.0;
}

void
evalNu(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct langmuir_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  fout[0] = 1.0e-4;
}

struct langmuir_ctx
create_ctx(void)
{
  double k0 = 0.5;
  double Lx = M_PI/k0; 
  int Nx = 32; 
  double dx = Lx/Nx;
  double t_end = 20.0;
  // initial dt guess so first step does not generate NaN
  double init_dt = (Lx/Nx)/(3.0);

  struct langmuir_ctx ctx = {
    .mass = 1.0,
    .charge = -1.0,
    .vt = 1.0,
    .perturb = 1.e-1, 
    .k0 = k0,
    .Lx = Lx,
    .Nx = Nx, 
    .t_end = t_end, 
    .init_dt = init_dt, 
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

  struct langmuir_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 32);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // electrons
  struct gkyl_pkpm_species elc = {
    .name = "elc",
    .charge = ctx.charge, .mass = ctx.mass,
    .lower = { -6.0*ctx.vt},
    .upper = { 6.0*ctx.vt}, 
    .cells = { NV },

    .ctx_dist = &ctx,
    .ctx_fluid = &ctx,
    .init_dist = evalDistFuncElc,
    .init_fluid = evalFluidElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNu,
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
    .name = "pkpm_landau_damping_p1",

    .cdim = 1, .vdim = 1,
    .lower = { -ctx.Lx },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .use_explicit_source = true, 

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 1,
    .species = { elc },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&pkpm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.t_end;
  double dt = ctx.init_dt;
  int nframe = 1;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };


  // initialize simulation
  gkyl_pkpm_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
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
    write_data(&io_trig, app, tcurr);

    step += 1;
  }

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
