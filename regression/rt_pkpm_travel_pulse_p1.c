#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct pkpm_sod_shock_ctx {
  double charge; // charge
  double mass; // mass
  double vt; // thermal velocity
  double Lx; // size of the box
};

static inline double sq(double x) { return x*x; }

static inline double
maxwellian(double n, double v, double u, double vth)
{
  double v2 = (v - u)*(v - u);
  return n/sqrt(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

void
evalDistFunc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_sod_shock_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double n = 1.0 + 0.2*sin(M_PI*x);
  double p = 1.0;
  fout[0] = maxwellian(n, v, 0.0, sqrt(p/n));
}

void 
evalFluidFunc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_sod_shock_ctx *app = ctx;
  double x = xn[0];
  double n = 1.0 + 0.2*sin(M_PI*x);
  double p = 1.0, u = 1.0;
  double gas_gamma = 5.0/3.0;
  
  fout[0] = n*u;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = p;
}

void
evalNu(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_sod_shock_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  fout[0] = 100.0;
}

struct pkpm_sod_shock_ctx
create_ctx(void)
{
  struct pkpm_sod_shock_ctx ctx = {
    .mass = 1.0,
    .charge = 0.0,
    .vt = 1.0,
    .Lx = 2.0,
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 16);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 16);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct pkpm_sod_shock_ctx ctx = create_ctx(); // context for init functions

  // PKPM fluid                                                                                      
  struct gkyl_vlasov_fluid_species fluid = {
    .name = "fluid",
    .num_eqn = 4,
    .pkpm_species = "neut",
    .ctx = &ctx,
    .init = evalFluidFunc,
  };  

  // Neutrals
  struct gkyl_vlasov_species neut = {
    .name = "neut",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid",
    .charge = ctx.charge, .mass = ctx.mass,
    .lower = { -6.0*ctx.vt},
    .upper = { 6.0*ctx.vt}, 
    .cells = { NV },

    .ctx = &ctx,
    .init = evalDistFunc,

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNu,
    },

    .num_diag_moments = 0,
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "pkpm_travel_pulse_p1",

    .cdim = 1, .vdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 1,
    .species = { neut },
    .num_fluid_species = 1,
    .fluid_species = { fluid },
    .skip_field = true,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 2.0;
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
