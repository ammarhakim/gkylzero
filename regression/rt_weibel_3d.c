#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct weibel_ctx {
  // parameters for plasma streams
  double nElc10, nElc20;
  double vthElc10, vthElc20;
  double uxElc10, uxElc20;
  double uyElc10, uyElc20;

  // perturbation parameters
  double k0;
  double perturb;
};

static inline double
maxwellian2D(double n, double vx, double vy, double ux, double uy, double vth)
{
  double v2 = (vx-ux)*(vx-ux) + (vy-uy)*(vy-uy);
  return n/(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

void
evalDistFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct weibel_ctx *app = ctx;

  double nElc10 = app->nElc10, nElc20 = app->nElc20;
  double uxElc10 = app->uxElc10, uxElc20 = app->uxElc20;
  double uyElc10 = app->uyElc10, uyElc20 = app->uyElc20;
  double vthElc10 = app->vthElc10, vthElc20 = app->vthElc20;
  
  double x = xn[0], vx = xn[1], vy = xn[2];
  
  fout[0] = maxwellian2D(nElc10, vx, vy, uxElc10, uyElc10, vthElc10) +
    maxwellian2D(nElc20, vx, vy, uxElc20, uyElc20, vthElc20);
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct weibel_ctx *app = ctx;

  double perturb = app->perturb, k0 = app->k0;
  double x = xn[0];
  double B_z = perturb*sin(k0*x);
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = B_z;
  fout[6] = 0.0; fout[7] = 0.0;
}

struct weibel_ctx
create_ctx(void)
{
  double ud = 0.3;
  double k0 = 1.0;
  double massElc = 1.0;
  double TElc10 = 0.01;
  double TElc20 = 0.01;
  double vthElc10 = sqrt(TElc10/massElc);
  double vthElc20 = sqrt(TElc20/massElc);  
  
  struct weibel_ctx ctx = {
    .nElc10 = 0.5,
    .nElc20 = 0.5,
    .uxElc10 = 0.0,
    .uyElc10 = 0.3,
    .uxElc20 = 0.0,
    .uyElc20 = -0.3,

    .vthElc10 = vthElc10,
    .vthElc20 = vthElc20,

    .k0 = 0.4,
    .perturb = 1e-3,
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct weibel_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 24);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 12);
  int VY = APP_ARGS_CHOOSE(app_args.vcells[1], 12);

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = -1.0, .mass = 1.0,
    .lower = { -1.0, -1.0 },
    .upper = { 1.0, 1.0 }, 
    .cells = { VX, VY },

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFunc,
      .ctx_func = &ctx,
    },

    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2" },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,

    .ctx = &ctx,
    .init = evalFieldFunc
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "weibel_3d",

    .cdim = 1, .vdim = 2,
    .lower = { 0.0 },
    .upper = { 2*M_PI/ctx.k0 },
    .cells = { NX },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 1,
    .species = { elc },
    .field = field,

    .use_gpu = app_args.use_gpu
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 80.0;
  double dt = tend-tcurr;

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  
  gkyl_vlasov_app_write(app, tcurr, 0);
  gkyl_vlasov_app_calc_mom(app);
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
    step += 1;
  }

  gkyl_vlasov_app_write(app, tcurr, 1);
  gkyl_vlasov_app_calc_mom(app);
  gkyl_vlasov_app_write_mom(app, tcurr, 0);
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
  printf("Field RHS calc took %g secs\n", stat.field_rhs_tm);
  printf("Current evaluation and accumulate took %g secs\n", stat.current_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
