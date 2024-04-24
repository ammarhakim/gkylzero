#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct weibel_ctx {
  // parameters for plasma streams
  double nElc10, nElc20;
  double vthElc10, vthElc20;
  double uxElc10, uxElc20;
  double uyElc10, uyElc20;

  // perturbation parameters
  double kx, ky;
  double alpha; // ratio of E_y/E_x
  double perturb_n;
  bool use_gpu;
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
  double kx = app->kx, ky = app->ky, perturb_n = app->perturb_n;  
  
  double x = xn[0], y = xn[1], vx = xn[2], vy = xn[3];
  
  double fv = maxwellian2D(nElc10, vx, vy, uxElc10, uyElc10, vthElc10) +
    maxwellian2D(nElc20, vx, vy, uxElc20, uyElc20, vthElc20);
    
  fout[0] = (1.0+app->perturb_n*cos(kx*x+ky*y))*fv;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct weibel_ctx *app = ctx;

  double perturb_n = app->perturb_n, alpha = app->alpha;
  double kx = app->kx, ky = app->ky;
  
  double x = xn[0], y = xn[1];
  
  double E_x = -perturb_n*sin(kx*x+ky*y)/(kx+ky*alpha);
  double E_y = alpha*E_x;
  double B_z = kx*E_y-ky*E_x;
  
  fout[0] = E_x; fout[1] = E_y, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = B_z;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalNu(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct weibel_ctx *app = ctx;
  fout[0] = 1.0e-4;
}

struct weibel_ctx
create_ctx(void)
{
  double ud = 0.3;
  double k0 = 1.0, theta = 45.0/180.0*M_PI;
  double kx = k0*cos(theta), ky = k0*sin(theta);

  double massElc = 1.0, R = 0.333333333333333;
  double TElc10 = massElc*R*ud*R*ud;
  double TElc20 = TElc10;
  double vthElc10 = sqrt(TElc10/massElc);
  double vthElc20 = sqrt(TElc20/massElc);  
  
  struct weibel_ctx ctx = {
    .nElc10 = 0.5,
    .nElc20 = 0.5,
    .uxElc10 = 0.0,
    .uxElc20 = 0.0,
    .uyElc10 = ud,
    .uyElc20 = -ud,
    .vthElc10 = vthElc10,
    .vthElc20 = vthElc20,

    .kx = kx, .ky = ky,
    .alpha = 1.18281106421231,
    .perturb_n = 1e-8,
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 8);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 8);
  
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 16);
  int VY = APP_ARGS_CHOOSE(app_args.vcells[1], 16);  

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = -1.0, .mass = 1.0,
    .lower = { -0.9, -0.9 },
    .upper = { 0.9, 0.9 }, 
    .cells = { VX, VY },

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFunc,
      .ctx_func = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalNu,
    },
    
    .num_diag_moments = 2,
    .diag_moments = { "M0", "M1i" },
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
    .name = "weibel_lbo_4d",

    .cdim = 2, .vdim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { 2*M_PI/ctx.kx, 2*M_PI/ctx.ky },
    .cells = { NX, NY },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1},

    .num_species = 1,
    .species = { elc },
    .field = field,
 
    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 5.0;
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
  printf("Species collisions took %g secs\n", stat.species_coll_mom_tm);  
  printf("Species collisions took %g secs\n", stat.species_coll_tm);
  printf("Field RHS calc took %g secs\n", stat.field_rhs_tm);
  printf("Current evaluation and accumulate took %g secs\n", stat.current_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
