#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct lbo_cross_ctx {
  double q1, q2; // charge
  double m1, m2; // mass
  double vt1, vt2; // thermal velocity
  double n1, n2;
  double u1, u2;
  double p1, p2;
  double K1, K2;
};

static inline double sq(double x) { return x*x; }

static inline double
maxwellian(double n, double vx, double vy, double ux, double vt)
{
  double v2 = (vx-ux)*(vx-ux) + vy*vy;
  return n/(2*M_PI*vt*vt)*exp(-v2/(2*vt*vt));
}

void
evalDistFunc1(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct lbo_cross_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2];

  fout[0] = maxwellian(app->n1, vx, vy, app->u1, app->vt1);
}

void
evalDistFunc2(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct lbo_cross_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2];

  fout[0] = maxwellian(app->n2, vx, vy, app->u2, app->vt2);;
}

void
evalNu1(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct lbo_cross_ctx *app = ctx;

  fout[0] = 1.0/0.01;
}

void
evalNu2(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct lbo_cross_ctx *app = ctx;

  fout[0] = sqrt(0.5/0.05)/0.01;
}

struct lbo_cross_ctx
create_ctx(void)
{
  double m1 = 1.0;
  double m2 = 0.05;
  double n1 = 1.0;
  double n2 = 1.0;
  double p1 = 1.0;
  double p2 = 0.50;
  double vt1 = sqrt(p1/(n1*m1));
  double vt2 = sqrt(p2/(n2*m2));
  struct lbo_cross_ctx ctx = {
    .m1 = m1,
    .m2 = m2,
    .n1 = n1,
    .n2 = n2,
    .u1 = 0.10,
    .u2 = 2.50,
    .q1 = 0.0,
    .q2 = 0.0,
    .p1 = p1,
    .p2 = p2,
    .vt1 = vt1,
    .vt2 = vt2,
    .K1 = 0.01/n1,
    .K2 = 0.01/n2,
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
  struct lbo_cross_ctx ctx = create_ctx(); // context for init functions

  struct gkyl_vlasov_species neut1 = {
    .name = "neut1",
    .charge = ctx.q1, .mass = ctx.m1,
    .lower = { -6.0*ctx.vt1, -6.0*ctx.vt1 },
    .upper = {  6.0*ctx.vt1,  6.0*ctx.vt1 }, 
    .cells = { 32, 32 },

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFunc1,
      .ctx_func = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalNu1,
      .num_cross_collisions = 1,
      .collide_with = { "neut2" },
    },
    
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2" },
  };

  struct gkyl_vlasov_species neut2 = {
    .name = "neut2",
    .charge = ctx.q2, .mass = ctx.m2,
    .lower = { -6.0*ctx.vt2, -6.0*ctx.vt2 },
    .upper = {  6.0*ctx.vt2,  6.0*ctx.vt2 }, 
    .cells = { 32, 32 },

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFunc2,
      .ctx_func = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalNu2,
      .num_cross_collisions = 1,
      .collide_with = { "neut1" },
    },
    
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2" },
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "lbo_vlasov_cross_1x2v_p2",

    .cdim = 1, .vdim = 2,
    .lower = { 0.0 },
    .upper = { 1.0 },
    .cells = { 16 },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { neut1, neut2 },
    .skip_field = true,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);
  // start, end and initial time-step
  double tcurr = 0.0, tend = 0.0025;
  double dt = tend-tcurr;

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  
  gkyl_vlasov_app_write(app, tcurr, 0);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 0);
  gkyl_vlasov_app_calc_integrated_mom(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    gkyl_vlasov_app_calc_integrated_mom(app, tcurr);
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
  gkyl_vlasov_app_write_integrated_mom(app);
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
