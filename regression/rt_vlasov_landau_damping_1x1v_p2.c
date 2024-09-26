#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct langmuir_ctx {
  double charge; // charge
  double mass; // mass
  double vt; // thermal velocity
  double Lx; // size of the box
  double k0; // wave number
  double perturb; // perturbation amplitude
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
  struct langmuir_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double k = app->k0, alpha = app->perturb;
  fout[0] = (1 + alpha*cos(k*x))*maxwellian(1.0, v, 0.0, 1.0);
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
evalNu(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct langmuir_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  fout[0] = 0.1;
}

struct langmuir_ctx
create_ctx(void)
{
  struct langmuir_ctx ctx = {
    .mass = 1.0,
    .charge = -1.0,
    .vt = 1.0,
    .Lx = M_PI/0.5,
    .k0 = 0.5,
    .perturb = 1e-4
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
  struct langmuir_ctx ctx = create_ctx(); // context for init functions

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = ctx.charge, .mass = ctx.mass,
    .lower = { -6.0*ctx.vt},
    .upper = { 6.0*ctx.vt}, 
    .cells = { 32 },

    .num_init = 1, 
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFunc,
      .ctx_func = &ctx,
    },

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNu,
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
    .name = "vlasov_landau_damping_1x1v_p2",

    .cdim = 1, .vdim = 1,
    .lower = { -ctx.Lx },
    .upper = { ctx.Lx },
    .cells = { 32 },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 1,
    .species = { elc },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
    },
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 20.;
  double dt = tend-tcurr;

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  
  gkyl_vlasov_app_write(app, tcurr, 0);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 0);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
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
  printf("Field RHS calc took %g secs\n", stat.field_rhs_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
