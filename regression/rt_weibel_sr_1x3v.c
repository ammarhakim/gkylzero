#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct weibel_ctx {
  // parameters for plasma streams
  double nElc10, nElc20;
  double TElc10, TElc20;
  double uxElc10, uxElc20;
  double uyElc10, uyElc20;
  double uzElc10, uzElc20;

  // perturbation parameters
  double k0;
  double perturb;
};

static inline double
maxwelljuttner3D(double n, double px, double py, double pz, double ux, double uy, double uz, double T, double K_2)
{
  double gamma = 1.0/sqrt(1 - ux*ux - uy*uy - uz*uz);
  return n/(4*M_PI*T*K_2)*exp(-(gamma/T)*(sqrt(1 + px*px + py*py + pz*pz) - ux*px - uy*py - uz*pz));
}

void
evalDistFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct weibel_ctx *app = ctx;

  double nElc10 = app->nElc10, nElc20 = app->nElc20;
  double uxElc10 = app->uxElc10, uxElc20 = app->uxElc20;
  double uyElc10 = app->uyElc10, uyElc20 = app->uyElc20;
  double uzElc10 = app->uzElc10, uzElc20 = app->uzElc20;
  double TElc10 = app->TElc10, TElc20 = app->TElc20;
  
  double x = xn[0], px = xn[1], py = xn[2], pz = xn[3];

  // modified Bessel function of the second kind evaluated for T = mc^2 (K_2(1))
  //double K_2 = 1.6248388986351774828107073822838437146593935281628733843345054697;
  // modified Bessel function of the second kind evaluated for T = 0.1 mc^2 (K_2(10))
  //double K_2 = 0.0000215098170069327687306645644239671272492068461808732468335569;
  // modified Bessel function of the second kind evaluated for T = 0.04 mc^2 (K_2(25))
  double K_2 = 3.7467838080691090570137658745889511812329380156362352887017e-12;
  
  fout[0] = maxwelljuttner3D(nElc10, px, py, pz, uxElc10, uyElc10, uzElc10, TElc10, K_2) +
    maxwelljuttner3D(nElc20, px, py, pz, uxElc20, uyElc20, uzElc20, TElc20, K_2);
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
  double ud = 0.9;
  double k0 = 0.4;
  double massElc = 1.0;
  double TElc10 = 0.04;
  double TElc20 = 0.04;
  
  struct weibel_ctx ctx = {
    .nElc10 = 0.5,
    .nElc20 = 0.5,

    .uxElc10 = 0.0,
    .uyElc10 = ud,
    .uzElc10 = 0.0,

    .uxElc20 = 0.0,
    .uyElc20 = -ud,
    .uzElc20 = 0.0,

    .TElc10 = TElc10,
    .TElc20 = TElc20,

    .k0 = k0,
    .perturb = 1e-3,
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
  struct weibel_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 24);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 12);
  int VY = APP_ARGS_CHOOSE(app_args.vcells[1], 12);
  int VZ = APP_ARGS_CHOOSE(app_args.vcells[2], 12);

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_SR,
    .charge = -1.0, .mass = 1.0,
    .lower = { -8.0, -8.0, -8.0 },
    .upper = { 8.0, 8.0, 8.0 }, 
    .cells = { VX, VY, VZ },

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFunc,
      .ctx_func = &ctx,
    },

    .num_diag_moments = 2,
    .diag_moments = { "M0", "M1i" },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,

    .ctx = &ctx,
    .init = evalFieldFunc
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "weibel_sr_1x3v",

    .cdim = 1, .vdim = 3,
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
  double tcurr = 0.0, tend = 10.0;
  double dt = tend-tcurr;
  int nframe = 10;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

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
    write_data(&io_trig, app, tcurr);

    step += 1;
  }

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
