#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct pulsar_atm_ctx {
  double T; // temperature
  double vdrift; // drift velocity
  double J_0; // doubles as curl(B) and gravity
  double grav; // gravitational scale height = g/T
  double Lx; // size of box
  double perturb; // size of perturbation for exciting instabilities
  double elc_pmax; // electron momentum space extents
  int Nv_elc; // number of electron momentum grid points
  double pos_pmax; // positron momentum space extents
  int Nv_pos; // number of positron momentum grid points
};

static inline double sq(double x) { return x*x; }

void
evalElcDistFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_atm_ctx *app = ctx;
  double x = xn[0], p = xn[1];
  double T = app->T, vdrift = app->vdrift;
  double grav = app->grav;

  // modified Bessel function of the second kind evaluated for T = mc^2 (K_2(1))
  //double K_2 = 1.6248388986351774828107073822838437146593935281628733843345054697;
  // modified Bessel function of the second kind evaluated for T = 0.1 mc^2 (K_2(10))
  double K_2 = 0.0000215098170069327687306645644239671272492068461808732468335569;
  // modified Bessel function of the second kind evaluated for T = 0.04 mc^2 (K_2(25))
  //double K_2 = 3.7467838080691090570137658745889511812329380156362352887017e-12;
  // Lorentz factor for drift velocity
  double gamma = 1.0/sqrt(1 - vdrift*vdrift);

  double n = 10000.0*exp(-grav*x); // scale height = grav
  double nElc = 0.0;
  /* if (n > 10.0) */
  /*   nElc = 1.0; // extra charge density to make E = 0 in atmosphere */
  double mc2_T = 1.0/T;

  double fv = (n+nElc)/K_2*exp(-mc2_T*gamma*(sqrt(1 + p*p) - vdrift*p));
  fout[0] = fv;
}

void
evalDistFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_atm_ctx *app = ctx;
  double x = xn[0], p = xn[1];
  double T = app->T, vdrift = app->vdrift;
  double grav = app->grav;

  // modified Bessel function of the second kind evaluated for T = mc^2 (K_2(1))
  //double K_2 = 1.6248388986351774828107073822838437146593935281628733843345054697;
  // modified Bessel function of the second kind evaluated for T = 0.1 mc^2 (K_2(10))
  double K_2 = 0.0000215098170069327687306645644239671272492068461808732468335569;
  // modified Bessel function of the second kind evaluated for T = 0.04 mc^2 (K_2(25))
  //double K_2 = 3.7467838080691090570137658745889511812329380156362352887017e-12;
  // Lorentz factor for drift velocity
  double gamma = 1.0/sqrt(1 - vdrift*vdrift);

  double n = 10000.0*exp(-grav*x); // scale height = grav
  double mc2_T = 1.0/T;

  double fv = n/K_2*exp(-mc2_T*gamma*(sqrt(1 + p*p) - vdrift*p));
  fout[0] = fv;
}

void
evalAccelFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_atm_ctx *app = ctx;
  double x = xn[0];
  double grav = app->grav;
  double T = app->T;
  fout[0] = -grav*T; fout[1] = 0.0, fout[2] = 0.0;
}

void
evalCurrentFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_atm_ctx *app = ctx;
  double x = xn[0];
  double J_0 = app->J_0;
  double Lx = app->Lx;
  double grav = app->grav;
  // profile adjusts from J_0 -> 1 over a scale height halfway through the box
  //double curr = J_0/2.0*tanh(grav*(x - Lx/2.0)) + (1.0 - J_0);
  fout[0] = J_0; fout[1] = 0.0, fout[2] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pulsar_atm_ctx *app = ctx;
  double x = xn[0];
  double grav = app->grav;

  pcg64_random_t rng = gkyl_pcg64_init(0); // RNG for use in IC
  double Lx = app->Lx;
  double kx = 2.0*M_PI/Lx;
  double perturb = app->perturb;

  // linear profile in E_x corresponds to potential drop
  double E_x = x;
  /* // add perturbations to Ex to excite instabilities */
  /* for (int i=0; i<64; ++i)  */
  /*   E_x += perturb*gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng)); */
  
  double n = 10000.0*exp(-grav*x); // scale height = grav
  /* // Ex equals 0 inside atmosphere where n = n_GJ */
  /* if (n > 10.0) */
  /*   E_x = 0.0; */

  fout[0] = E_x; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

struct pulsar_atm_ctx
create_default_ctx(void)
{
  return (struct pulsar_atm_ctx) {
    .T = 0.1,
    .vdrift = 0.0,
    .grav = 1.0,
    .J_0 = 0.5,
    .Lx = 100.0,
    .perturb = 1.0e-4,
    .elc_pmax = 512.0,
    .Nv_elc = 8192,
    .pos_pmax = 512.0,
    .Nv_pos = 8192,
  };
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 1024);

  struct pulsar_atm_ctx ctx;
  
  ctx = create_default_ctx(); // context for init functions

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = -1.0,
    .mass = 1.0,
    .model_id = GKYL_MODEL_SR,
    .lower = { -ctx.elc_pmax },
    .upper = { ctx.elc_pmax }, 
    .cells = { ctx.Nv_elc },

    .ctx = &ctx,
    .init = evalElcDistFunc,

    .accel = evalAccelFunc,
    .accel_ctx = &ctx,

    .bcx = { GKYL_SPECIES_FIXED_FUNC, GKYL_SPECIES_COPY },

    .num_diag_moments = 2,
    .diag_moments = { "M0", "M1i" },
  };

  // positrons
  struct gkyl_vlasov_species pos = {
    .name = "pos",
    .model_id = GKYL_MODEL_SR,
    .charge = 1.0,
    .mass = 1.0,
    .lower = { -ctx.pos_pmax },
    .upper = { ctx.pos_pmax }, 
    .cells = { ctx.Nv_pos },

    .ctx = &ctx,
    .init = evalDistFunc,

    .accel = evalAccelFunc,
    .accel_ctx = &ctx,

    .bcx = { GKYL_SPECIES_FIXED_FUNC, GKYL_SPECIES_COPY },

    .num_diag_moments = 2,
    .diag_moments = { "M0", "M1i" },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,

    .app_current = evalCurrentFunc,
    .app_current_ctx = &ctx,

  };

  // VM app
  struct gkyl_vm vm = {
    .name = "pulsar_atm_sr_sub_GJ",
    
    .cdim = 1, .vdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 2,
    .basis_type = app_args.basis_type,
    .cfl_frac = 0.5, //can increase stability
    
    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, pos },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };
  
  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 500.0;
  double dt = tend-tcurr;
  int nframe = 100;
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
      fprintf(stderr, "** Update method failed! Aborting simulation ....\n");
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

  // simulation complete, free objects
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
