#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

#define sq(x) ((x)*(x))

struct moment_ctx {
  // following will be set to 1 to define unit normalizations
  double mu0;
  double vA0;
  double rho0;
  double l0;
  double n0;
  double mi;

  // problem-specific mhd parameters
  double gas_gamma;
  double Bz0__B0;
  double T0_mhd;
  // set one and only one of MA0 and gravity to nonzero
  double MA0;  // vTheta0/vA0; set to nonzero to set a Er x Bz drift
  double gravity; // allows an equilibrium to develop

  // non-mhd parameters
  double lightSpeed__vA0;
  double mi__me;
  double Ti0__Te0;
  double di0__l0;

  // domain
  double r_inn;
  double r_out;
};

struct moment_ctx
moment_ctx(void)
{
  return (struct moment_ctx) {
    .mu0 = 1,
    .vA0 = 1,
    .rho0 = 1,
    .l0 = 1,
    .n0 = 1,
    .mi = 1,

    .gas_gamma = 1.4,
    .Bz0__B0 = 1,
    .T0_mhd = 0.002,
    .MA0 = 0,
    .gravity = 0.0125,

    .lightSpeed__vA0 = 10,
    .mi__me = 25,
    .Ti0__Te0 = 1,
    .di0__l0 = 1,

    .r_inn = 0.45,
    .r_out = 1.45,
  };
}

// map (r,theta) -> (x,y)
void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double r = xc[0], th = xc[1];
  xp[0] = r*cos(th); xp[1] = r*sin(th);
}

// initial conditions
double calc_vt(const double r, const struct moment_ctx *ctx)
{
  double vt0 = ctx->vA0 * ctx->MA0;
  double vt = vt0 * sin((r - ctx->r_inn) * M_PI / (ctx->r_out - ctx->r_inn));
  return vt;
}

void evalElcInit(
    double t,
    const double* GKYL_RESTRICT xn,
    double* GKYL_RESTRICT fout,
    void *_ctx)
{
  struct moment_ctx *ctx = _ctx;
  double gamma = ctx->gas_gamma;
  double rho0 = ctx->rho0;
  double mi__me = ctx->mi__me;
  double me = ctx->mi / mi__me;
  double Ti0__Te0 = ctx->Ti0__Te0;
  double T0_mhd = ctx->T0_mhd;

  double r = xn[0]; // radial coordinate

  double rhoe = rho0 / (1 + mi__me);;
  double vre = 0;
  double vte = calc_vt(r, ctx);
  double vze = 0;
  double ne = rhoe / me;
  double Te = T0_mhd / (1 + Ti0__Te0);
  double pe = ne * Te;
  double ere = pe / (gamma - 1) + 0.5 * rhoe * (sq(vre) + sq(vte) + sq(vze));

  fout[0] = rhoe;
  fout[1] = rhoe * vre;
  fout[2] = rhoe * vte;
  fout[3] = rhoe * vze;
  fout[4] = ere;
}

void evalIonInit(
    double t,
    const double* GKYL_RESTRICT xn,
    double* GKYL_RESTRICT fout,
    void *_ctx)
{
  struct moment_ctx *ctx = _ctx;
  double gamma = ctx->gas_gamma;
  double rho0 = ctx->rho0;
  double mi__me = ctx->mi__me;
  double mi = ctx->mi;
  double Ti0__Te0 = ctx->Ti0__Te0;
  double T0_mhd = ctx->T0_mhd;

  double r = xn[0]; // radial coordinate

  double rhoi = rho0 * mi__me / (1 + mi__me);;
  double vri = 0;
  double vti = calc_vt(r, ctx);
  double vzi = 0;
  double ni = rhoi / mi;
  double Ti = T0_mhd * Ti0__Te0 / (1 + Ti0__Te0);
  double pi = ni * Ti;
  double eri = pi / (gamma - 1) + 0.5 * rhoi * (sq(vri) + sq(vti) + sq(vzi));

  fout[0] = rhoi;
  fout[1] = rhoi * vri;
  fout[2] = rhoi * vti;
  fout[3] = rhoi * vzi;
  fout[4] = eri;
}

void
evalFieldInit(
    double t,
    const double* GKYL_RESTRICT xn,
    double* GKYL_RESTRICT fout,
    void *_ctx)
{
  struct moment_ctx *ctx = _ctx;
  double B0 = ctx->vA0 * sqrt(ctx->mu0 * ctx->rho0);
  double Bz0 = B0 * ctx->Bz0__B0;

  double r = xn[0]; // radial coordinate

  double vt = calc_vt(r, ctx);
  double Er =  -vt * Bz0;
  double Et = 0;
  double Ez = 0;
  double Br = 0;
  double Bt = 0;
  double Bz = 0;

  // electric field
  fout[0] = Er; fout[1] = Et; fout[2] = Ez;
  // magnetic field
  fout[3] = Br; fout[4] = Bt; fout[5] = Bz;
  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
}

void gravity_func(double t, const double *xn, double *fout, void *_ctx)
{
  struct moment_ctx *ctx = _ctx;
  fout[0] = ctx->gravity;
}

void
write_data(
    struct gkyl_tm_trigger *iot,
    const gkyl_moment_app *app,
    double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr))
  {
    printf("Writing frame %d at t=%g\n", iot->curr-1, tcurr);
    gkyl_moment_app_write(app, tcurr, iot->curr-1);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NR = APP_ARGS_CHOOSE(app_args.xcells[0], 64);
  int NT = APP_ARGS_CHOOSE(app_args.xcells[1], 3);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct moment_ctx ctx = moment_ctx();
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(ctx.gas_gamma);

  double lightSpeed = ctx.vA0 * ctx.lightSpeed__vA0;
  double epsilon0 = 1 / sq(lightSpeed) / ctx.mu0;
  double me = ctx.mi / ctx.mi__me;
  double di0 = ctx.l0 * ctx.di0__l0;
  double wpi0 = ctx.lightSpeed__vA0 / di0;
  double qi = wpi0 * sqrt(epsilon0 * ctx.mi / ctx.n0);
  double qe = -qi;

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = qe,
    .mass = me,
    .equation = euler,
    .evolve = true,
    .ctx = &ctx,
    .init = evalElcInit,
    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
    .bcy = { GKYL_SPECIES_WEDGE, GKYL_SPECIES_WEDGE },
    .app_accel_func = gravity_func,
  };

  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = qi,
    .mass = ctx.mi,
    .equation = euler,
    .evolve = true,
    .ctx = &ctx,
    .init = evalIonInit,
    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
    .bcy = { GKYL_SPECIES_WEDGE, GKYL_SPECIES_WEDGE },
    .app_accel_func = gravity_func,
  };

  double theta = 0.01;
  struct gkyl_moment app_inp = {
    .name = "5m_radial",

    .ndim = 2,
    .lower = { ctx.r_inn, -theta / 2 },
    .upper = { ctx.r_out, theta / 2 },
    .cells = { NR, NT },
    .mapc2p = mapc2p,

    .cfl_frac = 0.9,

    .num_species = 2,
    .species = { elc, ion },
    .field = {
      .epsilon0 = epsilon0,
      .mu0 = ctx.mu0,
      .mag_error_speed_fact = 1.0,
      .evolve = true,
      .ctx = &ctx,
      .init = evalFieldInit,
      .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
      .bcy = { GKYL_FIELD_WEDGE, GKYL_FIELD_WEDGE },
    }
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tA0 = ctx.l0 / ctx.vA0;
  double tcurr = 0.0, tend = 10 * tA0;
  int nframe = 10;
  
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  printf("%45s = %g\n", "tA0", tA0);
  printf("%45s = %g\n", "tend", tend);
  printf("%45s = %g\n", "ouput interval dt", tend/nframe);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    bool do_print = step % 1000 == 0;

    if (do_print)
      printf("Taking time-step %6ld at t = %8g ...", step, tcurr);

    struct gkyl_update_status status = gkyl_moment_update(app, dt);

    if (do_print)
      printf(" dt = %8g, next frame = %d\n", status.dt_actual, io_trig.curr);

    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  // simulation complete, free resources
  gkyl_wv_eqn_release(euler);
  gkyl_moment_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);

  return 0;
}

#undef sq
