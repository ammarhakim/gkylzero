#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

#define sq(x) ((x) * (x))

/*********************************************/
/* CONTEXT PARAMETERS                        */
/* Extended from MHD normalization           */
/*********************************************/

struct moment_ctx {
  // following will be set to 1 to define unit normalizations
  double mu0;  // vaccum permeability
  double l0;   // reference length
  double vA0;  // reference Alfven speed
  double rho0; // reference mass density
  double mi;   // ion mass

  // mhd parameters
  double gas_gamma; // gas gamma
  double Bz0__B0;   // background Bz / reference B0
  double T0;        // background temperature Te0 + Ti0
  double vt0__vA0;  // vTheta0/vA0 to set radial electric field
  double gravity;   // radial gravity FIXME unit

  // non-mhd parameters
  double lightSpeed__vA0; // lightSpeed / vA0
  double mi__me;          // mass ratio mi / me
  double Ti0__Te0;        // temperature ratio Ti0 / Te0
  double di0__l0;         // ion inertial length / l0

  // domain limits needed for ramping initial velocity when vt0=/=0
  double r_inn; // inner radial; in units of l0, which is 1
  double r_out; // outer radial; in units of l0, which is 1

  // other parameters collected here for convenience
  double tend; // time at end of simulation; in units of tA0=l0/vA0=1
  int nframe;  // number of output frames
  int NR;      // radial cell number
  double cfl;
};

struct moment_ctx moment_ctx(void) {
  return (struct moment_ctx){
      .mu0 = 1,
      .vA0 = 1,
      .rho0 = 1,
      .l0 = 1,
      .mi = 1,

      .gas_gamma = 1.4,
      .Bz0__B0 = 1,
      .T0 = 0.002,
      .vt0__vA0 = 0,
      .gravity = 0.0025,

      .lightSpeed__vA0 = 10,
      .mi__me = 25,
      .Ti0__Te0 = 1,
      .di0__l0 = 1e-2,

      .r_inn = 0.45,
      .r_out = 1.45,

      .tend = 10,
      .nframe = 10,
      .NR = 64,
      .cfl = 0.9,
  };
}

/*********************************************/
/* INITIAL CONDITIONS                        */
/*********************************************/

double calc_vt(const double r, const struct moment_ctx *ctx) {
  double vt0 = ctx->vA0 * ctx->vt0__vA0;
  double vt = vt0 * sin((r - ctx->r_inn) * M_PI / (ctx->r_out - ctx->r_inn));
  return vt;
}

void init_elc(double t, const double *GKYL_RESTRICT xn,
              double *GKYL_RESTRICT fout, void *_ctx) {
  struct moment_ctx *ctx = _ctx;
  double gamma = ctx->gas_gamma;
  double rho0 = ctx->rho0;
  double mi__me = ctx->mi__me;
  double me = ctx->mi / mi__me;
  double Ti0__Te0 = ctx->Ti0__Te0;
  double T0 = ctx->T0;

  double r = xn[0];

  double rhoe = rho0 / (1 + mi__me);
  double vre = 0;
  double vte = calc_vt(r, ctx);
  double vze = 0;
  double ne = rhoe / me;
  double Te = T0 / (1 + Ti0__Te0);
  double pe = ne * Te;
  double ere = pe / (gamma - 1) + 0.5 * rhoe * (sq(vre) + sq(vte) + sq(vze));

  fout[0] = rhoe;
  fout[1] = rhoe * vre;
  fout[2] = rhoe * vte;
  fout[3] = rhoe * vze;
  fout[4] = ere;
}

void init_ion(double t, const double *GKYL_RESTRICT xn,
              double *GKYL_RESTRICT fout, void *_ctx) {
  struct moment_ctx *ctx = _ctx;
  double gamma = ctx->gas_gamma;
  double rho0 = ctx->rho0;
  double mi__me = ctx->mi__me;
  double mi = ctx->mi;
  double Ti0__Te0 = ctx->Ti0__Te0;
  double T0 = ctx->T0;

  double r = xn[0];

  double rhoi = rho0 / (1 + 1 / mi__me);
  double vri = 0;
  double vti = calc_vt(r, ctx);
  double vzi = 0;
  double ni = rhoi / mi;
  double Ti = T0 / (1 + 1 / Ti0__Te0);
  double pi = ni * Ti;
  double eri = pi / (gamma - 1) + 0.5 * rhoi * (sq(vri) + sq(vti) + sq(vzi));

  fout[0] = rhoi;
  fout[1] = rhoi * vri;
  fout[2] = rhoi * vti;
  fout[3] = rhoi * vzi;
  fout[4] = eri;
}

void init_field(double t, const double *GKYL_RESTRICT xn,
                double *GKYL_RESTRICT fout, void *_ctx) {
  struct moment_ctx *ctx = _ctx;
  double B0 = ctx->vA0 * sqrt(ctx->mu0 * ctx->rho0);
  double Bz0 = B0 * ctx->Bz0__B0;

  double r = xn[0];

  double vt = calc_vt(r, ctx);
  double Er = -vt * Bz0;
  double Et = 0;
  double Ez = 0;
  double Br = 0;
  double Bt = 0;
  double Bz = 0;

  fout[0] = Er;
  fout[1] = Et;
  fout[2] = Ez;
  fout[3] = Br;
  fout[4] = Bt;
  fout[5] = Bz;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

/*********************************************/
/* MISC.                                     */
/*********************************************/

void accel_func(double t, const double *xn, double *accel, void *_ctx) {
  struct moment_ctx *ctx = _ctx;
  accel[0] = ctx->gravity;
  accel[1] = 0;
  accel[2] = 0;
}

/*********************************************/
/* GRID MAPPING                              */
/*********************************************/

// map (r,theta) -> (x,y)
void mapc2p(double t, const double *xc, double *GKYL_RESTRICT xp, void *ctx) {
  double r = xc[0], theta = xc[1];
  xp[0] = r * cos(theta);
  xp[1] = r * sin(theta);
}

/*********************************************/
/* IO                                        */
/*********************************************/

void write_data(struct gkyl_tm_trigger *iot, const gkyl_moment_app *app,
                double tcurr) {
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    printf("Writing frame %d at t=%g\n", iot->curr - 1, tcurr);
    gkyl_moment_app_write(app, tcurr, iot->curr - 1);
  }
}

/*********************************************/
/* MAIN PROGRAM                              */
/*********************************************/

int main(int argc, char **argv) {
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

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
  double wpi0 = lightSpeed / di0;
  double rhoi0 = ctx.rho0 / (1 + 1 / ctx.mi__me);
  double qi = ctx.mi * sqrt(epsilon0 / rhoi0) * wpi0;
  double qe = -qi;

  struct gkyl_moment_species elc = {
      .name = "elc",
      .charge = qe,
      .mass = me,
      .equation = euler,
      .evolve = true,
      .ctx = &ctx,
      .init = init_elc,
      .bcx = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT},
      .bcy = {GKYL_SPECIES_WEDGE, GKYL_SPECIES_WEDGE},
      .app_accel_func = accel_func,
  };

  struct gkyl_moment_species ion = {
      .name = "ion",
      .charge = qi,
      .mass = ctx.mi,
      .equation = euler,
      .evolve = true,
      .ctx = &ctx,
      .init = init_ion,
      .bcx = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT},
      .bcy = {GKYL_SPECIES_WEDGE, GKYL_SPECIES_WEDGE},
      .app_accel_func = accel_func,
  };

  int NR = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.NR);
  int NT = APP_ARGS_CHOOSE(app_args.xcells[1], 3);
  double theta = 0.01; // some finite value not too small

  struct gkyl_moment app_inp = {
      .name = "5m_radial",

      .ndim = 2,
      .lower = {ctx.r_inn, -theta / 2},
      .upper = {ctx.r_out, theta / 2},
      .cells = {NR, NT},
      .mapc2p = mapc2p,

      .cfl_frac = 0.9,

      .num_species = 2,
      .species = {elc, ion},
      .field = {
          .epsilon0 = epsilon0,
          .mu0 = ctx.mu0,
          .mag_error_speed_fact = 1.0,
          .evolve = true,
          .ctx = &ctx,
          .init = init_field,
          .bcx = {GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL},
          .bcy = {GKYL_FIELD_WEDGE, GKYL_FIELD_WEDGE},
      }};

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // print parameters for record
  if (true) {
    double de0 = di0 / sqrt(ctx.mi__me);
    double wpe0 = wpi0 * sqrt(ctx.mi__me);
    double B0 = ctx.vA0 * sqrt(ctx.mu0 * ctx.rho0);
    double Bz0 = B0 * ctx.Bz0__B0;
    double wce0 = -qe * Bz0 / me;
    double wci0 = qi * Bz0 / ctx.mi;
    double Te0 = ctx.T0 / (1 + ctx.Ti0__Te0);
    double Ti0 = ctx.T0 / (1 + 1 / ctx.Ti0__Te0);
    double vte0 = sqrt(Te0 / me);
    double vti0 = sqrt(Ti0 / ctx.mi);
    double rLamor_e0 = vte0 / wce0;
    double rLamor_i0 = vti0 / wci0;
    double dr = (ctx.r_out - ctx.r_inn) / NR;

    printf("%45s = %g\n", "gas_gamma", ctx.gas_gamma);
    printf("%45s = %g\n", "T0", ctx.T0);
    printf("%45s = %g\n", "vTheta0", ctx.vt0__vA0);
    printf("%45s = %g\n", "gravity", ctx.gravity);
    printf("%45s = %g, %g, %g\n", "lightSpeed, epsilon0, mu0", lightSpeed,
           epsilon0, ctx.mu0);
    printf("%45s = %g, %g\n", "mi, qi", ctx.mi, qi);
    printf("%45s = %g, %g\n", "me, qe", me, qe);
    printf("%45s = %g\n", "Ti0/Te0", ctx.Ti0__Te0);
    printf("%45s = %g, %g\n", "di0, de0", di0, de0);
    printf("%45s = %g, %g\n", "rLamor_i0, rLamor_e0", rLamor_i0, rLamor_e0);
    printf("%45s = %g, %g = 1/%g, 1/%g\n", "wpi0, wpe0", wpi0, wpe0, 1 / wpi0,
           1 / wpe0);
    printf("%45s = %g, %g = 1/%g, 1/%g\n", "wci0, wce0 (based on Bz0)", wci0,
           wce0, 1 / wci0, 1 / wce0);
    printf("%45s = %g, %g\n", "r_inn, r_out", ctx.r_inn, ctx.r_out);
    printf("%45s = %d, %g\n", "radial cell # and grid size NR, dr", NR, dr);
    printf("%45s = %g\n", "cfl", ctx.cfl);
    printf("%45s = %g\n", "tend", ctx.tend);
    printf("%45s = %g\n", "ouput time interval", ctx.tend / ctx.nframe);
  }

  // start and end simulation time
  double tcurr = 0.0, tend = ctx.tend;

  // create trigger for IO
  struct gkyl_tm_trigger io_trig = {.dt = tend / ctx.nframe};

  // initialize simulation
  printf("Setting initial conditions ... ");
  gkyl_moment_app_apply_ic(app, tcurr);
  printf("Done\n");
  write_data(&io_trig, app, tcurr);

  // estimate maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  // iterate over the simulation loop
  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    if (step % 1000 == 0)
      printf("Step %6ld, t = %8g, dt = %8g (frame %d)\n", step, tcurr, dt,
             io_trig.curr);

    struct gkyl_update_status status = gkyl_moment_update(app, dt);

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
