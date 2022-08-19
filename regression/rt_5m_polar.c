/* two-fluid five-moment simulation in the r-theta domain       */
/*                                                              */
/* Hakim and Shumlak (2007). Physics of Plasmas, 14(5), 055911. */
/* Section V.B, equations (45) and (46)                         */

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
/*********************************************/

struct moment_ctx {
  double gas_gamma;
  double light_speed;
  double mu0;         // vacuum permeability
  double me;          // electron mass
  double q;           // elementary charge
  double n0;          // the uniform number density
  double wpe0_wce0;   // used to determine B0 (B field changes from -B0 to B0)
  double betae_inn;   // electron plasma beta at inner wall
  double mi__me;      // ion-to-electron mass ratio mi / me
  double Ti0__Te_inn; // ion temperature / electron temperature at inner wall

  double r_0;   // radial center of the B field shear
  double delta; // radial width of the B field shear
  double NN;    // number of points for integration to compute electron pressure
  double pert;  // relative perturbation level

  double r_inn; // radius at inner wall
  double r_out; // radius at outer wall
  int NR;       // radial cell number
  int NT;       // azimuthal cell number

  double tend; // time at end of simulation
  int nframe;  // number of output frames
  double cfl;

  // derived parameters kept for convenience
  double epsilon0;
  double B0;     // reference B field
  double mi;     // ion mass
  double pe_inn; // electron pressure at inner wall
  double pi0;    // uniform ion pressure
};

struct moment_ctx moment_ctx(void) {
  struct moment_ctx ctx = {
      .gas_gamma = 5.0 / 3.0,
      .light_speed = 1.0,
      .mu0 = 1.0,
      .me = 1.0,
      .q = 1.0,
      .n0 = 1.0,
      .wpe0_wce0 = 10.0,
      .betae_inn = 10.0,
      .mi__me = 25.0,
      .Ti0__Te_inn = 1.0,

      .r_0 = 5.0,
      .delta = 0.25,
      .NN = 2560,
      .pert = 1e-2,

      .r_inn = 1.0,
      .r_out = 10.0,
      .NR = 180,
      .NT = 360,

      .tend = 100.0,
      .nframe = 10,
      .cfl = 0.9,
  };

  ctx.epsilon0 = 1.0 / sq(ctx.light_speed) / ctx.mu0;
  ctx.mi = ctx.me * ctx.mi__me;

  double wpe0 = sqrt(ctx.n0 * ctx.q * ctx.q / ctx.me / ctx.epsilon0);
  double wce0 = wpe0 / ctx.wpe0_wce0;
  ctx.B0 = fabs(wce0 * ctx.me / ctx.q);

  double Bz_inn = -ctx.B0;
  double Bz_out = ctx.B0;
  double pmag_inn = 0.5 * Bz_inn * Bz_inn / ctx.mu0;
  ctx.pe_inn = ctx.betae_inn * pmag_inn;
  ctx.pi0 = ctx.pe_inn * ctx.Ti0__Te_inn;

  return ctx;
}

/****************************************************************/
/* INITIAL CONDITIONS                                           */
/****************************************************************/

double calc_Bz(const double r, const struct moment_ctx *ctx) {
  // custom magnetic shear profile
  double B_inn = -ctx->B0;
  double B_out = ctx->B0;
  double r_0 = ctx->r_0;
  double delta = ctx->delta;

  return B_inn + 0.5 * (B_out - B_inn) * (tanh((r - r_0) / delta) + 1.0);
}

double calc_ut(const double r, const struct moment_ctx *ctx) {
  // Eq. (46): u_theta = (dBz/dr) * (e * mu0 * m * n)
  double B_inn = -ctx->B0;
  double B_out = ctx->B0;
  double r_0 = ctx->r_0;
  double delta = ctx->delta;
  double e = ctx->q;
  double mu_0 = ctx->mu0;
  double m = ctx->me;
  double n = ctx->n0;

  double t = tanh((r - r_0) / delta);
  return (1.0 - t * t) * 0.5 * (B_out - B_inn) / (delta * e * m * mu_0 * n);
}

struct gkyl_array *pressure;

void calc_pe_profile(const struct moment_ctx *ctx) {
  // numerically integrate pressure gradient to get the full radial electron
  // pressure profile
  // Eq. (45): dpdr = m * n * u_theta**2 / r - e * n * ut * Bz
  double *data = pressure->data;
  double dr = (ctx->r_out - ctx->r_inn) / ctx->NN;
  data[0] = ctx->pe_inn;
  for (int i = 1; i < ctx->NN; ++i) {
    double r = ctx->r_inn + dr * i;
    double Bz = calc_Bz(r, ctx);
    double ut = calc_ut(r, ctx);

    double e = ctx->q;
    double mu0 = ctx->mu0;
    double m = ctx->me;
    double n = ctx->n0;
    double dpdr = m * n * ut * ut / r - e * n * ut * Bz;

    data[i] = data[i - 1] + dpdr * dr;
  }
}

double calc_pe(const double r, const struct moment_ctx *ctx) {
  double dr = (ctx->r_out - ctx->r_inn) / ctx->NN;
  int idx = floor((r - ctx->r_inn) / dr);
  idx = idx >= 0 ? idx : 0;
  idx = idx <= ctx->NN ? idx : ctx->NN - 1;
  double *data = pressure->data;
  return data[idx];
}

void init_elc(double t, const double *GKYL_RESTRICT xn,
              double *GKYL_RESTRICT fout, void *_ctx) {
  struct moment_ctx *ctx = _ctx;
  double gamma = ctx->gas_gamma;
  double me = ctx->me;
  double n0 = ctx->n0;

  double r = xn[0], theta = xn[1];
  double cost = cos(theta), sint = sin(theta);

  double ure = 0.0;
  double ute = calc_ut(r, ctx);
  double uze = 0.0;

  double rhoe = n0 * me;
  double uxe = ure * cost - ute * sint;
  double uye = ure * sint + ute * cost;
  double pe = calc_pe(r, ctx);

  double r_0 = ctx->r_0;
  double delta = ctx->delta;
  pe *= 1.0 + ctx->pert * (1.0 - sq(tanh((r - r_0) / delta))) * cos(theta);

  fout[0] = rhoe;
  fout[1] = rhoe * uxe;
  fout[2] = rhoe * uye;
  fout[3] = rhoe * uze;
  fout[4] = pe / (gamma - 1.0) + 0.5 * rhoe * (sq(uxe) + sq(uye) + sq(uze));
}

void init_ion(double t, const double *GKYL_RESTRICT xn,
              double *GKYL_RESTRICT fout, void *_ctx) {
  struct moment_ctx *ctx = _ctx;
  double gamma = ctx->gas_gamma;
  double mi = ctx->mi;
  double n0 = ctx->n0;
  double pi0 = ctx->pi0;

  double r = xn[0], theta = xn[1];

  double rhoi = n0 * mi;
  double uxi = 0.0;
  double uyi = 0.0;
  double uzi = 0.0;
  double pi = pi0;

  fout[0] = rhoi;
  fout[1] = rhoi * uxi;
  fout[2] = rhoi * uyi;
  fout[3] = rhoi * uzi;
  fout[4] = pi / (gamma - 1.0) + 0.5 * rhoi * (sq(uxi) + sq(uyi) + sq(uzi));
}

void init_field(double t, const double *GKYL_RESTRICT xn,
                double *GKYL_RESTRICT fout, void *_ctx) {
  struct moment_ctx *ctx = _ctx;
  double r = xn[0], theta = xn[1];

  double Ex = 0.0;
  double Ey = 0.0;
  double Ez = 0.0;
  double Bx = 0.0;
  double By = 0.0;
  double Bz = calc_Bz(r, ctx);

  fout[0] = Ex;
  fout[1] = Ey;
  fout[2] = Ez;
  fout[3] = Bx;
  fout[4] = By;
  fout[5] = Bz;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

/*********************************************/
/* MISC.                                     */
/*********************************************/

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

  if (true) {
    printf("%45s = %g\n", "light_speed", ctx.light_speed);
    printf("%45s = %g\n", "mu0", ctx.mu0);
    printf("%45s = %g\n", "epsilon0", ctx.epsilon0);
    printf("%45s = %g\n", "q", ctx.q);
    printf("%45s = %g\n", "me", ctx.me);
    printf("%45s = %g\n", "n0", ctx.n0);
    printf("%45s = %g\n", "B0", ctx.B0);
    printf("%45s = %g\n", "pe_inn", ctx.pe_inn);
    printf("%45s = %g\n", "magnetic shear location r_0", ctx.r_0);
    printf("%45s = %g\n", "magnetic shear width delta", ctx.delta);
    printf("%45s = %g\n", "tend", ctx.tend);
    printf("%45s = %d\n", "nframe", ctx.nframe);
  }

  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(ctx.gas_gamma);

  struct gkyl_moment_species elc = {
      .name = "elc",
      .charge = -ctx.q,
      .mass = ctx.me,
      .equation = euler,
      .evolve = true,
      .ctx = &ctx,
      .init = init_elc,
      .bcx = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT},
  };

  struct gkyl_moment_species ion = {
      .name = "ion",
      .charge = ctx.q,
      .mass = ctx.mi,
      .equation = euler,
      .evolve = true,
      .ctx = &ctx,
      .init = init_ion,
      .bcx = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT},
  };

  int NR = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.NR);
  int NT = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.NT);

  struct gkyl_moment app_inp = {
      .name = "5m_polar",

      .ndim = 2,
      .lower = {ctx.r_inn, 0.0},
      .upper = {ctx.r_out, 2.0 * M_PI},
      .cells = {NR, NT},
      .mapc2p = mapc2p,
      .num_periodic_dir = 1,
      .periodic_dirs = {1},

      .cfl_frac = ctx.cfl,

      .num_species = 1,
      .species = {elc, ion},
      .field = {
          .epsilon0 = ctx.epsilon0,
          .mu0 = ctx.mu0,
          .mag_error_speed_fact = 1.0,
          .evolve = true,
          .ctx = &ctx,
          .init = init_field,
          .bcx = {GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL},
      }};

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start and end simulation time
  double tcurr = 0.0, tend = ctx.tend;

  // create trigger for IO
  struct gkyl_tm_trigger io_trig = {.dt = tend / ctx.nframe};

  // initialize simulation
  printf("Setting initial conditions ...\n");
  printf("    Computing equilibrium electron pressure ...\n");
  pressure = gkyl_array_new(GKYL_DOUBLE, 1, ctx.NN);
  calc_pe_profile(&ctx);
  gkyl_moment_app_apply_ic(app, tcurr);
  gkyl_array_release(pressure);
  printf("Done\n");
  write_data(&io_trig, app, tcurr);

  // estimate maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  // iterate over the simulation loop
  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    if (step % 10 == 0)
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
