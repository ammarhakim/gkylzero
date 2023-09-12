#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_iso_euler.h>
#include <rt_arg_parse.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <gkyl_mpi_comm.h>
#include <mpi.h>
#endif

#define RHO (0)
#define RHOVX (1)
#define RHOVY (2)
#define RHOVZ (3)
#define ENERGY (4)

#define EX (0)
#define EY (1)
#define EZ (2)
#define BX (3)
#define BY (4)
#define BZ (5)
#define PSI_E (6)
#define PSI_B (7)

#define sq(x) ((x) * (x))
#define cube(x) ((x) * (x) * (x))

/**************************************************************/
/** PARAMETERS                                               **/
/**************************************************************/

struct problem_context {
  double eV, kB, Te0, Ti0, mi, B, n0, me;
  double c, mu0, epsilon0, gamma, qi, qe;
  bool force_first_order;

  bool has_collision;
  double nu_base_ei, nu_base_ee, nu_base_ii;

  double pert;
  pcg64_random_t rng;

  double cs0, rhos0; // for additional density and temperature sources

  double Lperp, r_inn, r_out;
  double tend, cfl;
  int nframe;
};

void setup_problem_context(struct problem_context *ctx, const bool verbose) {
  // basic parameters in table 4.2; the orignal values are:
  // Te0=6eV, Ti0=1eV, mi=3.973mp, B=0.0398T, n0=2e18/m^3, mi/me=400
  // double eV = ctx->eV = GKYL_EV2KELVIN;
  // double kB = ctx->kB = GKYL_BOLTZMANN_CONSTANT;
  double eV = ctx->eV = GKYL_ELEMENTARY_CHARGE;
  double kB = ctx->kB = 1.0;
  double Te0 = ctx->Te0 = 6.0 * eV;
  double Ti0 = ctx->Ti0 = 1.0 * eV;
  double mi = ctx->mi = 3.973 * GKYL_PROTON_MASS;
  double B = ctx->B = 0.0398;
  double n0 = ctx->n0 = 2e18;
  double me = ctx->me = mi / 400.0;

  double c = GKYL_SPEED_OF_LIGHT / 30.0; // FIXME reduced value for larger dt
  double mu0 = ctx->mu0 = GKYL_MU0;
  double epsilon0 = ctx->epsilon0 = 1 / mu0 / sq(c);
  double gamma = ctx->gamma = 5.0 / 3.0; // FIXME 2 or 5/3?
  double qi = ctx->qi = GKYL_ELEMENTARY_CHARGE;
  double qe = ctx->qe = -qi;

  // simulation control parameters
  double cs0 = ctx->cs0 = sqrt(kB * Te0 / mi); // table 4.2, 1.2e4 m/s
  double wci0 = qi * B / mi;                   // table 4.2, 9.6e5 rad/s
  double rhos0 = ctx->rhos0 = cs0 / wci0;      // table 4.2, 1.25 cm
  double Lperp = ctx->Lperp = 100.0 * rhos0; // table 4.2, 100 * rhos0 ~ 1.25 m
  double r_inn = ctx->r_inn = Lperp * 0.025; // inner radial boundary location
  double r_out = ctx->r_out = Lperp * 0.5;   // outer radial boundary location
  double cfl = ctx->cfl = 0.9;    // Shi used cfl=0.1 for the DG GK run
  double tend = ctx->tend = 1e-3; // Shi figue 4.5, frames at t=0.4, 1, 1.5 ms
  int nframe = ctx->nframe = 100000;
  bool force_first_order = ctx->force_first_order = true;

  // more derived parameters for verification/diagnostics
  double de0 = sqrt(me / (n0 * sq(qe) * mu0));
  double wpe0 = c / de0;
  double wce0 = -qe * B / me;
  double vte0 = sqrt(kB * Te0 / me);
  double larmor_e0 = vte0 / wce0;
  double vti0 = sqrt(kB * Ti0 / mi);
  double larmor_i0 = vti0 / wci0;
  double di0 = sqrt(mi / (n0 * sq(qi) * mu0));
  double wpi0 = c / di0;
  double vAe0 = B / sqrt(mu0 * n0 * me);
  double vAi0 = B / sqrt(mu0 * n0 * mi);
  double vA0 = B / sqrt(mu0 * n0 * (mi + me));

  // reference collision frequencies at a given reference density
  ctx->has_collision = false;
  double nu_ei0 = wci0 * 0;
  double nu_ie0 = nu_ei0 * me / mi;
  // scaling factors for collision frequencies; during the simulation, one has
  // nu_sr=base*rho_r and nu_rs=base*rho_s to satisfy
  // nu_sr*rho_s=nu_rs*rho_r=base
  ctx->nu_base_ei = nu_ei0 / (n0 * mi);
  ctx->nu_base_ee = 0;
  ctx->nu_base_ii = 0;

  double pert = ctx->pert = 0.1; // random noise perturbation level
  ctx->rng = gkyl_pcg64_init(0); // random number generator

  if (verbose) {
    printf("%30s = %.1e\n", "cs0", cs0);
    printf("%30s = %.1e\n", "wci0", wci0);
    printf("%30s = %.2f cm\n", "rhos0 = sqrt(Te0/me) / wci0", rhos0 * 100.0);
    // following must match values in the section 4
    printf("%30s = %g\n", "mi/me", mi / me);
    printf("%30s = %gMHz, %gGHz\n", "wpi0, wpe0", wpi0 / 1e6, wpe0 / 1e9);
    printf("%30s = %gkHz, %gGHz\n", "wci0, wce0", wci0 / 1e3, wce0 / 1e9);
    printf("%30s = %gcm, %gcm\n", "larmor_i0, larmor_e0", larmor_i0 * 1e2,
      larmor_e0 * 1e2);
    printf("%30s = %gcm/s, %gcm/s\n", "vti0, vte0", vti0 * 1e2, vte0 * 1e2);
    printf("%30s = %gcm, %gcm\n", "di0, de0", di0 * 1e2, de0 * 1e2);
    printf("%30s = %gcm/s, %gcm/s, %gcm/s\n", "vA0, vAi0, vAe0", vA0 * 1e2,
      vAi0 * 1e2, vAe0 * 1e2);
    printf("%30s = %g\n", "Lperp/c", Lperp / c);
    printf(
      "%30s = %g, %g\n", "Lperp/vte0, Lperp/vti0", Lperp / vte0, Lperp / vti0);
    printf(
      "%30s = %.2f = %g = %g rhos0\n", "Lperp", Lperp, Lperp, Lperp / rhos0);
    // following values must be reasonable
    printf("%30s = %g de0 = %g di0\n", "", Lperp / de0, Lperp / di0);
    printf("%30s = %g vte0/wce0 = %g vti0/wci0\n", "", Lperp / larmor_e0,
      Lperp / larmor_i0);
    printf("%30s = %.2f = %g Lperp\n", "r_inn", r_inn, r_inn / Lperp);
    printf("%30s = %.2f = %g Lperp\n", "r_out", r_out, r_out / Lperp);
    printf("%30s = %g = %g real light speed = %g cs0 = %g vte0\n", "c", c,
      c / GKYL_SPEED_OF_LIGHT, c / cs0, c / vte0);
    printf("%30s = %g vAe0 = %g vAi0\n", "", c / vAe0, c / vAi0);
    printf("%30s = %g\n", "cfl", cfl);
    printf("%30s = %g = %g ms\n", "tend", tend, tend * 1e3);
    printf("%30s = %d\n", "nframe", nframe);
    printf("%30s = %g\n", "tend/nframe", tend / nframe);
    printf(
      "%30s = %s\n", "force_first_order", force_first_order ? "yes" : "no");
    printf("%30s = %g\n", "pert", pert);
    printf("%30s = %g kHz = wci0/%g = wpi0/%g\n", "nu_ei0", nu_ei0 / 1e3,
      wci0 / nu_ei0, wpi0 / nu_ei0);
    printf("%30s = %g kHz = wci0/%g = wpi0/%g\n", "nu_ie0", nu_ie0 / 1e3,
      wci0 / nu_ie0, wpi0 / nu_ie0);

    printf("\n");
    printf("%30s = %g, %g\n", "eV, kB", eV, kB);
    printf("%30s = %g\n", "gamma", gamma);
    printf("%30s = %g, %g\n", "mi, me", mi, me);
    printf("%30s = %g PROTON_MASS, %g ELEMENTARY_CHARGE\n", "",
      mi / GKYL_PROTON_MASS, me / GKYL_ELEMENTARY_CHARGE);
    printf("%30s = %g = %g ELEMENTARY_CHARGE\n", "qi", qi,
      qi / GKYL_ELEMENTARY_CHARGE);
    printf(
      "%30s = %g, %g = %geV, %geV\n", "Ti0, Te0", Ti0, Te0, Ti0 / eV, Te0 / eV);
    printf("%30s = %g, %g\n", "kB*Ti0, kB*Te0", kB * Ti0, kB * Te0);
    printf("%30s = %g = %gcc\n", "n0", n0, n0 / 1e6);
    printf("%30s = %g = %g G\n", "B", B, B * 1e4);
    printf("%30s = %g\n", "rhos0", rhos0);
  }
}

/**************************************************************/
/** GRID                                                     **/
/**************************************************************/

// Cornerific Tapered2 Mapping
// https://gist.github.com/liangwang0734/29592325d5cfc28f63e865b6abbda81f
void mapc2p(
  double t, const double *xc, double *GKYL_RESTRICT xp, void *ctx_raw) {
  struct problem_context *ctx = ctx_raw;

  double u = xc[0], v = xc[1];

  double tmp = sqrt((u * u + v * v - 2 * u * u * v * v) / (u * u + v * v) /
                    (1 - u * u * v * v));

  if ((fabs(u) < 1e-10) && (fabs(v) < 1e-10))
    tmp = 0.0;

  if (fabs(fabs(u * v) - 1) < 1e-10)
    tmp = sqrt(1.0 / 2.0);

  double x = u * tmp * ctx->Lperp / 2.0;
  double y = v * tmp * ctx->Lperp / 2.0;

  xp[0] = x;
  xp[1] = y;
}

/**************************************************************/
/** INITIAL CONDITIONS                                       **/
/**************************************************************/

// eq 4.2 for initial density and temperature, as well as electron temperature
// source
double A(double r, double cedge, struct problem_context *ctx) {
  double r0 = ctx->Lperp / 2.0;
  if (r < r0)
    return (1 - cedge) * cube(1 - sq(r) / sq(r0)) + cedge;
  else
    return cedge;
}

void init_elc(double t, const double *GKYL_RESTRICT xn,
  double *GKYL_RESTRICT fout, void *ctx_raw) {
  struct problem_context *ctx = ctx_raw;

  double xnp[3];
  mapc2p(t, xn, xnp, ctx_raw);
  double x = xnp[0], y = xnp[1];
  double r = sqrt(x * x + y * y);
  double phi = atan2(x, y);

  double n = ctx->n0 * A(r, 1.0 / 20.0, ctx);      // above eq 4.2
  double T = 5.7 * A(r, 1.0 / 5.0, ctx) * ctx->eV; // below eq 4.2
  double vx = 0.0, vy = 0.0, vz = 0.0;

  T = T * (1 + ctx->pert * gkyl_pcg64_rand_double(&ctx->rng)); // FIXME

  double rho = n * ctx->me;
  double p = n * ctx->kB * T;
  double E = p / (ctx->gamma - 1.0) + 0.5 * rho * (vx * vx + vy * vy + vz * vz);

  fout[RHO] = rho;
  fout[RHOVX] = rho * vx;
  fout[RHOVY] = rho * vy;
  fout[RHOVZ] = rho * vz;
  fout[ENERGY] = E;
}

void init_ion(double t, const double *GKYL_RESTRICT xn,
  double *GKYL_RESTRICT fout, void *ctx_raw) {
  struct problem_context *ctx = ctx_raw;

  double xnp[3];
  mapc2p(t, xn, xnp, ctx_raw);
  double x = xnp[0], y = xnp[1];
  double r = sqrt(x * x + y * y);
  double phi = atan2(x, y);

  double n = ctx->n0 * A(r, 1.0 / 20.0, ctx); // above eq 4.2
  double T = 1.0 * ctx->eV;                   // below eq 4.2
  double vx = 0.0, vy = 0.0, vz = 0.0;

  double rho = n * ctx->mi;
  double p = n * ctx->kB * T;
  double E = p / (ctx->gamma - 1.0) + 0.5 * rho * (vx * vx + vy * vy + vz * vz);

  fout[RHO] = rho;
  fout[RHOVX] = rho * vx;
  fout[RHOVY] = rho * vy;
  fout[RHOVZ] = rho * vz;
  fout[ENERGY] = E;
}

void init_emf(
  double t, const double *restrict xn, double *restrict fout, void *ctx_raw) {
  struct problem_context *ctx = ctx_raw;

  fout[EX] = 0.0;
  fout[EY] = 0.0;
  fout[EZ] = 0.0;
  fout[BX] = 0.0;
  fout[BY] = 0.0;
  fout[BZ] = ctx->B;
  fout[PSI_E] = 0.0;
  fout[PSI_B] = 0.0;
}

/**************************************************************/
/** SOURCES                                                  **/
/**************************************************************/

// eq 4.3, the S without the n0 term
double S(double r, struct problem_context *ctx) {
  double n0 = ctx->n0;
  double cs0 = ctx->cs0;
  double rs = 20.0 * ctx->rhos0;
  double Ls = 0.5 * ctx->rhos0;
  double Lz = 1440.0 * ctx->rhos0; // table 4.2
  return (1.08 * n0 * cs0 / Lz) *
         (0.01 + 0.99 * 0.5 * (1 - tanh((r - rs) / Ls)));
}

void calc_elc_nT_source(double t, const double *GKYL_RESTRICT xn,
  double *GKYL_RESTRICT fout, void *ctx_raw) {
  struct problem_context *ctx = ctx_raw;

  double r = xn[0], phi = xn[1];

  fout[0] = S(r, ctx); // number density source
  // source for scaled temperature
  double fact = ctx->kB / (ctx->gamma - 1);
  fout[1] = fact * 6.8 * A(r, 1.0 / 2.5, ctx) * ctx->eV;
}

void calc_ion_nT_source(double t, const double *GKYL_RESTRICT xn,
  double *GKYL_RESTRICT fout, void *ctx_raw) {
  struct problem_context *ctx = ctx_raw;

  double r = xn[0], phi = xn[1];

  fout[0] = S(r, ctx); // number density source
  // source for scaled temperature
  double fact = ctx->kB / (ctx->gamma - 1);
  fout[1] = fact * 1.0 * ctx->eV;
}

/*************************************************************/
/** BOUNDARY CONDITIIONS                                    **/
/*************************************************************/

/**************************************************************/
/** IO                                                       **/
/**************************************************************/

void write_data(
  struct gkyl_tm_trigger *iot, const gkyl_moment_app *app, double tcurr) {
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_moment_app_cout(
      app, stdout, ">>>> Writing output frame %d...", iot->curr - 1);
    gkyl_moment_app_write(app, tcurr, iot->curr - 1);
    gkyl_moment_app_cout(app, stdout, " Done.\n");
  }
}

/**************************************************************/
/** MAIN BODY OF THE PROGRAM                                 **/
/**************************************************************/

int main(int argc, char **argv) {
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 31);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 31);

  // create global range
  struct gkyl_range globalr;
  gkyl_create_global_range(2, (int[]){NX, NY}, &globalr);

#ifdef GKYL_HAVE_MPI

  // create decomposition
  int cuts[] = {1, 1};
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
    cuts[0] = app_args.cuts[0];
    cuts[1] = app_args.cuts[1];
  }

  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(2, cuts, &globalr);

  // construct communcator for use in app
  struct gkyl_comm *comm;
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new(&(struct gkyl_mpi_comm_inp){
      .mpi_comm = MPI_COMM_WORLD, .decomp = decomp});
  } else {
    comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){
      .decomp = decomp,
    });
  }

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_sz;
  gkyl_comm_get_size(comm, &comm_sz);

  int ncuts = cuts[0] * cuts[1];
  if (ncuts != comm_sz) {
    if (my_rank == 0)
      fprintf(stderr, "*** Number of ranks, %d, do not match total cuts, %d!\n",
        comm_sz, ncuts);
    goto mpifinalize;
  }

#else

  // create decomposition
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(2, (int[]){1, 1}, &globalr);

  // construct communcator for use in app
  struct gkyl_comm *comm;
  comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){
    .decomp = decomp,
  });

  int my_rank = 0;

#endif

  // set problem-specifc parameters
  struct problem_context ctx;
  setup_problem_context(&ctx, my_rank == 0 ? true : false);

  // create equation
  struct gkyl_wv_euler_inp inp = {
    .gas_gamma = ctx.gamma,
    .rp_type = WV_EULER_RP_HLLC,
  };

  // create electron(elc), ion species
  struct gkyl_wv_eqn *elc_eqn = gkyl_wv_euler_inew(&inp);
  struct gkyl_moment_species elc = {
    .name = "elc",
    .mass = ctx.me,
    .charge = ctx.qe,
    .equation = elc_eqn,
    .ctx = &ctx,
    .init = init_elc,
    .bcx = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT},
    .bcy = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT},
    .nT_source_func = calc_elc_nT_source,
  };

  struct gkyl_wv_eqn *ion_eqn = gkyl_wv_euler_inew(&inp);
  struct gkyl_moment_species ion = {
    .name = "ion",
    .mass = ctx.mi,
    .charge = ctx.qi,
    .equation = ion_eqn,
    .ctx = &ctx,
    .init = init_ion,
    .bcx = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT},
    .bcy = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT},
    .nT_source_func = calc_ion_nT_source,
  };

  // create app
  struct gkyl_moment app_inp = {.name = "lapd",

    .ndim = 2,
    .lower = {-1.0, -1.0}, // lower bounds of computational domain
    .upper = {1.0, 1.0},   // upper bounds of computational domain
    .cells = {NX, NY},     // # of cells in computational domain
    .mapc2p = mapc2p,      // mapping computational -> physical coords
    .c2p_ctx = &ctx,
    .num_periodic_dir = 0, // number of periodic directions
    .cfl_frac = ctx.cfl,   // cfl number

    .num_species = 2,
    .species = {elc, ion},

    .field =
      {
        .epsilon0 = ctx.epsilon0,
        .mu0 = ctx.mu0,
        .ctx = &ctx,
        .init = init_emf,
        .bcx = {GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL},
        .bcy = {GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL},
      },

    .has_collision = ctx.has_collision,
    // scaling factors for collision frequencies; nu_sr = nu_base_sr * rho_r
    // the matrix must be symmetrix
    .nu_base =
      {
        {0, ctx.nu_base_ei}, // ee, ei
        {ctx.nu_base_ei, 0}  // ie, ii
      },

    .has_low_inp = true,
    .low_inp = {.local_range = decomp->ranges[my_rank], .comm = comm}

  };
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);
  gkyl_moment_app_cout(app, stdout, "%30s = %d, %d\n", "NX, NY", NX, NY);

  // set start and end time
  double tcurr = 0.0, tend = ctx.tend;

  // number of output frames and output trigger
  int nframe = ctx.nframe;
  struct gkyl_tm_trigger io_trig = {.dt = tend / nframe};

  // initialize simulation and write the initial conditions
  gkyl_moment_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

  // compute estimate of maximum stable time-step for the first step;
  // FIXME perhaps wrong for mapped grid
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    bool show_log = step % 100 == 0 || step == 1;
    if (show_log)
      gkyl_moment_app_cout(
        app, stdout, "Frame %d step %ld t=%-10g", io_trig.curr, step, tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    if (show_log) {
      gkyl_moment_app_cout(app, stdout, " dt=%g\n", status.dt_actual);
    }

    if (!status.success) {
      gkyl_moment_app_cout(
        app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }
  if (io_trig.curr == 1)
    gkyl_moment_app_write(app, tcurr, 1);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);
  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(
    app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(
    app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(
    app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(
    app, stdout, "Total updates took %g secs\n", stat.total_tm);

  // simulation complete, free resources
  gkyl_wv_eqn_release(elc_eqn);
  gkyl_wv_eqn_release(ion_eqn);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);

mpifinalize:;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
