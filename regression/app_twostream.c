#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_vlasov.h>

struct twostream_ctx {
    double knumber; // wave-number
    double vth; // electron thermal velocity
    double vdrift; // drift velocity
    double perturbation;
};

static inline double sq(double x) { return x*x; }

void
evalDistFunc(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  struct twostream_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double alpha = app->perturbation, k = app->knumber, vt = app->vth, vdrift = app->vdrift;

  double fv = 1/sqrt(8*M_PI*sq(vt))*(exp(-sq(v-vdrift)/(2*sq(vt)))+exp(-sq(v+vdrift)/(2*sq(vt))));
  fout[0] = (1+alpha*cos(k*x))*fv;
}

void
evalFieldFunc(double t, const double* restrict xn, double* restrict fout, void *ctx)
{
  struct twostream_ctx *app = ctx;
  double x = xn[0];
  double alpha = app->perturbation, k = app->knumber;

  double E_x = -alpha*sin(k*x)/k;
  
  fout[0] = E_x; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

struct twostream_ctx
create_ctx(void)
{
  struct twostream_ctx ctx = {
    .knumber = 0.5,
    .vth = 0.2,
    .vdrift = 1.0,
    .perturbation = 1.0e-6
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct twostream_ctx ctx = create_ctx(); // context for init functions

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = -1.0, .mass = 1.0,
    .lower = { -6.0 },
    .upper = { 6.0 }, 
    .cells = { 32 },

    .evolve = 1,
    .ctx = &ctx,
    .init = evalDistFunc,

    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2" },
  };

  // field
  struct gkyl_em_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalFieldFunc
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "twostream",

    .cdim = 1, .vdim = 1,
    .lower = { -M_PI/ctx.knumber },
    .upper = { M_PI/ctx.knumber },
    .cells = { 64 },
    .poly_order = 2,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 1,
    .species = { elc },
    .field = field
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 40.0;
  double dt = tend-tcurr;

  // initialize simulation
  gkyl_vlasov_app_init(app, tcurr);
  
  gkyl_vlasov_app_write(app, tcurr, 0);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 0);

  int count = 0;
  while (tcurr < tend) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_vlasov_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    count++;
  }

  gkyl_vlasov_app_write(app, tcurr, 1);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 1);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  // simulation complete, free app
  gkyl_vlasov_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of forward-Euler calls %ld\n", stat.nfeuler);
  printf("Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  printf("Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  printf("Species RHS calc took %g secs\n", stat.species_rhs_tm);
  printf("Field RHS calc took %g secs\n", stat.field_rhs_tm);
  printf("Current evaluation and accumulate took %g secs\n", stat.current_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
