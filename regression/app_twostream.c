#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_vlasov.h>

struct twostream_ctx {
    gkyl_real knumber; // wave-number
    gkyl_real vth; // electron thermal velocity
    gkyl_real vdrift; // drift velocity
    gkyl_real perturbation;
};

static inline gkyl_real sq(gkyl_real x) { return x*x; }

void
evalDistFunc(gkyl_real t, const gkyl_real * restrict xn, gkyl_real* restrict fout, void *ctx)
{
  struct twostream_ctx *app = ctx;
  gkyl_real x = xn[0], v = xn[1];
  gkyl_real alpha = app->perturbation, k = app->knumber, vt = app->vth, vdrift = app->vdrift;

  gkyl_real fv = 1/sqrt(8*M_PI*sq(vt))*(exp(-sq(v-vdrift)/(2*sq(vt)))+exp(-sq(v+vdrift)/(2*sq(vt))));
  fout[0] = (1+alpha*cos(k*x))*fv;
}

void
evalFieldFunc(gkyl_real t, const gkyl_real* restrict xn, gkyl_real* restrict fout, void *ctx)
{
  struct twostream_ctx *app = ctx;
  gkyl_real x = xn[0];
  gkyl_real alpha = app->perturbation, k = app->knumber;

  gkyl_real E_x = -alpha*sin(k*x)/k;
  
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
  gkyl_real tcurr = 0.0, tend = 40.0;
  gkyl_real dt = tend-tcurr;

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  
  gkyl_vlasov_app_write(app, tcurr, 0);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 0);

  while (tcurr < tend) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
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
