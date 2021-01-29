#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_vlasov.h>

struct esshock_ctx {
    gkyl_real chargeElc; // electron charge
    gkyl_real massElc; // electron mass
    gkyl_real chargeIon; // ion charge
    gkyl_real massIon; // ion mass
    gkyl_real Te_Ti; // electron to ion temperature ratio
    gkyl_real vte; // electron thermal velocity
    gkyl_real vti; // ion thermal velocity
    gkyl_real cs; // sound speed
    gkyl_real uShock; // in-flow velocity
    gkyl_real Lx; // size of the box
};

static inline gkyl_real sq(gkyl_real x) { return x*x; }

void
evalDistFuncElc(gkyl_real t, const gkyl_real * restrict xn, gkyl_real* restrict fout, void *ctx)
{
  struct esshock_ctx *app = ctx;
  gkyl_real x = xn[0], v = xn[1];
  gkyl_real vt = app->vte, vdrift = app->uShock;
  gkyl_real fv = 0.0;
  if (x < 0)
    fv = 1.0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v-vdrift)/(2*sq(vt))));
  else
    fv = 1.0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v+vdrift)/(2*sq(vt))));
  fout[0] = fv;
}

void
evalDistFuncIon(gkyl_real t, const gkyl_real * restrict xn, gkyl_real* restrict fout, void *ctx)
{
  struct esshock_ctx *app = ctx;
  gkyl_real x = xn[0], v = xn[1];
  gkyl_real vt = app->vti, vdrift = app->uShock;
  gkyl_real fv = 0.0;
  if (x < 0)
    fv = 1.0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v-vdrift)/(2*sq(vt))));
  else
    fv = 1.0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v+vdrift)/(2*sq(vt))));
  fout[0] = fv;
}

void
evalFieldFunc(gkyl_real t, const gkyl_real* restrict xn, gkyl_real* restrict fout, void *ctx)
{
  struct esshock_ctx *app = ctx;
  gkyl_real x = xn[0];
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

struct esshock_ctx
create_ctx(void)
{
  struct esshock_ctx ctx = {
    .chargeElc = -1.0,
    .massElc = 1.0,
    .chargeIon = 1.0,
    .massIon = 1836.153,
    .Te_Ti = 4.0,
    .vte = 1.0,
    // .vti = 0.01,
    // .cs = 1.0/sqrt(1836.153),
    // .uShock = 0.05,
    .vti = ctx.vte/sqrt(ctx.Te_Ti*ctx.massIon),
    .cs = ctx.vte/sqrt(ctx.massIon),
    .uShock = 2.0*ctx.cs,
    .Lx = 128.0
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct esshock_ctx ctx = create_ctx(); // context for init functions

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -6.0 * ctx.vte},
    .upper = { 6.0 * ctx.vte}, 
    .cells = { 64 },

    .evolve = 1,
    .ctx = &ctx,
    .init = evalDistFuncElc,

    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2" },
  };

  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -16.0 * ctx.vti},
    .upper = { 16.0 * ctx.vti}, 
    .cells = { 64 },

    .evolve = 1,
    .ctx = &ctx,
    .init = evalDistFuncIon,

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
    .name = "esshock",

    .cdim = 1, .vdim = 1,
    .lower = { -ctx.Lx },
    .upper = { ctx.Lx },
    .cells = { 256 },
    .poly_order = 2,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, ion },
    .field = field
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(vm);

  // start, end and initial time-step
  gkyl_real tcurr = 0.0, tend = 20.0;
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
