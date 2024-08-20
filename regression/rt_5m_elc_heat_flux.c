#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

struct heat_flux_ctx {
  // parameters for plasma
  double gasGamma;
  double epsilon0;
  double mu0;
  double charge;
  double elcMass;
  double ionMass;
  double n;
  double T_high;
  double T_low;
  double tau;
  double lambda;
  double Lx;
};

struct heat_flux_ctx
create_ctx(void)
{
  double gasGamma = 5.0/3.0;
  double charge = 1.602176487e-19;
  double elcMass = 9.10938215e-31;
  double ionMass = 1.672621637e-27;
  double epsilon0 = 8.854187817620389850536563031710750260608e-12;
  double mu0 = 12.56637061435917295385057353311801153679e-7;

  // plasma density in particles per m^3
  double n = 1.0e18;
  // plasma temperature in eV
  double T_high = 100.0*charge;
  double T_low = 1.0*charge;
  double vte = sqrt(T_low/elcMass);
  
  // Box size is chosen to be multiple mean-free-paths long
  double tau = 6.0*sqrt(2.0*M_PI*elcMass*T_low*M_PI*T_low*M_PI*T_low)*epsilon0*epsilon0/(charge*charge*charge*charge*n);
  double lambda = vte*tau;
  double Lx = 1000.0*lambda;

  struct heat_flux_ctx ctx = {
    .gasGamma = gasGamma,
    .charge = charge,
    .elcMass = elcMass,
    .ionMass = ionMass,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .n = n,
    .T_high = T_high,
    .T_low = T_low,
    .tau = tau,
    .lambda = lambda,
    .Lx = Lx,
  };
  return ctx;
}

void
evalElcInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  struct heat_flux_ctx *app = ctx;

  double T = (app->T_high + (app->T_low-app->T_high)*x/app->Lx);

  double rhoe = app->n*app->elcMass;
  double ere = app->n*T/(app->gasGamma-1);

  fout[0] = rhoe;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = ere;  
}

void
evalIonInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  struct heat_flux_ctx *app = ctx;

  double T = app->T_low;

  double rhoi = app->n*app->ionMass;
  double eri = app->n*T/(app->gasGamma-1);

  fout[0] = rhoi;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = eri;   
}

void
evalFieldInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  // electric field
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // magnetic field
  fout[3] = 0.0, fout[4] = 0.0; fout[5] = 0.0;
  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalElcBcLower(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  struct heat_flux_ctx *app = ctx;  
  double T = app->T_high;

  double rhoe = app->n*app->elcMass;
  double ere = app->n*T/(app->gasGamma-1);

  ghost[0] = rhoe;
  ghost[1] = 0.0; ghost[2] = 0.0; ghost[3] = 0.0;
  ghost[4] = ere;   
}

void
evalElcBcUpper(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  struct heat_flux_ctx *app = ctx;  
  double T = app->T_low;

  double rhoe = app->n*app->elcMass;
  double ere = app->n*T/(app->gasGamma-1);

  ghost[0] = rhoe;
  ghost[1] = 0.0; ghost[2] = 0.0; ghost[3] = 0.0;
  ghost[4] = ere;   
}

void
write_data(struct gkyl_tm_trigger *iot, const gkyl_moment_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr))
    gkyl_moment_app_write(app, tcurr, iot->curr-1);
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);
  struct heat_flux_ctx ctx = create_ctx(); // context for init functions
  // electron/ion equations
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(ctx.gasGamma, app_args.use_gpu);
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(ctx.gasGamma, app_args.use_gpu);

  printf("Collision time = %lg\n", ctx.tau);
  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = -1.0*ctx.charge, .mass = ctx.elcMass,

    .equation = elc_euler,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalElcInit,
    .type_brag = GKYL_BRAG_UNMAG_FULL,    

    .bcx = { GKYL_SPECIES_FUNC, GKYL_SPECIES_FUNC },  
    .bcx_func = { evalElcBcLower, evalElcBcUpper}, 
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.charge, .mass = ctx.ionMass,

    .equation = ion_euler,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalIonInit,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };  

  // Field.
  struct gkyl_moment_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .mag_error_speed_fact = 1.0,
    
    .evolve = false,
    .init = evalFieldInit,
    .ctx = &ctx,
  };

  // VM app
  struct gkyl_moment app_inp = {
    .name = "5m_elc_heat_flux",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx }, 
    .cells = { 64 },

    .num_periodic_dir = 0,
    .periodic_dirs = { },
    .cfl_frac = 0.05,

    .has_braginskii = true, 
    .coll_fac = 1.0,  

    .field = field,    

    .num_species = 2,
    .species = { elc, ion },
  };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 1.0*ctx.tau;
  int nframe = 1;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  gkyl_moment_app_write(app, tcurr, 0);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1;
  while ((tcurr < tend) && (step <= app_args.num_steps)) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
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

  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Source updates took %g secs\n", stat.sources_tm);
  printf("Total updates took %g secs\n", stat.total_tm);

  // simulation complete, free resources
  gkyl_wv_eqn_release(elc_euler);
  gkyl_wv_eqn_release(ion_euler);
  gkyl_moment_app_release(app);  
  
  return 0;
}
