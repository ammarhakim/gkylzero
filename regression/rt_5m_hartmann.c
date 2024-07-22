#include <math.h>
#include <stdio.h>

#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>
#include <rt_arg_parse.h>

struct hartmann_ctx {
  // parameters for plasma
  double gasGamma;
  double charge;
  double elcMass, ionMass;
  double T_elc, T_ion;
  double beta;
  double B0;
  double n;
  double Lx;

  double epsilon0, mu0;

  double rel_fac, coll_fac;
  double tau, lambda;
  double eta_par, eta_perp;
  double hartmann_num;
  double accel, force_norm;
};

struct hartmann_ctx
create_ctx(void)
{
  double gasGamma = 5.0/3.0;
  double charge = 1.0;
  double ionMass = 1.0;
  double elcMass = ionMass/sqrt(1836.153);
  double epsilon0 = 1.0;
  double mu0 = 1.0;
  // double charge = 1.602176487e-19;
  // double elcMass = 9.10938215e-31;
  // double ionMass = 1.672621637e-27;
  // double epsilon0 = 8.854187817620389850536563031710750260608e-12;
  // double mu0 = 12.56637061435917295385057353311801153679e-7;
  double c = 1.0/sqrt(epsilon0*mu0);
  double vAe = c;
  // plasma density in particles per m^3
  double n = 1.0;
  // Compute magnetic field from specified parameters
  double B0 = vAe*sqrt(mu0*n*elcMass);

  double beta = 0.1;
  double T_ion = beta*(B0*B0)/(2.0*n*mu0);
  double T_elc = T_ion;

  double wpi = sqrt(charge*charge*n/(epsilon0*ionMass));
  double di = c/wpi;
  double Lx = 256.0*di;
  
  double rel_fac = ionMass*c*c/T_ion;
  double coll_fac = 1e5;
  double tau = coll_fac*6.0*sqrt(2.0*M_PI*elcMass*T_elc*M_PI*T_elc*M_PI*T_elc)*epsilon0*epsilon0/(charge*charge*charge*charge*n);
  double lambda = 1.0/mu0*(elcMass/(charge*charge*n*tau));

  double eta_par = 0.96*n*T_ion*sqrt(2.0)*tau;
  double eta_perp = 0.3*n*T_ion/(tau*charge*charge*B0*B0/(ionMass*ionMass));
  double mu_par = 1.5*eta_par - 15.0/2.0*eta_perp;
  double hartmann_num = B0*Lx/sqrt(mu0*lambda*mu_par);

  double accel = 0.01;
  double force_norm = mu0*n*accel*Lx/(B0*B0);
  struct hartmann_ctx ctx = {
    .gasGamma = gasGamma,
    .charge = charge,
    .elcMass = elcMass,
    .ionMass = ionMass,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .T_elc = T_elc,
    .T_ion = T_ion,
    .beta = beta,
    .B0 = B0,
    .n = n,
    .Lx = Lx,
    .rel_fac = rel_fac,
    .coll_fac = coll_fac,
    .tau = tau,
    .lambda = lambda,
    .eta_par = eta_par,
    .eta_perp = eta_perp,
    .hartmann_num = hartmann_num,
    .accel = accel,
    .force_norm = force_norm,
  };
  return ctx;
}

void
evalElcInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  struct hartmann_ctx *app = ctx;

  double rhoe = app->n*app->elcMass;
  double ere = app->n*app->T_elc/(app->gasGamma-1);

  fout[0] = rhoe;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = ere;  
}

void
evalIonInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  struct hartmann_ctx *app = ctx;

  double rhoi = app->n*app->ionMass;
  double eri = app->n*app->T_ion/(app->gasGamma-1);

  fout[0] = rhoi;
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = eri;   
}

void
evalFieldInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  struct hartmann_ctx *app = ctx;
  // 10 Tesla magnetic field in x direction, no other fields
  // electric field
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // magnetic field
  fout[3] = app->B0, fout[4] = 0.0; fout[5] = 0.0;

  // correction potentials
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalAppAccel(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  // Constant force
  struct hartmann_ctx *app = ctx;
  fout[0] = 0.0;
  fout[1] = app->accel;
  fout[2] = 0.0;
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
  struct hartmann_ctx ctx = create_ctx(); // context for init functions
  // electron/ion equations
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(ctx.gasGamma, app_args.use_gpu);
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(ctx.gasGamma, app_args.use_gpu);

  printf("Parallel viscosity = %lg\n", ctx.eta_par);
  printf("Perpendicular viscosity = %lg\n", ctx.eta_perp);
  printf("Ratio eta_par/eta_perp = %lg\n", ctx.eta_par/ctx.eta_perp);
  printf("Collision time = %lg\n", ctx.tau);
  printf("Hartmann number = %lg\n", ctx.hartmann_num);
  printf("Normalized force value = %lg\n", ctx.force_norm);
  printf("Estimated time step ratio = %lg\n", ctx.rel_fac/ctx.coll_fac);
  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = -1.0*ctx.charge, .mass = ctx.elcMass,

    .equation = elc_euler,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalElcInit,
    .app_accel_func = evalAppAccel,
    .type_brag = GKYL_BRAG_MAG_FULL,
    .coll_fac = ctx.coll_fac,  

    .bcx = { GKYL_SPECIES_NO_SLIP, GKYL_SPECIES_NO_SLIP },
  };
  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.charge, .mass = ctx.ionMass,

    .equation = ion_euler,
    .evolve = 1,
    .ctx = &ctx,
    .init = evalIonInit,
    .app_accel_func = evalAppAccel,
    .type_brag = GKYL_BRAG_MAG_FULL,
    .coll_fac = ctx.coll_fac,

    .bcx = { GKYL_SPECIES_NO_SLIP, GKYL_SPECIES_NO_SLIP },    
  };  

  // VM app
  struct gkyl_moment app_inp = {
    .name = "5m_hartmann",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx }, 
    .cells = { 256 },

    .num_periodic_dir = 0,
    .periodic_dirs = { },
    .cfl_frac = 0.001,

    .num_species = 2,
    .species = { elc, ion },

    .field = {
      .epsilon0 = ctx.epsilon0,
      .mu0 = ctx.mu0,
      
      .evolve = 1,
      .ctx = &ctx,
      .init = evalFieldInit,
      
      .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
    }
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
