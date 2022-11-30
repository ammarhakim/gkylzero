#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct pkpm_par_firehose_ctx {
  double epsilon0;
  double mu0;
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te_Ti; // electron to ion temperature ratio
  double n0;
  double vAe;
  double B0;
  double beta;
  double vtElc;
  double vtIon;
  double elcTemp;
  double ionTempPar;
  double ionTempPerp;
  double nuElc;
  double nuIon;
  double L;
  double tend;
  bool use_gpu;
};

static inline double sq(double x) { return x*x; }

void
evalDistFuncElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_par_firehose_ctx *app = ctx;
  double x = xn[0], v = xn[1];

  double vt = app->vtElc;
  double n0 = app->n0;

  fout[0] = n0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
}

void
evalDistFuncIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_par_firehose_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double mi = app->massIon;
  double Lx = app->L;
  double T_par = app->ionTempPar;
  double n0 = app->n0;
  double vt = sqrt(T_par/mi);

  fout[0] = n0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_par_firehose_ctx *app = ctx;
  double x = xn[0];
  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double Lx = app->L;
  double elcTemp = app->elcTemp;
  double n0 = app->n0;

  // flow velocity initialized as a tanh function to avoid issues with recovery of a step function
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = n0*me*elcTemp;
}

void
evalFluidIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_par_firehose_ctx *app = ctx;
  double x = xn[0];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double Lx = app->L;
  double T_perp = app->ionTempPerp;
  double n0 = app->n0;

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = n0*mi*T_perp;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_par_firehose_ctx *app = ctx;
  double x = xn[0];
  pcg64_random_t rng = gkyl_pcg64_init(0);

  double Bx = app->B0;
  double By = 0.0, Bz = 0.0;

  double alpha = 1.0e-6*app->B0;
  double Lx = app->L;
  double kx = 2.0*M_PI/Lx;
  for (int i=0; i<64; ++i) {
    By -= alpha*gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
    Bz -= alpha*gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
  }

  // electric field
  fout[0] = 0.0; fout[1] = 0.0; fout[2] = 0.0;
  // magnetic field
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;

}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_par_firehose_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_par_firehose_ctx *app = ctx;
  fout[0] = app->nuIon;
}

struct pkpm_par_firehose_ctx
create_ctx(void)
{
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // Permiability of free space.
  double lightSpeed = 1/sqrt(mu0*epsilon0); // Speed of light.
  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge
  double massIon = 1836.153; // ion mass
  double chargeIon = 1.0; // ion charge

  double Te_Ti = 1.0; // ratio of electron to ion temperature
  double n0 = 1.0; // initial number density

  double vAe = 0.0125;
  double B0 = vAe*sqrt(mu0*n0*massElc);
  double beta = 300.0/M_PI;

  double dbeta = 100.0; // dbeta = beta_parallel - beta_perp
  double betaPar = beta + 2.*dbeta/3.; // parallel proton plasma beta = 2 mu0 ni T_paralleli/B0^2
  double betaPerp = beta - dbeta/3.; // perp proton plasma beta = 2 mu0 ni T_perpi/B0^2  

  double vtElc = vAe*sqrt(beta);
  double vtIon = vtElc/sqrt(massIon*Te_Ti/massElc);
  double elcTemp = vtElc*vtElc/2.0;

  double ionTempPar = vAe*vAe*(betaPar*massElc/2);
  double ionTempPerp = vAe*vAe*(betaPerp*massElc/2);

  // ion cyclotron frequency and gyroradius
  double omegaCi = chargeIon*B0/massIon;
  double wpi = sqrt(n0*chargeIon*chargeIon/(epsilon0*massIon));
  double di = lightSpeed/wpi;

  // collision frequencies
  double nuElc = 0.1*omegaCi;
  double nuIon = 0.1*omegaCi/sqrt(massIon);

  // domain size and simulation time
  double L = 100.0*M_PI*di;
  double tend = 0.1/omegaCi;

  struct pkpm_par_firehose_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .chargeElc = chargeElc,
    .massElc = massElc,
    .chargeIon = chargeIon,
    .massIon = massIon,
    .Te_Ti = Te_Ti,
    .n0 = n0,
    .vAe = vAe,
    .B0 = B0,
    .beta = beta,
    .vtElc = vtElc,
    .vtIon = vtIon,
    .elcTemp = elcTemp,
    .ionTempPar = ionTempPar,
    .ionTempPerp = ionTempPerp,
    .nuElc = nuElc,
    .nuIon = nuIon,
    .L = L,
    .tend = tend,
  };
  return ctx;
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
  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 8192);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 32);  

  struct pkpm_par_firehose_ctx ctx = create_ctx(); // context for init functions

  // electron Tperp                                                                                              
  struct gkyl_vlasov_fluid_species fluid_elc = {
    .name = "fluid_elc",
    .num_eqn = 4,
    .pkpm_species = "elc",
    .ctx = &ctx,
    .init = evalFluidElc,
  };  
  
  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -6.0 * ctx.vtElc},
    .upper = { 6.0 * ctx.vtElc}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
    },    

    .num_diag_moments = 0,
  };

  // ion Tperp                                                                                              
  struct gkyl_vlasov_fluid_species fluid_ion = {
    .name = "fluid_ion",
    .num_eqn = 4,
    .pkpm_species = "ion",
    .ctx = &ctx,
    .init = evalFluidIon,
  };  
  
  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -6.0 * ctx.vtIon},
    .upper = { 6.0 * ctx.vtIon}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncIon,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
    },    

    .num_diag_moments = 0,
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "pkpm_par_firehose_p1",

    .cdim = 1, .vdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.L },
    .cells = { NX },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, ion },
    .num_fluid_species = 2,
    .fluid_species = { fluid_elc, fluid_ion },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.tend;
  double dt = tend-tcurr;
  int nframe = 100;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_vlasov_app_calc_field_energy(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);

    gkyl_vlasov_app_calc_field_energy(app, tcurr);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_vlasov_app_write_field_energy(app);
  gkyl_vlasov_app_stat_write(app);

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
