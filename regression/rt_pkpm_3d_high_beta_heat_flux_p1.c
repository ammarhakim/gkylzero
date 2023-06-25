#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct pkpm_heat_flux_ctx {
  double epsilon0;
  double mu0;
  double chargeElc; // electron charge
  double massElc; // electron mass
  double n0;
  double vAe;
  double B0;
  double beta;
  double elcTemp;
  double vtElc;
  double tempFac;
  double nuElc;
  double Lx;
  double Ly;
  double Lz;
  double tend;
  double min_dt;
  bool use_gpu;
};

static inline double
maxwellian(double n, double v, double temp)
{
  double v2 = v*v;
  // assumes mass = 1.0
  return n/sqrt(2*M_PI*temp)*exp(-v2/(2*temp));
}

void
evalDistFuncElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_heat_flux_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], z = xn[2], vx = xn[3];

  double Lx = app->Lx;
  double elcTemp = app->elcTemp;
  double tempFac = app->tempFac;

  // linear temperature gradient
  double temp = elcTemp*(1.0 - x/Lx*tempFac);
  double fv = maxwellian(app->n0, vx, temp);
    
  fout[0] = fv;
  fout[1] = temp*fv;
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_heat_flux_ctx *app = ctx;
  
  double x = xn[0], y = xn[1];

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_heat_flux_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  pcg64_random_t rng = gkyl_pcg64_init(0);

  double B_x = app->B0;
  double B_y = 0.0;
  double B_z = 0.0;

  double alpha = 1.0e-6*app->B0;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double Lz = app->Lz;
  double kx = 2.0*M_PI/Lx;
  double ky = 2.0*M_PI/Ly;
  double kz = 2.0*M_PI/Lz;
  // Modes depend on x & z for By, x and y for Bz
  for (int i=0; i<16; ++i) {
    for (int j=0; j<16; ++j) {
      B_y -= alpha*gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng))*sin(j*kz*z + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
      B_z -= alpha*gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng))*sin(j*ky*y + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
    }
  }
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = B_x; fout[4] = B_y; fout[5] = B_z;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_heat_flux_ctx *app = ctx;
  double x = xn[0];
  double Lx = app->Lx;
  double Ly = app->Ly/4.0; // Ly sets the length scale of the collisionality profile
  fout[0] = app->nuElc*(tanh((x - Lx)/Ly) - tanh(x/Ly) + 2.0);
}

struct pkpm_heat_flux_ctx
create_ctx(void)
{
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // pemiability of free space

  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge

  double n0 = 1.0; // initial number density

  double beta = 64.0;
  double elcTemp = 0.02;
  double vtElc = sqrt(2.0*elcTemp/massElc);
  double tempFac = 0.75; // Temperature is 4 times colder at right wall
  double vAe = vtElc/sqrt(beta);
  double B0 = vAe*sqrt(mu0*n0*massElc);

  // ion cyclotron frequency and gyroradius
  double omegaCe = fabs(chargeElc)*B0/massElc;
  double rhoe = vtElc/omegaCe;

  // collision frequencies
  double nuElc = omegaCe/5.0;

  // domain size and simulation time
  double Lx = 256.0*rhoe;
  double Ly = 32.0*rhoe;
  double Lz = Ly;
  double tend = 2000.0/omegaCe;
  
  struct pkpm_heat_flux_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .chargeElc = chargeElc,
    .massElc = massElc,
    .n0 = n0,
    .vAe = vAe,
    .B0 = B0,
    .beta = beta,
    .elcTemp = elcTemp,
    .vtElc = vtElc,
    .tempFac = tempFac,
    .nuElc = nuElc,
    .Lx = Lx,
    .Ly = Ly,
    .Lz = Lz,
    .tend = tend,
    .min_dt = 1.0e-4, 
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 256);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 32);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], 32);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 16);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
     
  struct pkpm_heat_flux_ctx ctx = create_ctx(); // context for init functions

  // electron momentum 
  struct gkyl_vlasov_fluid_species fluid_elc = {
    .name = "fluid_elc",
    .num_eqn = 3,
    .pkpm_species = "elc",
    .ctx = &ctx,
    .init = evalFluidElc,
    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
    .diffusion = {.D = 1.0e-2, .order=4},
  };  
  
  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -4.0 * ctx.vtElc},
    .upper = { 4.0 * ctx.vtElc}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
    },    

    .bcx = { GKYL_SPECIES_FIXED_FUNC, GKYL_SPECIES_FIXED_FUNC },

    .num_diag_moments = 0,
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,
    .bcx = { GKYL_FIELD_RESERVOIR, GKYL_FIELD_RESERVOIR }
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "pkpm_3d_high_beta_heat_flux_p1",

    .cdim = 3, .vdim = 1,
    .lower = { 0.0, 0.0, 0.0 },
    .upper = { ctx.Lx, ctx.Ly, ctx.Lz },
    .cells = { NX, NY, NZ },
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .cfl_frac = 0.5,
    
    .num_periodic_dir = 2,
    .periodic_dirs = { 1, 2 },

    .num_species = 1,
    .species = { elc },
    .num_fluid_species = 1,
    .fluid_species = { fluid_elc },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.tend;
  double dt = tend-tcurr;
  int nframe = 200;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    if (status.dt_actual < ctx.min_dt) {
      printf("** Time step crashing! Aborting simulation and writing out last output ....\n");
      gkyl_vlasov_app_write(app, tcurr, 1000);
      gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 1000);
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }

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
  printf("Species collisions took %g secs\n", stat.species_coll_mom_tm);
  printf("Species collisions took %g secs\n", stat.species_coll_tm);
  printf("Field RHS calc took %g secs\n", stat.field_rhs_tm);
  printf("Current evaluation and accumulate took %g secs\n", stat.current_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
