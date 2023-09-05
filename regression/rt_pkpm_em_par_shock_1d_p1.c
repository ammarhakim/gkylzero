#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct pkpm_em_par_shock_ctx {
  double epsilon0;
  double mu0;
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double n0;
  double B0;
  double noise_amp;
  int k_init;
  int k_final;
  double T_e;
  double T_i;
  double uShock;
  double vtElc;
  double vtIon;
  double nuElc;
  double nuIon;
  double Lx;
  double Ly;
  double tend;
  double min_dt;
  bool use_gpu;
};

static inline void
noise_init(double noise_amp, int k_init, int k_final, double Lx,  double x, double noise[2])
{
  pcg64_random_t rng = gkyl_pcg64_init(0);
  double kx = 2.0*M_PI/Lx;
  double rand_amp, rand_phase_x;
  for (int i = k_init; i < k_final; ++i) {
    for (int j = k_init; j < k_final; ++j) {
      rand_amp = gkyl_pcg64_rand_double(&rng);
      rand_phase_x = gkyl_pcg64_rand_double(&rng);
      noise[0] -= noise_amp*rand_amp*i*kx*cos(i*kx*x + 2.0*M_PI*rand_phase_x);
      noise[1] += noise_amp*rand_amp*(i*i*kx*kx*sin(i*kx*x + 2.0*M_PI*rand_phase_x));
    }
  }
}

static inline double
maxwellian(double n, double v, double temp, double mass)
{
  double v2 = v*v;
  return n/sqrt(2.0*M_PI*temp/mass)*exp(-v2/(2.0*temp/mass));
}

void
evalDistFuncElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_em_par_shock_ctx *app = ctx;
  double x = xn[0], vx = xn[1];

  double me = app->massElc;
  double mi = app->massIon;
  double T_e = app->T_e;
  double T_i = app->T_i;
  double Lx = app->Lx;
  double n0 = app->n0;
  
  double fv = maxwellian(n0, vx, T_e, me);
    
  fout[0] = fv;
  fout[1] = T_e/me*fv;
}

void
evalDistFuncIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_em_par_shock_ctx *app = ctx;
  double x = xn[0], vx = xn[1];

  double me = app->massElc;
  double mi = app->massIon;
  double T_e = app->T_e;
  double T_i = app->T_i;
  double Lx = app->Lx;
  double n0 = app->n0;
  
  double fv = maxwellian(n0, vx, T_i, mi);
    
  fout[0] = fv;
  fout[1] = T_i/mi*fv;
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_em_par_shock_ctx *app = ctx;
  double x = xn[0];
  double vdrift = app->uShock, n0 = app->n0;
  double mass = app->massElc;
  double charge = app->chargeElc;
  double Lx = app->Lx;
  double B0 = app->B0;
  double noise_amp = app->noise_amp;
  int k_init = app->k_init;
  int k_final = app->k_final;

  double noise[2] = {0.0};
  noise_init(noise_amp, k_init, k_final, Lx, x, noise);

  fout[0] = -n0*mass*vdrift;
  fout[1] = 0.0;
  fout[2] = mass*noise[1]/charge; // initial noise in Jz_elc to excite instabilities
}

void
evalFluidIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_em_par_shock_ctx *app = ctx;
  double x = xn[0];
  double vdrift = app->uShock, n0 = app->n0;
  double mass = app->massIon;
  fout[0] = -n0*mass*vdrift;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_em_par_shock_ctx *app = ctx;
  double x = xn[0];
  double Lx = app->Lx;
  double B0 = app->B0;
  double noise_amp = app->noise_amp;
  int k_init = app->k_init;
  int k_final = app->k_final;

  double noise[2] = {0.0};
  noise_init(noise_amp, k_init, k_final, Lx, x, noise);
  // corresponding noise to By
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = noise[0]; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}


void
evalExtEmFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_em_par_shock_ctx *app = ctx;
  double x = xn[0];
  double B_x = app->B0;
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = B_x; fout[4] = 0.0; fout[5] = 0.0;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_em_par_shock_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_em_par_shock_ctx *app = ctx;
  fout[0] = app->nuIon;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_vlasov_app_write(app, tcurr, iot->curr-1);
    gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, iot->curr-1);
  }
}

struct pkpm_em_par_shock_ctx
create_ctx(void)
{
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // pemiability of free space

  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge
  double massIon = 1836.153; // ion mass
  double chargeIon = 1.0; // ion charge

  double n0 = 1.0; // initial number density
  double vAe = 1.0/5.0;
  double vAi = vAe/sqrt(massIon);
  double beta = 1.0;

  double B0 = vAe*sqrt(mu0*n0*massElc);
  double T_e = beta*B0*B0/(2.0*mu0*n0);
  double T_i = T_e;
  double vtElc = sqrt(2.0*T_e/massElc);
  double vtIon = sqrt(2.0*T_i/massIon);
  double uShock = 2.0*vAi; 

  // ion cyclotron frequency and gyroradius
  double omegaCi = chargeIon*B0/massIon;
  double di = vAi/omegaCi;
  double rhoi = vtIon/omegaCi; // rhoi ~ 484 lambdaD at real mass ratio and vtElc/c = 1/8

  // noise levels for perturbation
  double noise_amp = 0.001*B0;
  int k_init = 1;            // first wave mode to perturb with noise, 1.0 correspond to box size
  int k_final = 32;          // last wave mode to perturb with noise

  // collision frequencies
  double nuElc = 0.01*omegaCi;
  double nuIon = 0.01*omegaCi/sqrt(massIon);

  // domain size and simulation time
  double Lx = 2.0*M_PI*rhoi;
  double tend = 10.0/omegaCi;

  struct pkpm_em_par_shock_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .chargeElc = chargeElc,
    .massElc = massElc,
    .chargeIon = chargeIon,
    .massIon = massIon,
    .T_e = T_e,
    .T_i = T_i,
    .n0 = n0,
    .B0 = B0,
    .noise_amp = noise_amp,
    .k_init = k_init,
    .k_final = k_final,
    .vtElc = vtElc,
    .vtIon = vtIon,
    .uShock = uShock,
    .nuElc = nuElc,
    .nuIon = nuIon,
    .Lx = Lx,
    .tend = tend,
    .min_dt = 1.0e-2, 
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 128);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 32);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
     
  struct pkpm_em_par_shock_ctx ctx = create_ctx(); // context for init functions

  int nrank = 1; // number of processors in simulation
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
#endif  

  // create global range
  int cells[] = { NX };
  struct gkyl_range globalr;
  gkyl_create_global_range(1, cells, &globalr);
  
  // create decomposition
  int cuts[] = { 1 };
#ifdef GKYL_HAVE_MPI  
  if (app_args.use_mpi) {
    cuts[0] = app_args.cuts[0];
  }
#endif 
    
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(1, cuts, &globalr);

  // construct communcator for use in app
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  }
  else
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu        
      }
    );
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp,
      .use_gpu = app_args.use_gpu      
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_sz;
  gkyl_comm_get_size(comm, &comm_sz);

  int ncuts = cuts[0];
  if (ncuts != comm_sz) {
    if (my_rank == 0)
      fprintf(stderr, "*** Number of ranks, %d, do not match total cuts, %d!\n", comm_sz, ncuts);
    goto mpifinalize;
  }

  // electron momentum
  struct gkyl_vlasov_fluid_species fluid_elc = {
    .name = "fluid_elc",
    .num_eqn = 3,
    .pkpm_species = "elc",
    .ctx = &ctx,
    .init = evalFluidElc,
    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_COPY },
    .diffusion = {.D = 1.0e-5, .order=4},
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

    .num_diag_moments = 0,
    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_COPY },
  };

  // ion momentum
  struct gkyl_vlasov_fluid_species fluid_ion = {
    .name = "fluid_ion",
    .num_eqn = 3,
    .pkpm_species = "ion",
    .ctx = &ctx,
    .init = evalFluidIon,
    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_COPY },
    .diffusion = {.D = 1.0e-5, .order=4},
  };  
  
  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .model_id = GKYL_MODEL_PKPM,
    .pkpm_fluid_species = "fluid_ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -16.0 * ctx.vtIon},
    .upper = { 16.0 * ctx.vtIon}, 
    .cells = { VX },

    .ctx = &ctx,
    .init = evalDistFuncIon,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
    },    

    .num_diag_moments = 0,
    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_COPY },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,
    // Plasma EM field BCs are PEC on left and copy on right, external field goes into wall
    .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_COPY }, 
    .ext_em = evalExtEmFunc,
    .ext_em_ctx = &ctx,
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "pkpm_em_par_shock_1d_p1",

    .cdim = 1, .vdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, ion },
    .num_fluid_species = 2,
    .fluid_species = { fluid_elc, fluid_ion },
    .field = field,

    .use_gpu = app_args.use_gpu,

    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp->ranges[my_rank],
      .comm = comm
    }
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
  gkyl_vlasov_app_calc_integrated_L2_f(app, tcurr);
  gkyl_vlasov_app_calc_integrated_mom(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    if (step % 100 == 0) {
      gkyl_vlasov_app_calc_field_energy(app, tcurr);
      gkyl_vlasov_app_calc_integrated_L2_f(app, tcurr);
      gkyl_vlasov_app_calc_integrated_mom(app, tcurr);
    }
    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    if (status.dt_actual < ctx.min_dt) {
      gkyl_vlasov_app_cout(app, stdout, "** Time step crashing! Aborting simulation and writing out last output ....\n");
      gkyl_vlasov_app_write(app, tcurr, 1000);
      gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 1000);
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }
  gkyl_vlasov_app_calc_field_energy(app, tcurr);
  gkyl_vlasov_app_calc_integrated_L2_f(app, tcurr);
  gkyl_vlasov_app_calc_integrated_mom(app, tcurr);
  gkyl_vlasov_app_write_field_energy(app);
  gkyl_vlasov_app_write_integrated_L2_f(app);
  gkyl_vlasov_app_write_integrated_mom(app);
  gkyl_vlasov_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  gkyl_vlasov_app_cout(app, stdout, "\n");
  gkyl_vlasov_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_vlasov_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_vlasov_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_vlasov_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_vlasov_app_cout(app, stdout, "Fluid Species RHS calc took %g secs\n", stat.fluid_species_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species PKPM Vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_vlasov_app_cout(app, stdout, "EM Variables (bvar) calculation took %g secs\n", stat.field_em_vars_tm);
  gkyl_vlasov_app_cout(app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm);

  gkyl_vlasov_app_cout(app, stdout, "Species BCs took %g secs\n", stat.species_bc_tm);
  gkyl_vlasov_app_cout(app, stdout, "Fluid Species BCs took %g secs\n", stat.fluid_species_bc_tm);
  gkyl_vlasov_app_cout(app, stdout, "Field BCs took %g secs\n", stat.field_bc_tm);
  
  gkyl_vlasov_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  
  // simulation complete, free app
  gkyl_vlasov_app_release(app);

  mpifinalize:
  ;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif  
  
  return 0;
}
