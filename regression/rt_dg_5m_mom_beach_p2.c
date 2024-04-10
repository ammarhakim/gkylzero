// Propagation into a plasma wave beach test for the 5-moment equations.
// Input parameters match the initial conditions found in entry JE8 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je8/je8-plasmabeach.html), adapted from Section III. A. of the article:
// D. N. Smithe (2007), "Finite-difference time-domain simulation of fusion plasmas at radiofrequency time scales",
// Physics of Plasmas, Volume 14 (5): 056104.
// https://pubs.aip.org/aip/pop/article/14/5/056104/929539/Finite-difference-time-domain-simulation-of-fusion

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct mom_beach_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using non-normalized physical units).
  double gas_gamma; // Adiabatic index.
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double J0; // Reference current density (Amps / m^3).
  
  // Derived physical quantities (using non-normalized physical units).
  double light_speed; // Speed of light.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  double Lx100; // Domain size over 100 (x-direction).
  double x_last_edge; // Location of center of last cell.
  double cfl_frac; // CFL coefficient.
  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.

  double deltaT; // Arbitrary constant, with units of time.
  double factor; // Numerical factor for calculation of electron number density.
  double omega_drive; // Drive current angular frequency.
  double init_dt; // Initial time step guess so first step does not generate NaN
};

struct mom_beach_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using non-normalized physical units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.
  double epsilon0 = 8.854187817620389850536563031710750260608e-12; // Permittivity of free space.
  double mu0 = 12.56637061435917295385057353311801153679e-7; // Permeability of free space.
  double mass_elc = 9.10938215e-31; // Electron mass.
  double charge_elc = -1.602176487e-19; // Electron charge.

  double J0 = 1.0; // Reference current density (Amps / m^3).
  
  // Derived physical quantities (using non-normalized physical units).
  double light_speed = 1.0 / sqrt(mu0 * epsilon0); // Speed of light.

  // Simulation parameters.
  int Nx = 400; // Cell count (x-direction).
  double Lx = 1.0; // Domain size (x-direction).
  double Lx100 = Lx / 100.0; // Domain size over 100 (x-direction).
  double x_last_edge = Lx / Nx; // Location of center of last cell.
  double cfl_frac = 1.0; // CFL coefficient.
  double t_end = 5.0e-10; // Final simulation time.
  int num_frames = 1; // Number of output frames.

  double deltaT = Lx100 / light_speed; // Arbitrary constant, with units of time.
  double factor = deltaT * deltaT * charge_elc * charge_elc / (mass_elc * epsilon0); // Numerical factor for calculation of electron number density.
  double omega_drive = pi / 10.0 / deltaT; // Drive current angular frequency.
  // initial dt guess so first step does not generate NaN
  double init_dt = ((Lx/Nx)/light_speed)/(5.0);

  struct mom_beach_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .J0 = J0,
    .light_speed = light_speed,
    .Nx = Nx,
    .Lx = Lx,
    .Lx100 = Lx100,
    .x_last_edge = x_last_edge,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .deltaT = deltaT,
    .factor = factor,
    .omega_drive = omega_drive,
    .init_dt = init_dt, 
  };

  return ctx;
}

void
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct mom_beach_ctx *app = ctx;

  double gas_gamma = app -> gas_gamma;
  double epsilon0 = app -> epsilon0;
  double mass_elc = app -> mass_elc;
  double charge_elc = app -> charge_elc;

  double light_speed = app -> light_speed;

  double Lx100 = app -> Lx100;

  double factor = app ->factor;

  double omegaPdt = 25.0 * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x); // Plasma frequency profile.
  double ne = omegaPdt * omegaPdt / factor; // Electron number density.

  // Set electron mass density.
  fout[0] = mass_elc * ne;
  // Set electron momentum density.
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  // Set electron total energy density.
  fout[4] = ne * (-charge_elc) / (gas_gamma - 1.0);  
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set electric field.
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // Set magnetic field.
  fout[3] = 0.0, fout[4] = 0.0; fout[5] = 0.0;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalAppCurrent(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct mom_beach_ctx *app = ctx;

  double pi = app -> pi;

  double J0 = app -> J0;

  double light_speed = app -> light_speed;

  double Lx = app -> Lx;
  double Nx = app -> Nx;
  double Lx100 = app -> Lx100;
  double x_last_edge = app -> x_last_edge;

  double omega_drive = app -> omega_drive;

  if (x > x_last_edge) {
    // Set applied current.
    fout[1] = -J0 * sin(omega_drive * t);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct mom_beach_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  // Electron equations.
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);

  struct gkyl_vlasov_fluid_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .equation = elc_euler,
    .init = evalElcInit,
    .ctx = &ctx,
    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };
  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,
    .limit_em = true, 

    .init = evalFieldInit,
    .ctx = &ctx,
    .app_current = evalAppCurrent,
    .app_current_ctx = &ctx,
    .app_current_evolve = true,
    .bcx = { GKYL_FIELD_COPY, GKYL_FIELD_COPY }, 
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NX };
  int dim = sizeof(cells) / sizeof(cells[0]);
  struct gkyl_range global_r;
  gkyl_create_global_range(dim, cells, &global_r);

  // Create decomposition.
  int cuts[dim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < dim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < dim; d++) {
    cuts[d] = 1;
  }
#endif

  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(dim, cuts, &global_r);

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp,
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < dim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // VM app
  struct gkyl_vm vm = {
    .name = "dg_5m_mom_beach_p2",

    .cdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx }, 
    .cells = { NX },
    .poly_order = 2,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 0,
    .species = {},

    .num_fluid_species = 1,
    .fluid_species = { elc },

    .field = field,

    .use_gpu = app_args.use_gpu,

    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp -> ranges[my_rank],
      .comm = comm
    }
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double t_curr = 0.0, t_end = ctx.t_end;
  double dt = ctx.init_dt;

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, t_curr);
  
  gkyl_vlasov_app_write(app, t_curr, 0);

  long step = 1, num_steps = app_args.num_steps;
  while ((t_curr < t_end) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", t_curr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      fprintf(stderr, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;
    step += 1;
  }

  gkyl_vlasov_app_write(app, t_curr, 1);
  gkyl_vlasov_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

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

  // Free resources after simulation completion.
  gkyl_wv_eqn_release(elc_euler);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_vlasov_app_release(app);
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
