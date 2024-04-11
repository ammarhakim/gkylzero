// Generalized Brio-Wu Riemann problem for the 5-moment equations.
// Input parameters match the initial conditions found in entry JE4 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je4/je4-twofluid-shock.html), adapted from Section 7.1 of the article:
// A. Hakim, J. Loverich and U. Shumlak (2006), "A high resolution wave propagation scheme for ideal Two-Fluid plasma equations",
// Journal of Computational Physics, Volume 219 (1): 418-442.
// https://www.sciencedirect.com/science/article/pii/S0021999106001707

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

struct riem_ctx
{
  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double rhol_ion; // Left ion mass density.
  double rhor_ion; // Right ion mass density;
  double pl; // Left electron/ion pressure.
  double pr; // Right electron/ion pressure.

  double Bx; // Total magnetic field (x-direction).
  double Bzl; // Left total magneic field (z-direction).
  double Bzr; // Right total magnetic field (z-direction).

  // Derived physical quantities (using normalized code units).
  double rhol_elc; // Left electron mass density.
  double rhor_elc; // Right electron mass density.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  double init_dt; // Initial time step guess so first step does not generate NaN
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct riem_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 1.0; // Proton mass.
  double charge_ion = 1.0; // Proton charge.
  double mass_elc = 1.0 / 1836.2; // Electron mass.
  double charge_elc = -1.0; // Electron charge.

  double rhol_ion = 1.0; // Left ion mass density.
  double rhor_ion = 0.125; // Right ion mass density;
  double pl = 5.0e-5; // Left electron/ion pressure.
  double pr = 5.0e-6; // Right electron/ion pressure.

  double Bx = 0.75e-2; // Total magnetic field (x-direction).
  double Bzl = 1.0e-2; // Left total magneic field (z-direction).
  double Bzr = -1.0e-2; // Right total magnetic field (z-direction).

  // Derived physical quantities (using normalized code units).
  double rhol_elc = rhol_ion * mass_elc / mass_ion; // Left electron mass density.
  double rhor_elc = rhor_ion * mass_elc / mass_ion; // Right electron mass density.

  // Simulation parameters.
  int Nx = 1024; // Cell count (x-direction).
  double Lx = 1.0; // Domain size (x-direction).
  double cfl_frac = 1.0; // CFL coefficient.
  double t_end = 0.05; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.  
  // initial dt guess so first step does not generate NaN (speed of light = 1.0)
  double init_dt = (Lx/Nx)/7.0;
  
  struct riem_ctx ctx = {
    .gas_gamma = gas_gamma,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .rhol_ion = rhol_ion,
    .rhor_ion = rhor_ion,
    .pl = pl,
    .pr = pr,
    .Bx = Bx,
    .Bzl = Bzl,
    .Bzr = Bzr,
    .rhol_elc = rhol_elc,
    .rhor_elc = rhor_elc,
    .Nx = Nx,
    .Lx = Lx,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .init_dt = init_dt, 
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct riem_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;

  double pl = app->pl;
  double pr = app->pr;

  double rhol_elc = app->rhol_elc;
  double rhor_elc = app->rhor_elc;

  double rho = 0.0;
  double p = 0.0;

  if (x < 0.5) {
    rho = rhol_elc; // Electron mass density (left).
    p = pl; // Electron pressure (left).
  }
  else {
    rho = rhor_elc; // Electron mass density (right).
    p = pr; // Electron pressure (right).
  }

  // Set electron mass density.
  fout[0] = rho;
  // Set electron momentum density.
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  // Set electron total energy density.
  fout[4] = p / (gas_gamma - 1.0);  
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct riem_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;

  double rhol_ion = app->rhol_ion;
  double rhor_ion = app->rhor_ion;
  double pl = app->pl;
  double pr = app->pr;

  double rho = 0.0;
  double p = 0.0;

  if (x < 0.5) {
    rho = rhol_ion; // Ion mass density (left).
    p = pl; // Ion pressure (left).
  }
  else {
    rho = rhor_ion; // Ion mass density (right).
    p = pr; // Ion pressure (right).
  }

  // Set ion mass density.
  fout[0] = rho;
  // Set ion momentum density.
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  // Set ion total energy density.
  fout[4] = p / (gas_gamma - 1.0);  
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct riem_ctx *app = ctx;

  double Bx = app->Bx;
  double Bzl = app->Bzl;
  double Bzr = app->Bzr;

  double Bz = 0.0;

  if (x < 0.5) {
    Bz = Bzl; // Total magnetic field (z-direction, left).
  }
  else {
    Bz = Bzr; // Total magnetic field (z-direction, right).
  }

  // Set electric field.
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // Set magnetic field.
  fout[3] = Bx, fout[4] = 0.0; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_vlasov_app_write(app, t_curr, frame);
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

  struct riem_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  // Electron/ion equations.
  struct gkyl_wv_eqn *elc_euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);
  struct gkyl_wv_eqn *ion_euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);

  struct gkyl_vlasov_fluid_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .equation = elc_euler,
    .init = evalElcInit,
    .ctx = &ctx,
    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };
  struct gkyl_vlasov_fluid_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .equation = ion_euler,
    .init = evalIonInit,
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
    .name = "dg_5m_riem_p3",

    .cdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx }, 
    .cells = { NX },
    .poly_order = 3,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 0,
    .species = {},

    .num_fluid_species = 2,
    .fluid_species = { elc, ion },

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

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // Initialize simulation.
  gkyl_vlasov_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr, false);

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_vlasov_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_vlasov_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_vlasov_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_vlasov_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_vlasov_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  write_data(&io_trig, app, t_curr, false);
  gkyl_vlasov_app_stat_write(app);

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
  gkyl_wv_eqn_release(ion_euler);
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
