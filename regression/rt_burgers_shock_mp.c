// Sinusoidal initial data test, with monotonicity-preserving reconstruction, for the inviscid Burgers' equation.
// Input parameters correspond to a sine wave, equivalent to the initial conditions in Section 4, Example 2, but with periodicity 1 rather than 2, from the article:
// X-d. Liu and S. Osher (1996), "Nonoscillatory High Order Accurate Self-Similar Maximum Principle Satisfying Shock Capturing Schemes I",
// SIAM Journal on Numerical Analysis, Volume 33 (2): 760-779.
// https://epubs.siam.org/doi/10.1137/0733038

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_burgers.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct burgers_shock_mp_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  double cfl_frac; // CFL coefficient.
  double t_end; // Final simulation time.
};

struct burgers_shock_mp_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = 3.141592653589793238462643383279502884;

  // Simulation parameters.
  int Nx = 128; // Cell count (x-direction).
  double Lx = 1.0; // Domain size (x-direction).
  double cfl_frac = 0.9; // CFL coefficient.
  double t_end = 0.4; // Final simulation time.

  struct burgers_shock_mp_ctx ctx = {
    .pi = pi,
    .Nx = Nx,
    .Lx = Lx,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
  };

  return ctx;
}

void
evalInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct burgers_shock_mp_ctx *app = ctx;

  double pi = app -> pi;

  // Set advected quantity.
  fout[0] = sin(2.0 * pi * x);
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
  {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem)
  {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct burgers_shock_mp_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  // Inviscid Burgers' equation.
  struct gkyl_wv_eqn *burgers = gkyl_wv_burgers_new();

  struct gkyl_moment_species fluid = {
    .name = "burgers",
    .equation = burgers,
    .evolve = true,
    .init = evalInit,
    .ctx = &ctx,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
  {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NX };
  struct gkyl_range globalr;
  gkyl_create_global_range(1, cells, &globalr);

  // Create decomposition.
  int cuts[] = { 1 };
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
  {
    cuts[0] = app_args.cuts[0];
  }
#endif

  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(1, cuts, &globalr);

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
  {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp)
      {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  }
  else
  {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp)
      {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gklyl_null_comm_inew( &(struct gkyl_null_comp_inp)
    {
      .decomp = decomp,
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = cuts[0];
  if (ncuts != comm_size)
  {
    if (my_rank == 0)
    {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Moment app.
  struct gkyl_moment app_inp = {
    .name = "burgers_shock_mp",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx }, 
    .cells = { NX },

    .cfl_frac = ctx.cfl_frac,

    .scheme_type = GKYL_MOMENT_MP,
    .mp_recon = app_args.mp_recon,

    .num_species = 1,
    .species = { fluid },

    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp -> ranges[my_rank],
      .comm = comm
    }
  };

  // Create app object.
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  gkyl_moment_app_apply_ic(app, t_curr);
  gkyl_moment_app_write(app, t_curr, 0);

  // Compute estimate of maximum stable time-step.
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps))
  {
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success)
    {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    step += 1;
  }

  gkyl_moment_app_write(app, t_curr, 1);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  // Free resources after simulation completion.
  gkyl_wv_eqn_release(burgers);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);  
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
  {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
