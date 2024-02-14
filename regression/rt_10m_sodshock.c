#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_ten_moment.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct sodshock_ctx
{
  // Physical constants (using normalized code units).
  double rhol; // Left fluid mass density.
  double ul; // Left fluid velocity.
  double pl; // Left fluid pressure.

  double rhor; // Right fluid mass density.
  double ur; // Right fluid velocity.
  double pr; // Right fluid pressure.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  double k0; // Closure parameter.
  double cfl_frac; // CFL coefficient.
  double t_end; // Final simulation time.
};

struct sodshock_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double rhol = 3.0; // Left fluid mass density.
  double ul = 0.0; // Left fluid velocity.
  double pl = 3.0; // Left fluid pressure.

  double rhor = 1.0; // Right fluid mass density.
  double ur = 0.0; // Right fluid velocity.
  double pr = 1.0; // Right fluid pressure.

  // Simulation parameters.
  int Nx = 512; // Cell count (x-direction).
  long Lx = 1.0; // Domain size (x-direction).
  double k0 = 0.0; // Closure parameter.
  double cfl_frac = 0.9; // CFL coefficient.
  double t_end = 0.1; // Final simulation time.

  struct sodshock_ctx ctx = {
    .rhol = rhol,
    .ul = ul,
    .pl = pl,
    .rhor = rhor,
    .ur = ur,
    .pr = pr,
    .Nx = Nx,
    .Lx = Lx,
    .k0 = k0,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
  };

  return ctx;
}

void
eval10mInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct sodshock_ctx *app = ctx;

  double rhol = app -> rhol;
  double ul = app -> ul;
  double pl = app -> pl;

  double rhor = app -> rhor;
  double ur = app -> ur;
  double pr = app -> pr;

  double rho = 0.0;
  double u = 0.0;
  double p = 0.0;

  if (x < 0.75)
  {
    rho = rhol; // Fluid momentum density (left).
    u = ul; // Fluid velocity (left).
    p = pl; // Fluid pressure (left).
  }
  else
  {
    rho = rhor; // Fluid momentum density (right).
    u = ur; // Fluid velocity (right).
    p = pr; // Fluid pressure (right).
  }

  // Set fluid mass density.
  fout[0] = rho;
  // Set fluid momentum density.
  fout[1] = rho * u; fout[2] = 0.0; fout[3] = 0.0;
  // Set fluid pressure tensor.
  fout[4] = p + 0.5 * rho * u * u; fout[5] = 0.0; fout[6] = 0.0;
  fout[7] = p; fout[8] = 0.0; fout[9] = p;
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

  struct sodshock_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  // Fluid equations.
  struct gkyl_wv_eqn *ten_moment = gkyl_wv_ten_moment_new(ctx.k0);

  struct gkyl_moment_species fluid = {
    .name = "10m",
    .equation = ten_moment,
    .evolve = true,
    .init = eval10mInit,
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
    .name = "10m_sodshock",

    .ndim = 1,
    .lower = { 0.25 },
    .upper = { 0.25 + ctx.Lx }, 
    .cells = { NX },

    .cfl_frac = ctx.cfl_frac,

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
  gkyl_wv_eqn_release(ten_moment);
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
