// Supersonic flow over a blunt body test for the Euler equations.
// Input parameters match the initial conditions found in Entry JE24 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je24/je24-euler-embedded-bc.html#supersonic-flow-over-blunt-body).
// Problem setup is similar to Section V. C. of the article:
// P. V. Tota and Z. J. Wang (2007), "Meshfree Euler Solver using local Radial Basis Functions for inviscid Compressible Flows",
// 18th AIAA Computational Fluid Dynamics Conference (25th June 2007 - 28th June 2007, Miami, Florida).
// https://arc.aiaa.org/doi/10.2514/6.2007-4581

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct bump_in_channel_ctx
{
  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.

  double rhol; // Left fluid mass density.
  double pl; // Left fluid pressure.

  double rhor; // Right fluid mass density.
  double pr; // Right fluid pressure.

  // Derived physical quantities (using normalized code units).
  double cs; // Fluid sound speed;

  double ul; // Left fluid velocity.
  double ur; // Right fluid velocity;

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  int Ny; // Cell count (y-direction).
  double cfl_frac; // CFL coefficient.
  double t_end; // Final simulation time.

  double inlet_xlen; // Length of inlet (x-direction).
  double bump_xlen; // Length of bump (x-direction).
  double height; // Height of channel.
  double R; // Radius of bump.

  double half_xlen; // Mid-point of bump (x-direction).

  double Lx; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
};

struct bump_in_channel_ctx
create_ctx(void)
{
    // Physical constants (using normalized code units).
  double gas_gamma = 1.4; // Adiabatic index.

  double rhol = 1.0; // Left fluid mass density.
  double pl = 1.0; // Left fluid pressure.

  double rhor = 1e-5; // Right fluid mass density.
  double pr = 1e-5; // Right fluid pressure.

  // Derived physical quantities (using normalized code units).
  double cs = sqrt(gas_gamma); // Fluid sound speed;

  double ul = 2.0 * cs; // Left fluid velocity;
  double ur = 0.0; // Right fluid velocity;

  // Simulation parameters.
  int Nx = 150; // Cell count (x-direction).
  int Ny = 75; // Cell count (y-direction).
  double cfl_frac = 0.9; // CFL coefficient.
  double t_end = 10.0; // Final simulation time.

  double inlet_xlen = 1.0; // Length of inlet (x-direction).
  double bump_xlen = 2.0; // Length of bump (x-direction).
  double height = 2.0; // Height of channel.
  double R = 5.0; // Radius of bump.

  double half_xlen = inlet_xlen + 0.5 * bump_xlen; // Mid-point of bump (x-direction).

  double Lx = 2.0 * half_xlen; // Domain size (x-direction).
  double Ly = height; // Domain size (y-direction).

  struct bump_in_channel_ctx ctx = {
    .gas_gamma = gas_gamma,
    .rhol = rhol,
    .pl = pl,
    .rhor = rhor,
    .pr = pr,
    .cs = cs,
    .ul = ul,
    .ur = ur,
    .Nx = Nx,
    .Ny = Ny,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .inlet_xlen = inlet_xlen,
    .bump_xlen = bump_xlen,
    .height = height,
    .R = R,
    .half_xlen = half_xlen,
    .Lx = Lx,
    .Ly = Ly,
  };

  return ctx;
}

void
evalEulerInit(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = zc[0];
  struct bump_in_channel_ctx *app = ctx;

  double gas_gamma = app -> gas_gamma;

  double rhol = app -> rhol;
  double ul = app -> ul;
  double pl = app -> pl;

  double rhor = app -> rhor;
  double ur = app -> ur;
  double pr = app -> pr;

  double half_xlen = app -> half_xlen;

  double rho = 0.0;
  double u = 0.0;
  double p = 0.0;

  if (x < -half_xlen * (1.0 - 0.25))
  {
    rho = rhol; // Fluid mass density (left).
    u = ul; // Fluid velocity (left).
    p = pl; // Fluid pressure (left).
  }
  else
  {
    rho = rhor; // Fluid mass density (right).
    u = ur; // Fluid velocity (right).
    p = pr; // Fluid pressure (right).
  }

  // Set fluid mass density.
  fout[0] = rho;
  // Set fluid momenutm density.
  fout[1] = rho * u; fout[2] = 0.0; fout[3] = 0.0;
  // Set fluid total energy density.
  fout[4] = p / (gas_gamma - 1.0) + 0.5 * rho * u * u;
}

void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  double x = zc[0], y = zc[1];
  struct bump_in_channel_ctx *app = ctx;
  
  double bump_xlen = app -> bump_xlen;
  double height = app -> height;
  double R = app -> R;
  
  double zeta_min = sqrt((R * R) - ((0.5 * bump_xlen) * (0.5 * bump_xlen)));
  
  //Set physical coordinates (x, y) from computational coordinates (x, y).
  xp[0] = x;
  xp[1] = y;

  if (fabs(x) < 0.5 * bump_xlen)
  {
    double eta = x;
    double zeta = sqrt((R * R) - (eta * eta));
    double yb = zeta - zeta_min;

    xp[1] = (height - yb) * y / height + yb;
  }
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

  struct bump_in_channel_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  // Fluid equations.
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(ctx.gas_gamma, app_args.use_gpu);

  struct gkyl_moment_species fluid = {
    .name = "euler",
    .equation = euler,
    .evolve = true,
    .init = evalEulerInit,
    .ctx = &ctx,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_COPY },
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
  {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NX, NY };
  struct gkyl_range globalr;
  gkyl_create_global_range(2, cells, &globalr);

  // Create decomposition.
  int cuts[] = { 1, 1 };
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
  {
    cuts[0] = app_args.cuts[0];
    cuts[1] = app_args.cuts[1];
  }
#endif

  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &globalr);

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
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp)
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

  int ncuts = cuts[0] * cuts[1];
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
    .name = "euler_bump_in_channel",

    .ndim = 2,
    .lower = { -0.5 * ctx.Lx, 0.0 },
    .upper = { 0.5 * ctx.Lx, ctx.Ly },
    .cells = { NX, NY },

    .mapc2p = mapc2p,
    .c2p_ctx = &ctx,
    
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
  gkyl_wv_eqn_release(euler);
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