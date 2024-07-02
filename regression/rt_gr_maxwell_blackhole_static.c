// 2D plane electromagnetic wave propagating onto a (Schwarzschild) black hole, for the general relativistic Maxwell equations.
// Input parameters describe a plane electromagnetic wave being lensed by a non-rotating black hole.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_gr_maxwell.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct maxwell_blackhole_static_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double light_speed; // Speed of light.
  double e_fact; // Factor of speed of light for use in electric field correction.
  double b_fact; // Factor of speed of light for use in magnetic field correction.

  double k_wave_x; // Wave number (x-direction).
  double k_wave_y; // Wave number (y-direction).
  double E0; // Reference electromagnetic field strength.

  // Derived physical quantities (using normalized code units).
  double k_norm; // Wave number normalization.

  double k_xn; // Normalized wave number (x-direction).
  double k_yn; // Normalized wave number (y-direction).

  // Spacetime parameters (using geometric units).
  double mass; // Mass of the black hole.
  double spin; // Spin of the black hole.

  double pos_x; // Position of the black hole (x-direction).
  double pos_y; // Position of the black hole (y-direction).
  double pos_z; // Position of the black hole (z-direction).

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime;

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  int Ny; // Cell count (y-direction).
  double Lx; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct maxwell_blackhole_static_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double light_speed = 1.0; // Speed of light.
  double e_fact = 0.0; // Factor of speed of light for use in electric field correction.
  double b_fact = 1.0; // Factor of speed of light for use in magnetic field correction.

  double k_wave_x = 2.0; // Wave number (x-direction).
  double k_wave_y = 2.0; // Wave number (y-direction).
  double E0 = 1.0 / sqrt(3.0); // Reference electromagnetic field strength.

  // Derived physical quantities (using normalized code units).
  double k_norm = sqrt((k_wave_x * k_wave_x) + (k_wave_y * k_wave_y)); // Wave number normalization.

  double k_xn = k_wave_x / k_norm; // Normalized wave number (x-direction).
  double k_yn = k_wave_y / k_norm; // Normalized wave number (y-direction).

  // Spacetime parameters (using geometric units).
  double mass = 0.05; // Mass of the black hole.
  double spin = 0.0; // Spin of the black hole.

  double pos_x = 0.5; // Position of the black hole (x-direction).
  double pos_y = 0.5; // Position of the black hole (y-direction).
  double pos_z = 0.0; // Position of the black hole (z-direction).

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, mass, spin, pos_x, pos_y, pos_z);

  // Simulation parameters.
  int Nx = 256; // Cell count (x-direction).
  int Ny = 256; // Cell count (y-direction).
  double Lx = 1.0; // Domain size (x-direction).
  double Ly = 1.0; // Domain size (y-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 2.0; // Final simulation time.
  int num_frames = 100; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct maxwell_blackhole_static_ctx ctx = {
    .pi = pi,
    .light_speed = light_speed,
    .e_fact = e_fact,
    .b_fact = b_fact,
    .k_wave_x = k_wave_x,
    .k_wave_y = k_wave_y,
    .E0 = E0,
    .k_norm = k_norm,
    .k_xn = k_xn,
    .k_yn = k_yn,
    .mass = mass,
    .spin = spin,
    .pos_x = pos_x,
    .pos_y = pos_y,
    .pos_z = pos_z,
    .spacetime = spacetime,
    .Nx = Nx,
    .Ny = Ny,
    .Lx = Lx,
    .Ly = Ly,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .dt_failure_tol = dt_failure_tol,
  };

  return ctx;
}

void
evalGRMaxwellInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct maxwell_blackhole_static_ctx *app = ctx;

  double pi = app->pi;

  double k_wave_x = app->k_wave_x;
  double k_wave_y = app->k_wave_y;
  double E0 = app->E0;

  double k_xn = app->k_xn;
  double k_yn = app->k_yn;

  struct gkyl_gr_spacetime *spacetime = app->spacetime;
  
  double Lx = app->Lx;

  double phi = (2.0 * pi / Lx) * ((k_wave_x * x) + (k_wave_y * y));
  
  double spatial_det, lapse;
  bool in_excision_region;

  spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
  spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

  double Ex = -E0 * cos(phi); // Electric field (x-direction).
  double Ey = E0 * cos(phi); // Electric field (y-direction).
  double Ez = E0 * cos(phi); // Electric field (z-direction).

  double Bx = (E0 * cos(phi) * 2.0 * pi / Lx) * k_yn; // Magnetic field (x-direction).
  double By = -(E0 * cos(phi) * 2.0 * pi / Lx) * k_xn; // Magnetic field (y-direction).
  double Bz = (E0 * cos(phi) * 2.0 * pi / Lx) * (-k_xn - k_yn); // Magnetic field (z-direction).

  // Set electric field.
  fout[0] = Ex; fout[1] = Ey; fout[2] = Ez;
  // Set magnetic field.
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;

  // Set spatial metric determinant.
  fout[8] = spatial_det;
  // Set lapse gauge variable.
  fout[9] = lapse;

  // Set excision boundary conditions.
  if (in_excision_region) {
    for (int i = 0; i < 10; i++) {
      fout[i] = 0.0;
    }

    fout[0] = -5.0;

    fout[10] = -1.0;
  }
  else {
    fout[10] = 1.0;
  }
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_moment_app_write(app, t_curr, frame);
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

  struct maxwell_blackhole_static_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  // Field equations.
  struct gkyl_wv_eqn *gr_maxwell = gkyl_wv_gr_maxwell_new(ctx.light_speed, ctx.e_fact, ctx.b_fact, ctx.spacetime, app_args.use_gpu);

  struct gkyl_moment_species field = {
    .name = "gr_maxwell",
    .equation = gr_maxwell,
    .evolve = true,
    .init = evalGRMaxwellInit,
    .force_low_order_flux = true, // Use Lax fluxes.
    .ctx = &ctx,
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NX, NY };
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

  // Moment app.
  struct gkyl_moment app_inp = {
    .name = "gr_maxwell_blackhole_static",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { ctx.Lx, ctx.Ly },
    .cells = { NX, NY },

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },

    .scheme_type = GKYL_MOMENT_WAVE_PROP,
    .mp_recon = app_args.mp_recon,

    .cfl_frac = ctx.cfl_frac,

    .num_species = 1,
    .species = { field },

    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp->ranges[my_rank],
      .comm = comm
    }
  };

  // Create app object.
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // Initialize simulation.
  gkyl_moment_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr, false);

  // Compute estimate of maximum stable time-step.
  double dt = gkyl_moment_app_max_dt(app);

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
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

      gkyl_moment_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_moment_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_moment_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_moment_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_moment_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  write_data(&io_trig, app, t_curr, false);
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
  gkyl_wv_eqn_release(gr_maxwell);
  gkyl_gr_spacetime_release(ctx.spacetime);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);  
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
