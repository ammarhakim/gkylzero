#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_euler_mixture.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct fedkiw_shock_ctx
{
  // Physical constants (using normalized code units).
  double gas_gamma1; // First species adiabatic index.
  double gas_gamma2; // Second species adiabatic index.

  double rhol; // Left fluid mass density.
  double ul; // Left fluid velocity.
  double pl; // Left fluid pressure.
  double alpha1_l; // Left fluid volume fraction (first species).

  double rhoc; // Central fluid mass density.
  double uc; // Central fluid velocity.
  double pc; // Central fluid pressure.
  double alpha1_c; // Central fluid volume fraction (first species).

  double rhor; // Right fluid mass density.
  double ur; // Right fluid velocity.
  double pr; // Right fluid pressure.
  double alpha1_r; // Central fluid volume fraction (first species).

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct fedkiw_shock_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double gas_gamma1 = 1.4; // First species adiabatic index.
  double gas_gamma2 = 1.67; // Second species adiabatic index.

  double rhol = 1.3333; // Left fluid mass density.
  double ul = 0.3535 * sqrt(pow(10.0, 5.0)); // Left fluid velocity.
  double pl = 1.5 * pow(10.0, 5.0); // Left fluid pressure.
  double alpha1_l = 0.99999; // Left fluid volume fraction (first species).

  double rhoc = 1.0; // Central fluid mass density.
  double uc = 0.0; // Central fluid velocity.
  double pc = 1.0 * pow(10.0, 5.0); // Central fluid pressure.
  double alpha1_c = 0.99999; // Central fluid volume fraction (first species).

  double rhor = 0.1379; // Right fluid mass density.
  double ur = 0.0; // Right fluid velocity.
  double pr = 1.0 * pow(10.0, 5.0); // Right fluid pressure.
  double alpha1_r = 0.00001; // Right fluid volume fraction (first species).

  // Simulation parameters.
  int Nx = 2048; // Cell count (x-direction).
  double Lx = 1.0; // Domain size (x-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 0.0012; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct fedkiw_shock_ctx ctx = {
    .gas_gamma1 = gas_gamma1,
    .gas_gamma2 = gas_gamma2,
    .rhol = rhol,
    .ul = ul,
    .pl = pl,
    .alpha1_l = alpha1_l,
    .rhoc = rhoc,
    .uc = uc,
    .pc = pc,
    .alpha1_c = alpha1_c,
    .rhor = rhor,
    .ur = ur,
    .pr = pr,
    .alpha1_r = alpha1_r,
    .Nx = Nx,
    .Lx = Lx,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalEulerMixtureInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct fedkiw_shock_ctx *app = ctx;

  double gas_gamma1 = app->gas_gamma1;
  double gas_gamma2 = app->gas_gamma2;

  double rhol = app->rhol;
  double ul = app->ul;
  double pl = app->pl;
  double alpha1_l = app->alpha1_l;

  double rhoc = app->rhoc;
  double uc = app->uc;
  double pc = app->pc;
  double alpha1_c = app->alpha1_c;

  double rhor = app->rhor;
  double ur = app->ur;
  double pr = app->pr;
  double alpha1_r = app->alpha1_r;

  double rho1 = 0.0;
  double rho2 = 0.0;
  double alpha1 = 0.0;

  double vx_total = 0.0;
  double vy_total = 0.0;
  double vz_total = 0.0;
  double p_total = 0.0;

  if (x < 0.05) {
    rho1 = rhol; // First species fluid mass density (left).
    rho2 = rhor; // Second species fluid mass density (right).
    alpha1 = alpha1_l; // First species volume fraction (left).

    vx_total = ul; // Total mixture velocity (left).
    p_total = pl; // Total mixture pressure (left).
  }
  else if (x < 0.5) {
    rho1 = rhoc; // First species fluid mass density (central).
    rho2 = rhor; // Second species fluid mass density (right).
    alpha1 = alpha1_c; // First species volume fraction (central).

    vx_total = uc; // Total mixture velocity (central).
    p_total = pc; // Total mixture pressure (central).
  }
  else {
    rho1 = rhoc; // First species fluid mass density (central).
    rho2 = rhor; // Second species fluid mass density (right).
    alpha1 = alpha1_r; // First species volume fraction (right).

    vx_total = ur; // Total mixture velocity (right).
    p_total = pr; // Total mixture pressure (right).
  }
  double rho_total = (alpha1 * rho1) + ((1.0 - alpha1) * rho2); // Total mixture density.

  double momx_total = rho_total * vx_total; // Total mixture momentum density (x-direction).
  double momy_total = rho_total * vy_total; // Total mixture momentum density (y-direction).
  double momz_total = rho_total * vz_total; // Total mixture momentum density (z-direction).

  double E1 = (p_total / (gas_gamma1 - 1.0)) + (0.5 * rho1 * (vx_total * vx_total)); // First species total energy.
  double E2 = (p_total / (gas_gamma2 - 1.0)) + (0.5 * rho2 * (vx_total * vx_total)); // Second species total energy.
  double E_total = (alpha1 * E1) + ((1.0 - alpha1) * E2); // Total mixture energy.

  double vol_frac1 = rho_total * alpha1; // Mixture weighted volume fraction (first species).
  double mass_frac1 = alpha1 * rho1; // Mixture volume-weighted mass density (first species).
  double mass_frac2 = (1.0 - alpha1) * rho2; // Mixture volume-weighted mass density (second species).

  // Set fluid mixture total mass density.
  fout[0] = rho_total;
  // Set fluid mixture total momentum density.
  fout[1] = momx_total; fout[2] = momy_total; fout[3] = momz_total;
  // Set fluid mixture total energy density.
  fout[4] = E_total;
  // Set fluid mixture weighted volume fraction (first species).
  fout[5] = vol_frac1;
  // Set fluid mixture volume-weighted mass densities (first and second species).
  fout[6] = mass_frac1; fout[7] = mass_frac2;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_moment_app_write(app, t_curr, frame);
    gkyl_moment_app_write_field_energy(app);
    gkyl_moment_app_write_integrated_mom(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_moment_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_moment_app_calc_integrated_mom(app, t_curr);
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

  struct fedkiw_shock_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  // Fluid equations.
  double *gas_gamma_s = gkyl_malloc(sizeof(double[2]));
  gas_gamma_s[0] = ctx.gas_gamma1;
  gas_gamma_s[1] = ctx.gas_gamma2;
  struct gkyl_wv_eqn *euler_mixture = gkyl_wv_euler_mixture_new(2, gas_gamma_s, app_args.use_gpu);

  struct gkyl_moment_species fluid = {
    .name = "euler_mixture",
    .equation = euler_mixture,
    
    .init = evalEulerMixtureInit,
    .ctx = &ctx,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
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

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
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
    .name = "euler_mixture_fedkiw_shock",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx }, 
    .cells = { NX },

    .cfl_frac = ctx.cfl_frac,

    .num_species = 1,
    .species = { fluid },

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_moment_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_moment_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_moment_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_moment_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_moment_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr, false);

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = t_curr, .curr = frame_curr };

  write_data(&io_trig, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

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

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
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

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);

        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
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

freeresources:
  // Free resources after simulation completion.
  gkyl_wv_eqn_release(euler_mixture);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);

  gkyl_free(gas_gamma_s);
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
