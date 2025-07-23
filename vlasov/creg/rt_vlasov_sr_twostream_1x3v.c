#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <gkyl_util.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct twostream_sr_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double n_elc1; // First electron number density.
  double n_elc2; // Second electron number density.
  double ux_elc1; // First electron velocity (x-direction).
  double ux_elc2; // Second electron velocity (x-direction).
  double uy_elc1; // First electron velocity (y-direction).
  double uy_elc2; // Second electron velocity (y-direction).
  double uz_elc1; // First electron velocity (z-direction).
  double uz_elc2; // Second electron velocity (z-direction).
  double T_elc1; // First electron temperature (units of mc^2).
  double T_elc2; // Second electron temperature (units of mc^2).

  double alpha; // Applied perturbation amplitude.
  double kx; // Perturbed wave number (x-direction).

  // Derived physical quantities (using normalized code units).
  double gamma_elc1; // First electron gamma factor.
  double gamma_elc2; // Second electron gamma factor.

  double ux_elc1_sr; // First electron relativistic velocity (x-direction).
  double ux_elc2_sr; // Second electron relativistic velocity (x-direction).
  double uy_elc1_sr; // First electron relativistic velocity (y-direction).
  double uy_elc2_sr; // Second electron relativistic velocity (y-direction).
  double uz_elc1_sr; // First electron relativistic velocity (z-direction).
  double uz_elc2_sr; // Second electron relativistic velocity (z-direction).

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (velocity space: vx-direction).
  int Nvy; // Cell count (velocity space: vy-direction).
  int Nvz; // Cell count (velocity space: vz-direction).
  double Lx; // Domain size (configuration space: x-direction).
  double vx_max; // Domain boundary (velocity space: vx-direction).
  double vy_max; // Domain boundary (velocity space: vy-direction).
  double vz_max; // Domain boundary (velocity space: vz-direction).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct twostream_sr_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.

  double n_elc1 = 0.5; // First electron number density.
  double n_elc2 = 0.5; // Second electron number density.
  double ux_elc1 = 0.9; // First electron velocity (x-direction).
  double ux_elc2 = -0.9; // Second electron velocity (x-direction).
  double uy_elc1 = 0.0; // First electron velocity (y-direction).
  double uy_elc2 = 0.0; // Second electron velocity (y-direction).
  double uz_elc1 = 0.0; // First electron velocity (z-direction).
  double uz_elc2 = 0.0; // Second electron velocity (z-direction).
  double T_elc1 = 0.04; // First electron temperature (units of mc^2).
  double T_elc2 = 0.04; // Second electron temperature (units of mc^2).

  double alpha = 1.0e-3; // Applied perturbation amplitude.
  double kx = 0.5; // Perturbed wave number (x-direction).

  // Derived physical quantities (using normalized code units).
  double gamma_elc1 = 1.0 / sqrt(1.0 - (ux_elc1 * ux_elc1) - (uy_elc1 * uy_elc1) - (uz_elc1 * uz_elc1)); // First electron gamma factor.
  double gamma_elc2 = 1.0 / sqrt(1.0 - (ux_elc2 * ux_elc2) - (uy_elc2 * uy_elc2) - (uz_elc2 * uz_elc2)); // Second electron gamma factor.

  double ux_elc1_sr = gamma_elc1 * ux_elc1; // First electron relativistic velocity (x-direction).
  double ux_elc2_sr = gamma_elc2 * ux_elc2; // Second electron relativistic velocity (x-direction).
  double uy_elc1_sr = gamma_elc1 * uy_elc1; // First electron relativistic velocity (y-direction).
  double uy_elc2_sr = gamma_elc2 * uy_elc2; // Second electron relativistic velocity (y-direction).
  double uz_elc1_sr = gamma_elc1 * uz_elc1; // First electron relativistic velocity (z-direction).
  double uz_elc2_sr = gamma_elc2 * uz_elc2; // Second electron relativistic velocity (z-direction).

  // Simulation parameters.
  int Nx = 32; // Cell count (configuration space: x-direction).
  int Nvx = 16; // Cell count (velocity space: vx-direction).
  int Nvy = 16; // Cell count (velocity space: vy-direction).
  int Nvz = 16; // Cell count (velocity space: vz-direction).
  double Lx = 2.0 * pi / kx; // Domain size (configuration space: x-direction).
  double vx_max = 8.0; // Domain boundary (velocity space: vx-direction).
  double vy_max = 8.0; // Domain boundary (velocity space: vy-direction).
  double vz_max = 8.0; // Domain boundary (velocity space: vz-direction).
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 1.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs = INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct twostream_sr_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .n_elc1 = n_elc1,
    .n_elc2 = n_elc2,
    .ux_elc1 = ux_elc1,
    .ux_elc2 = ux_elc2,
    .uy_elc1 = uy_elc1,
    .uy_elc2 = uy_elc2,
    .uz_elc1 = uz_elc1,
    .uz_elc2 = uz_elc2,
    .T_elc1 = T_elc1,
    .T_elc2 = T_elc2,
    .alpha = alpha,
    .kx = kx,
    .gamma_elc1 = gamma_elc1,
    .gamma_elc2 = gamma_elc2,
    .ux_elc1_sr = ux_elc1_sr,
    .ux_elc2_sr = ux_elc2_sr,
    .uy_elc1_sr = uy_elc1_sr,
    .uy_elc2_sr = uy_elc2_sr,
    .uz_elc1_sr = uz_elc1_sr,
    .uz_elc2_sr = uz_elc2_sr,
    .Nx = Nx,
    .Nvx = Nvx,
    .Nvy = Nvy,
    .Nvz = Nvz,
    .Lx = Lx,
    .vx_max = vx_max,
    .vy_max = vy_max,
    .vz_max = vz_max,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .integrated_L2_f_calcs = integrated_L2_f_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalDensityLInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct twostream_sr_ctx *app = ctx;
  double x = xn[0];

  double alpha = app->alpha;
  double kx = app->kx;
  double n_elc1 = app->n_elc1;

  // Set left-going total number density.
  fout[0] = (1.0 + alpha * cos(kx * x)) * n_elc1;
}

void
evalDensityRInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct twostream_sr_ctx *app = ctx;
  double x = xn[0];

  double alpha = app->alpha;
  double kx = app->kx;
  double n_elc2 = app->n_elc2;

  // Set right-going total number density.
  fout[0] = (1.0 + alpha * cos(kx * x)) * n_elc2;
}

void
evalTempLInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct twostream_sr_ctx *app = ctx;

  double T_elc1 = app->T_elc1;

  // Set left-going temperature.
  fout[0] = T_elc1;
}

void
evalTempRInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct twostream_sr_ctx *app = ctx;

  double T_elc2 = app->T_elc2;

  // Set right-going temperature.
  fout[0] = T_elc2;
}

void
evalVDriftLInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct twostream_sr_ctx *app = ctx;

  double ux_elc1_sr = app->ux_elc1_sr;
  double uy_elc1_sr = app->uy_elc1_sr;
  double uz_elc1_sr = app->uz_elc1_sr;

  // Set left-going relativistic drift velocity.
  fout[0] = ux_elc1_sr; fout[1] = uy_elc1_sr; fout[2] = uz_elc1_sr;
}

void
evalVDriftRInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct twostream_sr_ctx *app = ctx;

  double ux_elc2_sr = app->ux_elc2_sr;
  double uy_elc2_sr = app->uy_elc2_sr;
  double uz_elc2_sr = app->uz_elc2_sr;

  // Set right-going relativistic drift velocity.
  fout[0] = ux_elc2_sr; fout[1] = uy_elc2_sr; fout[2] = uz_elc2_sr;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct twostream_sr_ctx *app = ctx;
  double x = xn[0];

  double alpha = app->alpha;
  double kx = app->kx;

  double Ex = -alpha * sin(kx * x) / kx; // Total electric field (x-direction).
  double Ey = 0.0; // Total electric field (y-direction).
  double Ez = 0.0; // Total electric field (z-direction).

  double Bx = 0.0; // Total magnetic field (x-direction).
  double By = 0.0; // Total magnetic field (y-direction).
  double Bz = 0.0; // Total magnetic field (z-direction).
  
  // Set electric field.
  fout[0] = Ex; fout[1] = Ey, fout[2] = Ez;
  // Set magnetic field.
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_vlasov_app_write(app, t_curr, frame);
    gkyl_vlasov_app_write_field_energy(app);
    gkyl_vlasov_app_write_integrated_mom(app);
    gkyl_vlasov_app_write_integrated_L2_f(app);

    gkyl_vlasov_app_calc_mom(app);
    gkyl_vlasov_app_write_mom(app, t_curr, frame);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  }
}

void
calc_integrated_L2_f(struct gkyl_tm_trigger* l2t, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_L2_f(app, t_curr);
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

  struct twostream_sr_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);
  int NVY = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nvy);
  int NVZ = APP_ARGS_CHOOSE(app_args.vcells[2], ctx.Nvz);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif  

  int ccells[] = { NX };
  int cdim = sizeof(ccells) / sizeof(ccells[0]);

  int cuts[cdim];
#ifdef GKYL_HAVE_MPI  
  for (int d = 0; d < cdim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < cdim; d++) {
    cuts[d] = 1;
  }
#endif  
    
  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_gpu && app_args.use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
#else
    printf(" Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (app_args.use_mpi) {
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
  for (int d = 0; d < cdim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Electrons.
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_SR,
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { -ctx.vx_max, -ctx.vy_max, -ctx.vz_max, },
    .upper = { ctx.vx_max, ctx.vy_max, ctx.vz_max }, 
    .cells = { NVX, NVY, NVZ },

    .num_init = 2, 
    // Two counter-streaming Maxwellians.
    .projection[0] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
      .density = evalDensityLInit,
      .ctx_density = &ctx,
      .temp = evalTempLInit,
      .ctx_temp = &ctx,
      .V_drift = evalVDriftLInit,
      .ctx_V_drift = &ctx,
      .correct_all_moms = true,
      .use_last_converged = true, 
    },
    .projection[1] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
      .density = evalDensityRInit,
      .ctx_density = &ctx,
      .temp = evalTempRInit,
      .ctx_temp = &ctx,
      .V_drift = evalVDriftRInit,
      .ctx_V_drift = &ctx,
      .correct_all_moms = true,
      .use_last_converged = true, 
    },

    .num_diag_moments = 2,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1 },
  };

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,

    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .init = evalFieldInit,
    .ctx = &ctx,
  };

  // Vlasov-Maxwell app.
  struct gkyl_vm app_inp = {
    .name = "vlasov_sr_twostream_1x3v",
    
    .cdim = 1, .vdim = 3,
    .lower = { -0.5 * ctx.Lx },
    .upper = { 0.5 * ctx.Lx },
    .cells = { NX },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 1,
    .species = { elc },

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };
  
  // Create app object.
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_vlasov_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_vlasov_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_vlasov_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_vlasov_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_vlasov_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr, false);

  // Create trigger for integrated L2 norm of the distribution function.
  int integrated_L2_f_calcs = ctx.integrated_L2_f_calcs;
  struct gkyl_tm_trigger l2f_trig = { .dt = t_end / integrated_L2_f_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);

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
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
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

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        calc_integrated_L2_f(&l2f_trig, app, t_curr, true);
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
  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
  gkyl_vlasov_app_stat_write(app);

  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  gkyl_vlasov_app_cout(app, stdout, "\n");
  gkyl_vlasov_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_vlasov_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_vlasov_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_vlasov_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

freeresources:
  // Free resources after simulation completion.
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
