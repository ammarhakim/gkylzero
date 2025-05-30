 #include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct escreen_ctx {
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double n; // reference number density
  double injRate;
  double T; // electron temperature
  double vdrift_begin; // drift velocity
  double vdrift_end; // drift velocity 
  double E0; // electric field amplitude
  double noise_amp; // Noise level for perturbation.
  int mode_init; // Initial wave mode to perturb with noise.
  int mode_final; // Final wave mode to perturb with noise.
  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  int Np; // Cell countr (p-direction).
  double Lx; // Domain size (x-direction).
  double pmax; // Momentum space extents for electrons and positrons
  double nonuniform_p_pow; // How nonuniform is momentum (four-velocity)-space? 
                           // momentum (four-velocity)-space grid points are spaced as p^n power.
  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

void
evalDensityPerturbInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct escreen_ctx *app = ctx;
  double x = xn[0];

  // Perturbation to initial density to satisfy div(E) = rho_c
  double E_x = app->E0;
  double width = app->Lx/4;
  fout[0] = E_x * (2 * x/(width*width)) * exp(-(x/width)*(x/width));
}

void
evalDensityPerturbPos(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; 
}

void
evalVDriftInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct escreen_ctx *app = ctx;

  // Set drift (four-) velocity.
  fout[0] = 0.0;
}

void
evalTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct escreen_ctx *app = ctx;

  double T = app->T;

  // Set temperature.
  fout[0] = T;
}

void
evalDensitySource(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct escreen_ctx *app = ctx;
  double x = xn[0], p = xn[1];
  double t_end = app->t_end;
  double T = app->T, vdriftB = app->vdrift_begin, vdriftE=app->vdrift_end;
  double vdrift = vdriftB+(vdriftE-vdriftB)*t/t_end;
  double inj = app->injRate;
  double center = 0.075*app->Lx+t;
  double width = app->Lx/30;
  double xLow = center-3*width;
  double xHigh = center+3*width;

  // Lorentz factor for drift velocity
  double gamma = 1.0/sqrt(1 - vdrift*vdrift);

  double n = app->n;
  double fv = 0.0;
  if (x>xLow && x<xHigh){
    fv = (1/gamma)*n*exp(-((x-center)/width)*((x-center)/width));
  }

  // Set density.
  fout[0] = inj*fv;
}

void
evalVDriftSource(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct escreen_ctx *app = ctx;
  double x = xn[0], p = xn[1];
  double t_end = app->t_end;
  double T = app->T, vdriftB = app->vdrift_begin, vdriftE=app->vdrift_end;
  double vdrift = vdriftB+(vdriftE-vdriftB)*t/t_end;
  double inj = app->injRate;
  double center = 0.075*app->Lx+t;
  double width = app->Lx/30;
  double xLow = center-3*width;
  double xHigh = center+3*width;

  // Lorentz factor for drift velocity
  double gamma = 1.0/sqrt(1 - vdrift*vdrift);

  double fv = 0.0;
  if (x>xLow && x<xHigh){
    fv = gamma*vdrift;
  }

  // Set drift (four-) velocity.
  fout[0] = fv;
}

void
evalTempSource(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct escreen_ctx *app = ctx;
  double x = xn[0], p = xn[1];
  double t_end = app->t_end;
  double T = app->T, vdriftB = app->vdrift_begin, vdriftE=app->vdrift_end;
  double vdrift = vdriftB+(vdriftE-vdriftB)*t/t_end;
  double inj = app->injRate;
  double center = 0.075*app->Lx+t;
  double width = app->Lx/30;
  double xLow = center-3*width;
  double xHigh = center+3*width;

  // Lorentz factor for drift velocity
  double gamma = 1.0/sqrt(1 - vdrift*vdrift);

  double fv = 0.0;
  if (x>xLow && x<xHigh){
    fv = T;
  }

  // Set temperature.
  fout[0] = fv;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct escreen_ctx *app = ctx;
  double x = xn[0];
  double pi = app->pi;

  double width = app->Lx/4;

  double E_x = app->E0 * exp(-(x/width)*(x/width));
  double Lx = app->Lx;

  fout[0] = E_x; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

static void
mapc2p_vel(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct escreen_ctx *app = ctx;
  double pmax = app->pmax;
  double nonuniform_p_pow = app->nonuniform_p_pow;

  if (vc[0] < 0.0) {
    vp[0] = -pmax*pow(vc[0], nonuniform_p_pow);
  }
  else {
    vp[0] =  pmax*pow(vc[0], nonuniform_p_pow);  
  }
}

struct escreen_ctx
create_default_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 1.0; // Positron mass.
  double charge_ion = 1.0; // Positron charge.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double injRate = 0.035;//0.0175;  // injection rate

  double n = 1.0; // Reference density
  double vdrift_begin = 0.995;
  double vdrift_end = 0.995; //0.9428;
  //double vdrift = 0.995; //0.9428; //0.995;//0.99; // Initial drift velocity
  double T = 0.25; // Reference temperature (in units of mc^2)

  double E0 = 2.4; // Electric field amplitude
  double noise_amp = 0.0; //5.0e-5 * E0; // Noise level for perturbation.
  int mode_init = 1; // Initial wave mode to perturb with noise.
  int mode_final = 4; // Final wave mode to perturb with noise.

  // Simulation parameters.
  int Nx = 10240;// 40960; // Cell count (x-direction).
  int Np = 5120; // 512; // Cell count (p-direction).
  double Lx = 1000.0; // Domain size (x-direction).
  double pmax = 150.0; // Momentum space extents for electrons and positrons  
  double nonuniform_p_pow = 2.0; // How nonuniform is momentum space: grid points are spaced as p^n power
  double t_end = 1000.0; // Final simulation time.
  int num_frames = 100; // Number of output frames.
  int field_energy_calcs = num_frames*100; // Number of times to calculate field energy.
  int integrated_mom_calcs = num_frames*100; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs = num_frames*100; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  return (struct escreen_ctx) {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .n = n,
    .injRate = injRate,
    .vdrift_begin = vdrift_begin,
    .vdrift_end = vdrift_end,  
    .T = T,
    .E0 = E0, 
    .noise_amp = noise_amp,
    .mode_init = mode_init,
    .mode_final = mode_final,
    .Nx = Nx,
    .Np = Np, 
    .Lx = Lx,
    .pmax = pmax, 
    .nonuniform_p_pow = nonuniform_p_pow, 
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .integrated_L2_f_calcs = integrated_L2_f_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
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

  struct escreen_ctx ctx;
  ctx = create_default_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Np);

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

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .model_id = GKYL_MODEL_SR,
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { -1.0 },
    .upper = { 1.0 }, 
    .cells = { VX },
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ABSORB, },
      .upper = { .type = GKYL_SPECIES_ABSORB, },
    },

    .num_init = 1, 
    .projection[0] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
      .density = evalDensityPerturbInit,
      .ctx_density = &ctx,
      .temp = evalTempInit,
      .ctx_temp = &ctx,
      .V_drift = evalVDriftInit,
      .ctx_V_drift = &ctx,
      .correct_all_moms = true,
      .use_last_converged = true, 
    },

    .source = {
      .source_id = GKYL_PROJ_ADAPT_DENSITY_SOURCE,
      .write_source = true, 

      .num_cross_source = 1, 
      .source_with = { "pos" }, 
      .source_with_v_thresh = { 10.0 }, // threshold velocity for partial moment
      .source_with_upper_half = { true }, // is the integral over the upper-half plane?
      .source_with_proj = { 0 }, 

      .num_sources = 1, 
      .projection[0] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalDensitySource,
        .ctx_density = &ctx,
        .temp = evalTempSource,
        .ctx_temp = &ctx,
        .V_drift = evalVDriftSource,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
        .max_iter = 10, 
        .iter_eps = 1.0e-6, 
        .use_last_converged = true, 
      },
    },

    .mapc2p_vel[0] = {
      .mapc2p_vel_func = mapc2p_vel, 
      .mapc2p_vel_ctx = &ctx, 
    },

    .num_diag_moments = 2,
    .diag_moments = { "M0", "M1i" },
  };

  // positrons
  struct gkyl_vlasov_species pos = {
    .name = "pos",
    .model_id = GKYL_MODEL_SR,
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .lower = { -1.0 },
    .upper = { 1.0 }, 
    .cells = { VX },
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ABSORB, },
      .upper = { .type = GKYL_SPECIES_ABSORB, },
    },

    .num_init = 1, 
    .projection[0] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
      .density = evalDensityPerturbPos,
      .ctx_density = &ctx,
      .temp = evalTempInit,
      .ctx_temp = &ctx,
      .V_drift = evalVDriftInit,
      .ctx_V_drift = &ctx,
      .correct_all_moms = true,
      .use_last_converged = true, 
    },

    .source = {
      .source_id = GKYL_PROJ_ADAPT_DENSITY_SOURCE,
      .write_source = true, 

      .num_cross_source = 1, 
      .source_with = { "elc" }, 
      .source_with_v_thresh = { 10.0 }, // threshold velocity for partial moment
      .source_with_upper_half = { true }, // is the integral over the upper-half plane?
      .source_with_proj = { 0 }, 

      .num_sources = 1, 
      .projection[0] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = evalDensitySource,
        .ctx_density = &ctx,
        .temp = evalTempSource,
        .ctx_temp = &ctx,
        .V_drift = evalVDriftSource,
        .ctx_V_drift = &ctx,
        .correct_all_moms = true,
        .max_iter = 10, 
        .iter_eps = 1.0e-6, 
        .use_last_converged = true, 
      },
    },

    .mapc2p_vel[0] = {
      .mapc2p_vel_func = mapc2p_vel, 
      .mapc2p_vel_ctx = &ctx, 
    },

    .num_diag_moments = 2,
    .diag_moments = { "M0", "M1i" },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,
    .bcx = { GKYL_FIELD_COPY, GKYL_FIELD_COPY },
  };

  // Vlasov-Maxwell app.
  struct gkyl_vm app_inp = {
    .name = "escreen_sr",
    
    .cdim = 1, .vdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    //.num_periodic_dir = 1,
    //.periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, pos },
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
