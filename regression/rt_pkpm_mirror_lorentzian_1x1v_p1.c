#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm.h>
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

static inline double sq(double x) { return x*x; }
static inline double cu(double x) { return x*x*x; }

struct mirror_ctx {
  double epsilon0;
  double mu0;

  double charge_elc; // electron charge
  double mass_elc; // electron mass
  double charge_ion; // ion charge
  double mass_ion; // ion mass

  double vte; // electron thermal velocity
  double vti; // ion thermal velocity
  double n0; // initial number density

  double lambdaD;
  double wpe;
  double Lambda;

  double gamma; // FWHM of Lorentzian
  double loc; // location of Lorentzian
  double mag; // magnitude of Lorentzian

  double b_z0; // magnetic field at z = 0, for parameter checking
  double b_zL; // magnetic field at z = Lx, for parameter checking

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (velocity space: vx-direction).
  double Lx; // Domain size (configuration space: x-direction).
  double vx_max_elc; // Domain boundary (electron velocity space: vx-direction).
  double vx_max_ion; // Domain boundary (ion velocity space: vx-direction).
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

struct mirror_ctx
create_ctx(void)
{
  double epsilon0 = 8.8541878128e-12;
  double mass_elc = 9.109384e-31;
  double charge = 1.602177e-19;
  double n0 = 3e19;
  double vte = 1.0e7; //~1 kev electrons
  double vti = 7.0e5; //~10 kev deuterium
  double lambdaD = sqrt(epsilon0*vte*vte*mass_elc/(charge*charge*n0));

  double loc = 1.0;
  double gamma = 0.1;
  double mag = 1.0;
  double b_z0 = mag/(gamma*(1.0 + sq((loc)/gamma))) + mag/(gamma*(1.0 + sq((loc)/gamma)));
  double Lx = 2.0;
  double b_zL = mag/(gamma*(1.0 + sq((Lx + loc)/gamma))) + mag/(gamma*(1.0 + sq((Lx - loc)/gamma)));

  // Simulation parameters.
  int Nx = 128; // Cell count (configuration space: x-direction).
  int Nvx = 32; // Cell count (velocity space: vx-direction)..
  double vx_max_elc = 8.0 * vte; // Domain boundary (electron velocity space: vx-direction).
  double vx_max_ion = 8.0 * vti; // Domain boundary (ion velocity space: vx-direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 0.1; // CFL coefficient.

  double t_end = 1.0e-9; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = 1000; // Number of times to calculate field energy.
  int integrated_mom_calcs = 1000; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs = 1000; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct mirror_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = 1.256637062e-6,
    .n0 = n0,
    .charge_elc = -charge,
    .mass_elc = mass_elc,
    .charge_ion = charge,
    .mass_ion = 3.34449469e-27, // Deuterium mass

    .vte = vte, 
    .vti = vti, 
    // total length is 4 m (2 meters on each side with mirror throat at plus-minus 1 meter)
    // Debye length is 0.03 mm 
    // Plasma parameter ~ 1e6
    // Electron plasma frequency ~ 3e11/s
    .lambdaD = lambdaD,
    .wpe = vte/lambdaD,
    .Lambda = n0*lambdaD*lambdaD*lambdaD,

    .gamma = gamma,
    .loc = loc,
    .mag = mag,
    .b_z0 = b_z0,
    .b_zL = b_zL, 

    .Nx = Nx,
    .Nvx = Nvx,
    .Lx = Lx,
    .vx_max_elc = vx_max_elc,
    .vx_max_ion = vx_max_ion,
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

double calcJacobian(const double * GKYL_RESTRICT xn, void *ctx)
{
  struct mirror_ctx *app = ctx;
  double x = xn[0];
  double gamma = app->gamma;
  double loc = app->loc;
  double mag = app->mag;
  double magB = mag/(gamma*(1.0 + sq((x-loc)/gamma))) + mag/(gamma*(1.0 + sq((x+loc)/gamma)));
  // Jacobian for field-line-following coordinate system is 1/B
  return 1.0/magB;
}

void
evalJacobian(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mirror_ctx *app = ctx;
  double x = xn[0];
  double jacobgeo = calcJacobian(xn,ctx);
  fout[0] = jacobgeo;
}

void
evalMinusdBdzOverB(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mirror_ctx *app = ctx;
  double x = xn[0];
  double gamma = app->gamma;
  double loc = app->loc;
  double jacobgeo = calcJacobian(xn,ctx);
  double mag = -1.0*app->mag*jacobgeo; // overall prefactor is -1/B
  fout[0] = -mag*(2.0*(-loc + x))/(cu(gamma)*sq(1.0 + sq(-loc + x)/sq(gamma))) 
  - mag*(2.0*(loc + x))/(cu(gamma)*sq(1 + sq(loc + x)/sq(gamma)));
}

void
evalElcDistInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mirror_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vte, n0 = app->n0;
  double gamma = app->gamma;
  double loc = app->loc;
  // initial density profile is a tanh with similar FWHM to Lorentzian magnetic field
  double n = n0*(-tanh((x-loc)/gamma) + tanh((x+loc)/gamma));
  double fv = n/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  double jacobgeo = calcJacobian(xn,ctx);
  fout[0] = fv*jacobgeo;
  fout[1] = vt*vt*fv*jacobgeo;
}

void
evalElcFluidInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double mom_x = 0.0; // Electron total momentum density (x-direction).
  double mom_y = 0.0; // Electron total momentum density (y-direction).
  double mom_z = 0.0; // Electron total momentum density (z-direction).

  // Set electron total momentum density.
  fout[0] = mom_x; fout[1] = mom_y; fout[2] = mom_z;
}

void
evalIonDistInit(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct mirror_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vti, n0 = app->n0;
  double gamma = app->gamma;
  double loc = app->loc;
  // initial density profile is a tanh with similar FWHM to Lorentzian magnetic field
  double n = n0*(-tanh((x-loc)/gamma) + tanh((x+loc)/gamma));
  double fv = n/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  double jacobgeo = calcJacobian(xn,ctx);
  fout[0] = fv*jacobgeo;
  fout[1] = vt*vt*fv*jacobgeo;
}

void
evalIonFluidInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double mom_x = 0.0; // Ion total momentum density (x-direction).
  double mom_y = 0.0; // Ion total momentum density (y-direction).
  double mom_z = 0.0; // Ion total momentum density (z-direction).

  // Set ion total momentum density.
  fout[0] = mom_x; fout[1] = mom_y; fout[2] = mom_z;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double Ex = 0.0; // Total electric field (x-direction).
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
evalExternalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double Ex = 0.0; // External electric field (x-direction).
  double Ey = 0.0; // External electric field (y-direction).
  double Ez = 0.0; // External electric field (z-direction).

  // Dummy external B field for making magnetic field unit vector unity everywhere
  double Bx = 1.0; // External magnetic field (x-direction).
  double By = 0.0; // External magnetic field (y-direction).
  double Bz = 0.0; // External magnetic field (z-direction).

  // Set external electric field.
  fout[0] = Ex; fout[1] = Ey; fout[2] = Ez;
  // Set external magnetic field.
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;
}

void
evalElcNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct mirror_ctx *app = ctx;
  fout[0] = app->wpe/app->Lambda*log(app->Lambda);
}

void
evalIonNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct mirror_ctx *app = ctx;
  double Te_Ti = app->vte*app->vte*app->mass_elc/(app->vti*app->vti*app->mass_ion);
  double nu_ee = app->wpe/app->Lambda*log(app->Lambda);
  fout[0] = nu_ee/sqrt(app->mass_ion/app->mass_elc)*(Te_Ti*sqrt(Te_Ti));
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_pkpm_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_pkpm_app_write(app, t_curr, frame);
    gkyl_pkpm_app_write_field_energy(app);
    gkyl_pkpm_app_write_integrated_mom(app);
    gkyl_pkpm_app_write_integrated_L2_f(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_pkpm_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr)) {
    gkyl_pkpm_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_pkpm_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr)) {
    gkyl_pkpm_app_calc_integrated_mom(app, t_curr);
  }
}

void
calc_integrated_L2_f(struct gkyl_tm_trigger* l2t, gkyl_pkpm_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr)) {
    gkyl_pkpm_app_calc_integrated_L2_f(app, t_curr);
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

  struct mirror_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);

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
  struct gkyl_pkpm_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { -ctx.vx_max_elc },
    .upper = { ctx.vx_max_elc },
    .cells = { NVX },

    .init_dist = evalElcDistInit,
    .ctx_dist = &ctx,
    .init_fluid = evalElcFluidInit,
    .ctx_fluid = &ctx,

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalElcNu,
      .ctx = &ctx,
    },

    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
  };

  // Ions.
  struct gkyl_pkpm_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .lower = { -ctx.vx_max_ion },
    .upper = { ctx.vx_max_ion },
    .cells = { NVX },

    .init_dist = evalIonDistInit,
    .ctx_dist = &ctx,
    .init_fluid = evalIonFluidInit,
    .ctx_fluid = &ctx,

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalIonNu,
      .ctx = &ctx,
    },

    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB },
  };

  // Field.
  struct gkyl_pkpm_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .init = evalFieldInit,
    .ctx = &ctx,

    .ext_em = evalExternalFieldInit,
    .ext_em_ctx = &ctx,

    .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
  };

  // PKPM app.
  struct gkyl_pkpm app_inp = {
    .name = "pkpm_mirror_lorentzian_1x1v_p1",

    .cdim = 1, .vdim = 1,
    .lower = { -ctx.Lx },
    .upper = { ctx.Lx },
    .cells = { NX },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,
    .use_explicit_source = true, 

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .jacobgeo_flf_ctx = &ctx, 
    .jacobgeo_flf = evalJacobian, 
    .minus_dBdz_over_B_ctx = &ctx, 
    .minus_dBdz_over_B = evalMinusdBdzOverB, 

    .num_species = 2,
    .species = { elc, ion },

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&app_inp);

  // Print some parameters for reference
  gkyl_pkpm_app_cout(app, stdout, "Debye length = %lg\n", ctx.lambdaD);
  gkyl_pkpm_app_cout(app, stdout, "Plasma frequency = %lg\n", ctx.wpe);
  gkyl_pkpm_app_cout(app, stdout, "Plasma Parameter = %lg\n", ctx.Lambda);
  gkyl_pkpm_app_cout(app, stdout, "Electron-Electron collision frequency = %lg\n", ctx.wpe/ctx.Lambda*log(ctx.Lambda));
  // electron beta at z = 0
  gkyl_pkpm_app_cout(app, stdout, "Electron beta = %lg\n", 2.0*ctx.mu0*ctx.n0*ctx.vte*ctx.vte*ctx.mass_elc/(ctx.b_z0*ctx.b_z0));
  // electron gyroradius at z = 0
  gkyl_pkpm_app_cout(app, stdout, "Electron gyroradius = %lg\n", ctx.vte/(ctx.charge_ion*ctx.b_z0/ctx.mass_elc));

  // magnetic field at z = Lx
  gkyl_pkpm_app_cout(app, stdout, "Magnetic field at z = Lx, %lg\n", ctx.b_zL);
  // Debye length at z = Lx
  gkyl_pkpm_app_cout(app, stdout, "Debye length at z = Lx, %lg\n", ctx.lambdaD/sqrt(exp(-sq(ctx.Lx))));

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_pkpm_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_pkpm_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_pkpm_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_pkpm_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_pkpm_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr);

  // Create trigger for integrated L2 norm of the distribution function.
  int integrated_L2_f_calcs = ctx.integrated_L2_f_calcs;
  struct gkyl_tm_trigger l2f_trig = { .dt = t_end / integrated_L2_f_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_L2_f(&l2f_trig, app, t_curr);

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
    gkyl_pkpm_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_pkpm_update(app, dt);
    gkyl_pkpm_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_pkpm_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr);
    calc_integrated_mom(&im_trig, app, t_curr);
    calc_integrated_L2_f(&l2f_trig, app, t_curr);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_pkpm_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_pkpm_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_pkpm_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_pkpm_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_pkpm_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr);
  calc_integrated_mom(&im_trig, app, t_curr);
  calc_integrated_L2_f(&l2f_trig, app, t_curr);
  write_data(&io_trig, app, t_curr, false);
  gkyl_pkpm_app_stat_write(app);

  struct gkyl_pkpm_stat stat = gkyl_pkpm_app_stat(app);

  gkyl_pkpm_app_cout(app, stdout, "\n");
  gkyl_pkpm_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_pkpm_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_pkpm_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_pkpm_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_pkpm_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_pkpm_app_cout(app, stdout, "Fluid species RHD calc took %g secs\n", stat.fluid_species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species PKPM vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_pkpm_app_cout(app, stdout, "EM variables (bvar) calc took %g secs\n", stat.field_em_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm);
  gkyl_pkpm_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_pkpm_app_cout(app, stdout, "Number of write calls %ld\n", stat.nio);
  gkyl_pkpm_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_pkpm_app_release(app);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}