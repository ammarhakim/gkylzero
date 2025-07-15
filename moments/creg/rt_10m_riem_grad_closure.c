// Generalized Brio-Wu Riemann problem, with gradient-based closure, for the 10-moment equations.
// Input parameters match the initial conditions found in entry JE4 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je4/je4-twofluid-shock.html), adapted from Section 7.1 of the article:
// A. Hakim, J. Loverich and U. Shumlak (2006), "A high resolution wave propagation scheme for ideal Two-Fluid plasma equations",
// Journal of Computational Physics, Volume 219 (1): 418-442.
// https://www.sciencedirect.com/science/article/pii/S0021999106001707

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

struct riem_grad_closure_ctx
{
  // Physical constants (using normalized code units).
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

  bool has_collision; // Whether to include collisions.
  double nu_base_ei; // Base electron-ion collision frequency.

  // Derived physical quantities (using normalized code units).
  double rhol_elc; // Left electron mass density.
  double rhor_elc; // Right electron mass density.

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  double Lx; // Domain size (x-direction).
  double k0; // Closure parameter.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  // Neural network parameters.
  bool use_nn_closure; // Use neural network-based closure?
  int poly_order; // Polynomial order of learned DG coefficients.
  const char* nn_closure_file; // File path of neural network to use.
};

struct riem_grad_closure_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
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

  bool has_collision = false; // Whether to include collisions.
  double nu_base_ei = 0.5; // Base electron-ion collision frequency.

  // Derived physical quantities (using normalized code units).
  double rhol_elc = rhol_ion * mass_elc / mass_ion; // Left electron mass density.
  double rhor_elc = rhor_ion * mass_elc / mass_ion; // Right electron mass density.

  // Simulation parameters.
  int Nx = 1024; // Cell count (x-direction).
  double Lx = 1.0; // Domain size (x-direction).
  double k0 = 500.0; // Closure parameter.
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 10.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  // Neural network parameters.
  bool use_nn_closure = false; // Use neural network-based closure?
  int poly_order = 1; // Polynomial order of learned DG coefficients.
  const char* nn_closure_file = "moments/data/neural_nets/pkpm_periodic_es_shock_p1_moms_nn_1"; // File path of neural network to use.
  
  struct riem_grad_closure_ctx ctx = {
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
    .has_collision = has_collision,
    .nu_base_ei = nu_base_ei,
    .rhol_elc = rhol_elc,
    .rhor_elc = rhor_elc,
    .Nx = Nx,
    .Lx = Lx,
    .k0 = k0,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .use_nn_closure = use_nn_closure,
    .poly_order = poly_order,
    .nn_closure_file = nn_closure_file,
  };

  return ctx;
}

void
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct riem_grad_closure_ctx *app = ctx;

  double rhol_elc = app->rhol_elc;
  double rhor_elc = app->rhor_elc;

  double pl = app->pl;
  double pr = app->pr;

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

  double mom_x = 0.0; // Electron momentum density (x-direction).
  double mom_y = 0.0; // Electron momentum density (y-direction).
  double mom_z = 0.0; // Electron momentum density (z-direction).

  double pr_xx = p; // Electron pressure tensor (xx-component).
  double pr_xy = 0.0; // Electron pressure tensor (xy-component).
  double pr_xz = 0.0; // Electron pressure tensor (xz-component).
  double pr_yy = p; // Electron pressure tensor (yy-component).
  double pr_yz = 0.0; // Electron pressure tensor (yz-component).
  double pr_zz = p; // Electron pressure tensor (zz-component).

  // Set electron mass density.
  fout[0] = rho;
  // Set electron momentum density.
  fout[1] = mom_x; fout[2] = mom_y; fout[3] = mom_z;
  // Set electron pressure tensor.
  fout[4] = pr_xx; fout[5] = pr_xy; fout[6] = pr_xz;
  fout[7] = pr_yy; fout[8] = pr_yz; fout[9] = pr_zz;
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct riem_grad_closure_ctx *app = ctx;

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

  double mom_x = 0.0; // Ion momentum density (x-direction).
  double mom_y = 0.0; // Ion momentum density (y-direction).
  double mom_z = 0.0; // Ion momentum density (z-direction).

  double pr_xx = p; // Ion pressure tensor (xx-component).
  double pr_xy = 0.0; // Ion pressure tensor (xy-component).
  double pr_xz = 0.0; // Ion pressure tensor (xz-component).
  double pr_yy = p; // Ion pressure tensor (yy-component).
  double pr_yz = 0.0; // Ion pressure tensor (yz-component).
  double pr_zz = p; // Ion pressure tensor (zz-component).

  // Set ion mass density.
  fout[0] = rho;
  // Set ion momentum density.
  fout[1] = mom_x; fout[2] = mom_y; fout[3] = mom_z;
  // Set ion pressure tensor.
  fout[4] = pr_xx; fout[5] = pr_xy; fout[6] = pr_xz;
  fout[7] = pr_yy; fout[8] = pr_yz; fout[9] = pr_zz;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct riem_grad_closure_ctx *app = ctx;

  double Bzl = app->Bzl;
  double Bzr = app->Bzr;

  double Ex = 0.0; // Total electric field (x-direction).
  double Ey = 0.0; // Total electric field (y-direction).
  double Ez = 0.0; // Total electric field (z-direction).
  
  double Bx = app->Bx; // Total magnetic field (x-direction).
  double By = 0.0; // Total magnetic field (y-direction).
  double Bz = 0.0;

  if (x < 0.5) {
    Bz = Bzl; // Total magnetic field (z-direction, left).
  }
  else {
    Bz = Bzr; // Total magnetic field (z-direction, right).
  }

  // Set electric field.
  fout[0] = Ex, fout[1] = Ey; fout[2] = Ez;
  // Set magnetic field.
  fout[3] = Bx, fout[4] = By; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
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

  struct riem_grad_closure_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);

  kann_t **ann = gkyl_malloc(sizeof(kann_t*) * 2);
  if (ctx.use_nn_closure) {
    const char *fmt_elc = "%s-%s.dat";
    int sz_elc = gkyl_calc_strlen(fmt_elc, ctx.nn_closure_file, "elc");
    char fileNm_elc[sz_elc + 1];
    snprintf(fileNm_elc, sizeof fileNm_elc, fmt_elc, ctx.nn_closure_file, "elc");
    FILE *file_elc = fopen(fileNm_elc, "r");
    if (file_elc != NULL) {
      ann[0] = kann_load(fileNm_elc);
      fclose(file_elc);
    }
    else {
      ann[0] = 0;
      ctx.use_nn_closure = false;
      fprintf(stderr, "Neural network for elc species not found! Disabling NN-based closure.\n");
    }

    const char *fmt_ion = "%s-%s.dat";
    int sz_ion = gkyl_calc_strlen(fmt_ion, ctx.nn_closure_file, "ion");
    char fileNm_ion[sz_ion + 1];
    snprintf(fileNm_ion, sizeof fileNm_ion, fmt_ion, ctx.nn_closure_file, "ion");
    FILE *file_ion = fopen(fileNm_ion, "r");
    if (file_ion != NULL) {
      ann[1] = kann_load(fileNm_ion);
      fclose(file_ion);
    }
    else {
      ann[1] = 0;
      ctx.use_nn_closure = false;
      fprintf(stderr, "Neural network for ion species not found! Disabling NN-based closure.\n");
    }
  }

  // Electron/ion equations.
  struct gkyl_wv_eqn *elc_ten_moment = gkyl_wv_ten_moment_new(ctx.k0, true, ctx.use_nn_closure, ctx.poly_order, ann[0], app_args.use_gpu);
  struct gkyl_wv_eqn *ion_ten_moment = gkyl_wv_ten_moment_new(ctx.k0, true, ctx.use_nn_closure, ctx.poly_order, ann[1], app_args.use_gpu);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .equation = elc_ten_moment,
    
    .init = evalElcInit,
    .ctx = &ctx,
  };

  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .equation = ion_ten_moment,
    
    .init = evalIonInit,
    .ctx = &ctx,
  };

  // Field.
  struct gkyl_moment_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    
    .init = evalFieldInit,
    .ctx = &ctx,
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
    .name = "10m_riem_grad_closure",

    .ndim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx }, 
    .cells = { NX },

    .cfl_frac = ctx.cfl_frac,

    .num_species = 2,
    .species = { elc, ion },

    .has_collision = ctx.has_collision,
    .nu_base = {
      { 0.0, ctx.nu_base_ei },
      { ctx.nu_base_ei, 0.0 }
    },

    .field = field,

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
  if (ctx.use_nn_closure) {
    for (int i = 0; i < app_inp.num_species; i++) {
      kann_delete(ann[i]);
    }
  }
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
  gkyl_wv_eqn_release(elc_ten_moment);
  gkyl_wv_eqn_release(ion_ten_moment);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);  
  gkyl_free(ann);
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
