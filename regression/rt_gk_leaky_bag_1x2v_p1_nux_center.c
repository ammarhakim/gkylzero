#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_dynvec.h>
#include <gkyl_gyrokinetic.h>
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

struct boundary_ctx
{ 
  // Physical constants (using non-normalized physical units).
  double epsilon0; // Permittivity of free space.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.

  double Ti; // Ion temperature.
  double n0; // Reference number density (1 / m^3).

  // Derived physical quantities (using non-normalized physical units).
  double B0; // Reference magnetic field strength (Tesla).

  double vti; // Ion thermal velocity.

  // Simulation parameters.
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  double Lz; // Domain size (configuration space: z-direction).
  double vpar_max_ion; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion; // Domain boundary (ion velocity space: magnetic moment direction).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct boundary_ctx
create_ctx(void)
{
  int poly_order = 1; // Polynomial order.

  // Physical constants (using non-normalized physical units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mass_ion = 1.0;
  double charge_ion = 1.0; // Proton charge.

  double Ti = 1; // Ion temperature.
  double n0 = 1; //  Reference number density (1 / m^3).

  // Derived physical quantities (using non-normalized physical units).
  double B0 = 1; // Reference magnetic field strength (Tesla).
  
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.

  // Simulation parameters.
  int Nz = 16; // Cell count (configuration space: z-direction).
  int Nvpar = 32; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 32; // Cell count (velocity space: magnetic moment direction).
  double Lz = 1.0; // Domain size (configuration space: z-direction).
  double vpar_max_ion = 4.0 * vti; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion = (3.0 / 2.0) * 0.5 * mass_ion * pow(4.0 * vti, 2.0) / (2.0 * B0); // Domain boundary (ion velocity space: magnetic moment direction).
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 3; // Final simulation time.
  int num_frames = 5; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct boundary_ctx ctx = {
    .epsilon0 = epsilon0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .Ti = Ti,
    .n0 = n0,
    .B0 = B0,
    .vti = vti,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .Lz = Lz,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .poly_order = poly_order,
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


void compareToAnalytics(const struct gkyl_gk *app_inp, void* ctx )
{
  struct boundary_ctx *app = ctx;
  char filename[256];
  sprintf(filename, "%s-ion_integrated_moms.gkyl", app_inp->name);

  struct gkyl_dynvec_etype_ncomp enc = { };
  enc = gkyl_dynvec_read_ncomp(filename);

  gkyl_dynvec integrated_moms = gkyl_dynvec_new(enc.type, enc.ncomp);
  gkyl_dynvec_read(integrated_moms, filename);

  size_t units = gkyl_dynvec_size(integrated_moms);

  for (int i = 0; i < units; i++) {
    double data[enc.ncomp];
    gkyl_dynvec_get(integrated_moms, i, data);
    double time = gkyl_dynvec_get_tm(integrated_moms, i);
    double N = data[0];

    double vt = sqrt(app->Ti / app->mass_ion);
    double L = app->Lz;
    double erf_arg = L / (sqrt(2) * time * vt);

    double N_analytic = 2 * app->n0 / (sqrt(2 * M_PI));
    N_analytic *= sqrt(M_PI / 2) * L * erf(erf_arg) + time * vt * (exp(-erf_arg*erf_arg) - 1);

    double diff = fabs(N - N_analytic);
    int check = gkyl_compare_double(N, N_analytic, 1e-2);
    if (check != 1) {
      printf("Error: N and N_analytic do not match within tolerance.\n");
      printf("N - N_analytic: %g\n", N - N_analytic);
    }
  }

  sprintf(filename, "%s-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl", app_inp->name);
  enc = gkyl_dynvec_read_ncomp(filename);
  gkyl_dynvec bflux_moms = gkyl_dynvec_new(enc.type, enc.ncomp);
  gkyl_dynvec_read(bflux_moms, filename);
  units = gkyl_dynvec_size(bflux_moms);

  for (int i = 1; i < units; i++) {
    double data[enc.ncomp];
    gkyl_dynvec_get(bflux_moms, i, data);
    double time = gkyl_dynvec_get_tm(bflux_moms, i);
    double dNdt = data[0];

    double vt = sqrt(app->Ti / app->mass_ion);
    double L = app->Lz;
    double erf_arg = L / (sqrt(2) * time * vt);

    double dNdt_analytic = app->n0 / (sqrt(2 * M_PI));
    dNdt_analytic *= 1 - exp(-erf_arg*erf_arg);

    // printf("dNdt: %g, dNdt_analytic: %g, difference: %g\n", dNdt, dNdt_analytic, fabs(dNdt - dNdt_analytic));

    double diff = fabs(dNdt - dNdt_analytic);
    int check = gkyl_compare_double(dNdt, dNdt_analytic, 1e-1);
    if (check != 1) {
      printf("Error: dNdt and dNdt_analytic do not match within tolerance.\n");
      printf("N - N_analytic: %g\n", dNdt - dNdt_analytic);
    }
  }

  sprintf(filename, "%s-ion_bflux_xupper_integrated_HamiltonianMoments.gkyl", app_inp->name);
  enc = gkyl_dynvec_read_ncomp(filename);
  gkyl_dynvec bflux_moms_upper = gkyl_dynvec_new(enc.type, enc.ncomp);
  gkyl_dynvec_read(bflux_moms_upper, filename);
  units = gkyl_dynvec_size(bflux_moms_upper);

  for (int i = 1; i < units; i++) {
    double data[enc.ncomp];
    gkyl_dynvec_get(bflux_moms_upper, i, data);
    double time = gkyl_dynvec_get_tm(bflux_moms_upper, i);
    double dNdt = data[0];

    double vt = sqrt(app->Ti / app->mass_ion);
    double L = app->Lz;
    double erf_arg = L / (sqrt(2) * time * vt);

    double dNdt_analytic = app->n0 / (sqrt(2 * M_PI));
    dNdt_analytic *= 1 - exp(-erf_arg*erf_arg);

    double diff = fabs(dNdt - dNdt_analytic);
    int check = gkyl_compare_double(dNdt, dNdt_analytic, 1e-1);
    if (check != 1) {
      printf("Error: dNdt and dNdt_analytic do not match within tolerance.\n");
      printf("N - N_analytic: %g\n", dNdt - dNdt_analytic);
    }
  }

  printf("If there are no errors above, then the test passed.\n");

  gkyl_dynvec_release(integrated_moms);
  gkyl_dynvec_release(bflux_moms);
  gkyl_dynvec_release(bflux_moms_upper);
}

void
evalIonDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct boundary_ctx *app = ctx;
  fout[0] = app->n0;
}

void
evalIonTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct boundary_ctx *app = ctx;
  fout[0] = app->Ti;
}

void
evalIonUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = 0.0;
}

static inline void
nonuniform_position_map_z(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  struct boundary_ctx *app = ctx;
  double z = zc[0];
  xp[0] = z - 0.1 * sin(z * 2 * M_PI/(app->Lz));
}

static inline void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  // Set physical coordinates (X, Y, Z) from computational coordinates (x, y, z).
  xp[0] = zc[0]; xp[1] = zc[1]; xp[2] = zc[2];
}

void
bmag_func(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void* ctx)
{
  struct boundary_ctx *app = ctx;
  fout[0] = app->B0;
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_gyrokinetic_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_gyrokinetic_app* app, double t_curr, double dt, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);

    if ( !(dt < 0.0) )
      gkyl_gyrokinetic_app_save_dt(app, t_curr, dt);
  }
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_gyrokinetic_app_write(app, t_curr, frame);
    gkyl_gyrokinetic_app_write_field_energy(app);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
    gkyl_gyrokinetic_app_write_dt(app);
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

  struct boundary_ctx ctx = create_ctx(); // Context for initialization functions.

  int NZ = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nz);
  int NVPAR = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvpar);
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nmu);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int ccells[] = { NZ };
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

  // Ions.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .lower = { -ctx.vpar_max_ion, 0.0 },
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion },
    .cells = { NVPAR, NMU },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .density = evalIonDensityInit,
      .ctx_density = &ctx,
      .temp = evalIonTempInit,
      .ctx_temp = &ctx,
      .upar = evalIonUparInit,
      .ctx_upar = &ctx,
    },

    .bcx = {
      .lower = { .type = GKYL_SPECIES_ABSORB, },
      .upper = { .type = GKYL_SPECIES_ABSORB, },
    },
    
    .num_diag_moments = 6,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_MAXWELLIAN },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    }
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .zero_init_field = true, // Don't compute the field at t = 0.
    .is_static = true, // Don't evolve the field in time.
  };

  // Gyrokinetic app.
  struct gkyl_gk app_inp = {
    .name = "gk_leaky_bag_1x2v_p1_nux_center",

    .cdim = 1, .vdim = 2,
    .lower = { -ctx.Lz/2.0 },
    .upper = {  ctx.Lz/2.0 },
    .cells = { NZ },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .world = { 0.0, 0.0 },
      .mapc2p = mapc2p,
      .c2p_ctx = &ctx,
      .bmag_func = bmag_func,
      .bmag_ctx = &ctx,
      .position_map_info = {
        .maps[2] = nonuniform_position_map_z,
        .ctxs[2] = &ctx,
      },
    },


    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 1,
    .species = { ion },

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr, -1.0, false);

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
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, status.dt_actual, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, status.dt_actual, true);
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
  calc_integrated_mom(&im_trig, app, t_curr, -1.0, false);
  write_data(&io_trig, app, t_curr, false);

  gkyl_gyrokinetic_app_stat_write(app);
  
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  compareToAnalytics(&app_inp, &ctx);

freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  gkyl_gyrokinetic_comms_release(comm);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
