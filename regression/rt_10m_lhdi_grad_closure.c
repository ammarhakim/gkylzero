// Lower-Hybrid Drift Instability (LHDI) test, with gradient-closure, for the 10-moment equations.
// Input parameters match the initial conditions in Section 5.2, from the article:
// J. Ng, A. Hakim, J. Juno and A. Bhattacharjee (2019), "Drift Instabilities in Thin Current Sheets Using a Two-Fluid Model With Pressure Tensor Effects",
// Journal of Geophysical Research: Space Physics, Volume 124 (5): 3331-3346.
// https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JA026313

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

struct lhdi_grad_closure_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_ion; // Ion mass.
  double charge_ion; // Ion charge.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double Te_over_Ti; // Electron temperature / ion temperature.

  double n0; // Reference number density.
  double nb_over_n0; // Background number density / reference number density.
  double vte; // Electron thermal velocity.

  double beta; // Electron plasma beta.

  double noise_amp; // Noise level for perturbation.
  double mode; // Wave mode to perturb with noise.

  // Derived physical quantities (using normalized code units).
  double Te; // Electron temperature.
  double Ti; // Ion temperature.
  double vti; // Ion thermal velocity.

  double vAe; // Electron Alfven velocity.
  double B0; // Reference magnetic field strength (derived from normalization of mass_elc and n0).
  double vAi; // Ion Alfven velocity.

  double omega_ci; // Ion cyclotron frequency.
  double omega_ce; // Electron cyclotron frequency.

  double larmor_ion; // Ion Larmor radius.
  double larmor_elc; // Electron Larmor radius.

  double l; // Current sheet width.
  
  double nb; // Background number density.
  double Te_frac; // Fraction of total temperature from electrons.
  double Ti_frac; // Fraction of total temperature from ions.
  double ix; // Current (x-direction).
  double iy; // Current (y-direction).

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  int Ny; // Cell count (y-direction).
  double Lx; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
  double k0_elc; // Closure parameter for electrons.
  double k0_ion; // Closure parameter for ions.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct lhdi_grad_closure_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 36.0; // Ion mass.
  double charge_ion = 1.0; // Ion charge.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double Te_over_Ti = 0.1; // Electron temperature / ion temperature.

  double n0 = 1.0; // Reference number density.
  double nb_over_n0 = 0.001; // Background number density / reference number density.
  double vte = 0.06; // Electron thermal velocity.

  double beta = 1.0 / 11.0; // Electron plasma beta.

  double noise_amp = 0.0001; // Noise level for perturbation.
  double mode = 8.0; // Wave mode to perturb with noise.

  // Derived physical quantities (using normalized code units).
  double Te = vte * vte * mass_elc / 2.0; // Electron temperature.
  double Ti = Te / Te_over_Ti; // Ion temperature.
  double vti = sqrt(2.0 * Ti / mass_ion); // Ion thermal velocity.

  double vAe = vte / sqrt(beta); // Electron Alfven velocity.
  double B0 = vAe; // Reference magnetic field strength (derived from normalization of mass_elc and n0).
  double vAi = vAe / sqrt(mass_ion); // Ion Alfven velocity.

  double omega_ci = charge_ion * B0 / mass_ion; // Ion cyclotron frequency.
  double omega_ce = charge_ion * B0 / mass_elc; // Electron cyclotron frequency.

  double larmor_ion = vti / omega_ci; // Ion Larmor radius.
  double larmor_elc = vte / omega_ce; // Electron Larmor radius.

  double l = larmor_ion; // Current sheet width.
  
  double nb = n0 * nb_over_n0; // Background number density.
  double Te_frac = Te / (Te + Ti); // Fraction of total temperature from electrons.
  double Ti_frac = 1.0 - Te_frac; // Fraction of total temperature from ions.
  double ix = mode; // Current (x-direction).
  double iy = 1.0; // Current (y-direction).

  // Simulation parameters.
  int Nx = 64; // Cell count (x-direction).
  int Ny = 128; // Cell count (y-direction).
  double Lx = 6.4 * l; // Domain size (x-direction).
  double Ly = 12.8 * l; // Domain size (y-direction).
  double k0_elc = 1.0; // Closure parameter for electrons.
  double k0_ion = 1.0 / 6.0; // Closure parameter for ions.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 1100.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct lhdi_grad_closure_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .Te_over_Ti = Te_over_Ti,
    .n0 = n0,
    .nb_over_n0 = nb_over_n0,
    .vte = vte,
    .beta = beta,
    .noise_amp = noise_amp,
    .mode = mode,
    .Te = Te,
    .Ti = Ti,
    .vti = vti,
    .vAe = vAe,
    .B0 = B0,
    .vAi = vAi,
    .omega_ci = omega_ci,
    .omega_ce = omega_ce,
    .larmor_ion = larmor_ion,
    .larmor_elc = larmor_elc,
    .l = l,
    .nb = n0 * nb_over_n0,
    .Te_frac = Te_frac,
    .Ti_frac = Ti_frac,
    .ix = ix,
    .iy = iy,
    .Nx = Nx,
    .Ny = Ny,
    .Lx = Lx,
    .Ly = Ly,
    .k0_elc = k0_elc,
    .k0_ion = k0_ion,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct lhdi_grad_closure_ctx *app = ctx;

  double pi = app -> pi;

  double mass_elc = app -> mass_elc;
  double charge_elc = app -> charge_elc;

  double n0 = app -> n0;

  double noise_amp = app -> noise_amp;
  double mode = app -> mode;
  
  double Te = app -> Te;

  double B0 = app -> B0;
  
  double l = app -> l;

  double Te_frac = app -> Te_frac;
  double ix = app -> ix;
  double iy = app -> iy;

  double Lx = app -> Lx;
  double Ly = app -> Ly;

  double sech_sq = (1.0 / cosh(y / l)) * (1.0 / cosh(y / l)); // Hyperbolic secant squared.

  double n = n0 * sech_sq; // Total number density.
  double Jx_noise = -noise_amp * (iy * pi / Ly) * sin(iy * pi * y / Ly) * sin(ix * 2.0 * pi* x / Lx) / mode; // Current density noise (x-direction).
  double Jy_noise = -noise_amp * (ix * 2.0 * pi / Lx) * cos(iy * pi * y / Ly) * cos(ix * 2.0 * pi * x / Lx) / mode; // Current density noise (y-direction).

  double Jx  = (B0 / l) * (-sech_sq) + Jx_noise; // Total current density, with noise (x-direction).
  double Jy  = Jy_noise; // Total current density, with noise (y-direction).

  double rhoe = n * mass_elc; // Electron mass density.
  double momxe = (mass_elc / charge_elc) * Jx * Te_frac; // Electron momentum density (x-direction).
  double momye = (mass_elc / charge_elc) * Jy * Te_frac; // Electron momentum density (y-direction).
  double pre = n * Te; // Electron pressure (scalar).

  // Set electron mass density.
  fout[0] = rhoe;
  // Set electron momentum density.
  fout[1] = momxe; fout[2] = momye; fout[3] = 0.0;
  // Set electron pressure tensor.
  fout[4] = pre + momxe * momxe / rhoe; fout[5] = momxe * momye / rhoe; fout[6] = 0.0;  
  fout[7] = pre + momye * momye / rhoe; fout[8] = 0.0; fout[9] = pre;  
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct lhdi_grad_closure_ctx *app = ctx;

  double pi = app -> pi;

  double mass_ion = app -> mass_ion;
  double charge_ion = app -> charge_ion;

  double n0 = app -> n0;

  double noise_amp = app -> noise_amp;
  double mode = app -> mode;
  
  double Ti = app -> Ti;

  double B0 = app -> B0;
  
  double l = app -> l;

  double Ti_frac = app -> Ti_frac;
  double ix = app -> ix;
  double iy = app -> iy;

  double Lx = app -> Lx;
  double Ly = app -> Ly;

  double sech_sq = (1.0 / cosh(y / l)) * (1.0 / cosh(y / l)); // Hyperbolic secant squared.

  double n = n0 * sech_sq; // Total number density.
  double Jx_noise = -noise_amp * (iy * pi / Ly) * sin(iy * pi * y / Ly) * sin(ix * 2.0 * pi* x / Lx) / mode; // Current density noise (x-direction).
  double Jy_noise = -noise_amp * (ix * 2.0 * pi / Lx) * cos(iy * pi * y / Ly) * cos(ix * 2.0 * pi * x / Lx) / mode; // Current density noise (y-direction).

  double Jx  = (B0 / l) * (-sech_sq) + Jx_noise; // Total current density, with noise (x-direction).
  double Jy  = Jy_noise; // Total current density, with noise (y-direction).

  double rhoi = n * mass_ion; // Ion mass density.
  double momxi = (mass_ion / charge_ion) * Jx * Ti_frac; // Ion momentum density (x-direction).
  double momyi = (mass_ion / charge_ion) * Jy * Ti_frac; // Ion momentum density (y-direction).
  double pri = n * Ti; // Ion pressure (scalar).

  // Set ion mass density.
  fout[0] = rhoi;
  // Set ion momentum density.
  fout[1] = momxi; fout[2] = momyi; fout[3] = 0.0;
  // Set ion pressure tensor.
  fout[4] = pri + momxi * momxi / rhoi; fout[5] = momxi * momyi / rhoi; fout[6] = 0.0;  
  fout[7] = pri + momyi * momyi / rhoi; fout[8] = 0.0; fout[9] = pri;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct lhdi_grad_closure_ctx *app = ctx;

  double pi = app -> pi;

  double noise_amp = app -> noise_amp;
  double mode = app -> mode;

  double B0 = app -> B0;

  double l = app -> l;

  double ix = app -> ix;
  double iy = app -> iy;
  
  double Lx = app -> Lx;
  double Ly = app -> Ly;

  double Bz_noise = noise_amp * cos(iy * pi * y / Ly) * sin(ix * 2.0 * pi * x / Lx) / mode;
  double Bx = 0.0;
  double By = 0.0;
  double Bz = -B0 * tanh(y / l) + Bz_noise;

  // Set electric field.
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // Set magnetic field.
  fout[3] = Bx, fout[4] = By; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    gkyl_moment_app_write(app, t_curr, iot -> curr - 1);
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

  struct lhdi_grad_closure_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);  

  // Electron/ion equations.
  struct gkyl_wv_eqn *elc_ten_moment = gkyl_wv_ten_moment_new(ctx.k0_elc);
  struct gkyl_wv_eqn *ion_ten_moment = gkyl_wv_ten_moment_new(ctx.k0_ion);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .equation = elc_ten_moment,
    .has_grad_closure = true, // Include gradient closure.
    .evolve = true,
    .init = evalElcInit,
    .ctx = &ctx,

    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
  };

  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .equation = ion_ten_moment,
    .has_grad_closure = true, // Include gradient closure.
    .evolve = true,
    .init = evalIonInit,
    .ctx = &ctx,

    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },    
  };

  // Field.
  struct gkyl_moment_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .mag_error_speed_fact = 1.0,
    
    .evolve = true,
    .init = evalFieldInit,
    .ctx = &ctx,
    
    .bcy = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
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
    .name = "10m_lhdi_grad_closure",

    .ndim = 2,
    .lower = { -0.5 * ctx.Lx, -0.5 * ctx.Ly},
    .upper = { 0.5 * ctx.Lx, 0.5 * ctx.Ly},
    .cells = { NX, NY },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },
    .cfl_frac = ctx.cfl_frac,

    .num_species = 2,
    .species = { elc, ion },

    .field = field,

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

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // Initialize simulation.
  gkyl_moment_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr);

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

    write_data(&io_trig, app, t_curr);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (dt < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_moment_app_cout(app, stdout, "WARNING: Time-step dt = %g", dt);
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

  write_data(&io_trig, app, t_curr);
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
  gkyl_wv_eqn_release(elc_ten_moment);
  gkyl_wv_eqn_release(ion_ten_moment);
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
