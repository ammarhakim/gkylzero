#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_gyrokinetic.h>

#include <rt_arg_parse.h>

struct cross_relax_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  int cdim, vdim; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double epsilon0; // Permittivity of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.

  double B0; // Reference magnetic field strength (Tesla).
  double n0; // Reference number density (1 / m^3).
  double Te; // Electron temperature.
  double Ti; // Ion temperature.
  double alpha; //Ratio of perpendicular to parallel temperature.
  double TePar;
  double TePerp;
  double TiPar;
  double TiPerp;
  double TeIC; // IC elctron temperature.
  double TiIC; // IC ion temperature.

  double log_lambda_elc; // Electron Coulomb logarithm.
  double log_lambda_ion; // Ion Coulomb logarithm.
  double nu_elc; // Electron collision frequency.
  double nu_ion; // Ion collision frequency.

  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.
  double vteIC; // IC electron thermal velocity.
  double vtiIC; // IC ion thermal velocity.
  double uPare; // Electron drift velocity.
  double uPari; // Ion drift velocity.

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double vpar_max_elc; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc; // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion; // Domain boundary (ion velocity space: magnetic moment direction).

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct cross_relax_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  int cdim = 1, vdim = 2; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double epsilon0 = GKYL_EPSILON0; // Permittivity of free space.
  double mass_elc = GKYL_ELECTRON_MASS; // Electron mass.
  double charge_elc = -GKYL_ELEMENTARY_CHARGE; // Electron charge.
  double mass_ion = 2.014 * GKYL_PROTON_MASS; // Deuterium ion mass.
  double charge_ion = GKYL_ELEMENTARY_CHARGE; // Proton charge.

  double B0 = 1.0; // Reference magnetic field strength (Tesla).
  double n0 = 7.0e19; //  Reference number density (1 / m^3).
  double Te = 300.0 * GKYL_ELEMENTARY_CHARGE; // Electron temperature.
  double Ti = 200.0 * GKYL_ELEMENTARY_CHARGE; // Ion temperature.
  double alpha = 1.3; //Ratio of perpendicular to parallel temperature.
  double TePar = Te;
  double TePerp = alpha*Te;
  double TiPar = Ti;
  double TiPerp = alpha*Ti;
  double TeIC = (2.0*TePerp+TePar)/3.0;
  double TiIC = (2.0*TiPerp+TiPar)/3.0;

  // Coulomb logarithms.
  double log_lambda_elc = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Te / charge_ion);
  double log_lambda_ion = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Ti / charge_ion);

  // Collision frequencies.
  double nu_elc = log_lambda_elc * pow(charge_ion,4) * n0 /
    (6.0 * sqrt(2.0) * pow(pi,3.0/2.0) * pow(epsilon0,2) * sqrt(mass_elc) * pow(Te,3.0/2.0));
  double nu_ion = log_lambda_ion * pow(charge_ion,4) * n0 /
    (12.0 * pow(pi,3.0/2.0) * pow(epsilon0,2) * sqrt(mass_ion) * pow(Ti,3.0/2.0));
  
  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.
  double vteIC = sqrt(TeIC / mass_elc); // IC electron thermal velocity.
  double vtiIC = sqrt(TiIC / mass_ion); // IC ion thermal velocity.
  double uPare = 0.5*sqrt(mass_elc/mass_ion)*vteIC; // Electron drift velocity.
  double uPari = 50.0*(mass_elc/mass_ion)*vtiIC; // Ion drift velocity.

  // Simulation parameters.
  int Nx = 1; // Cell count (configuration space: x-direction).
  int Nvpar = 16; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 16; // Cell count (velocity space: magnetic moment direction).
  double vpar_max_elc = 5.0 * vte; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc = mass_elc * pow(5.0*vte,2) / (2.0*B0); // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion = 5.0 * vti; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion = mass_ion * pow(5.0*vti,2) / (2.0*B0); // Domain boundary (ion velocity space: magnetic moment direction).

  double t_end = 0.1/nu_ion; // Final simulation time.
  int num_frames = 10; // Number of output frames.
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct cross_relax_ctx ctx = {
    .pi = pi,
    .cdim = cdim,
    .vdim = vdim,
    .epsilon0 = epsilon0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .B0 = B0,
    .n0 = n0,
    .Te = Te,
    .Ti = Ti, 
    .alpha = alpha, 
    .TePar = TePar,
    .TePerp = TePerp,
    .TiPar = TiPar,
    .TiPerp = TiPerp,
    .TeIC = TeIC, 
    .TiIC = TiIC, 
    .log_lambda_elc = log_lambda_elc,
    .nu_elc = nu_elc,
    .log_lambda_ion = log_lambda_ion,
    .nu_ion = nu_ion,
    .vte = vte,
    .vti = vti,
    .vteIC = vteIC,
    .vtiIC = vtiIC,
    .uPare = uPare,
    .uPari = uPari,
    .Nx = Nx,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nvpar, Nmu},
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .t_end = t_end,
    .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
eval_density_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;
  double n = app->n0;
  fout[0] = n;
}
void
eval_upar_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;
  double uPar = app->uPare;
  fout[0] = uPar;
}
void
eval_temp_par_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;
  double Tpar = app->TePar;
  fout[0] = Tpar;
}
void
eval_temp_perp_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;
  double Tperp = app->TePerp;
  fout[0] = Tperp;
}
void
eval_density_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;
  double n = app->n0;
  fout[0] = n;
}
void
eval_upar_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;
  double uPar = app->uPari;
  fout[0] = uPar;
}
void
eval_temp_par_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;
  double Tpar = app->TiPar;
  fout[0] = Tpar;
}
void
eval_temp_perp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct cross_relax_ctx *app = ctx;
  double Tperp = app->TiPerp;
  fout[0] = Tperp;
}

void
evalNuElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct cross_relax_ctx *app = ctx;

  double nu_elc = app->nu_elc;

  // Set electron collision frequency.
  fout[0] = nu_elc;
}

void
evalNuIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct cross_relax_ctx *app = ctx;

  double nu_ion = app->nu_ion;

  // Set ion collision frequency.
  fout[0] = nu_ion;
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
  struct cross_relax_ctx *app = ctx;

  double B0 = app->B0;

  // Set magnetic field strength.
  fout[0] = B0;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
  }
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_gyrokinetic_app_write(app, t_curr, frame);

    gkyl_gyrokinetic_app_calc_mom(app);
    gkyl_gyrokinetic_app_write_mom(app, t_curr, frame);
    gkyl_gyrokinetic_app_write_source_mom(app, t_curr, frame);

    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_write_field_energy(app);

    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
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

  struct cross_relax_ctx ctx = create_ctx(); // Context for initialization functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Create decomposition.
  struct gkyl_rect_decomp *decomp = gkyl_gyrokinetic_comms_decomp_new(ctx.cdim, cells_x, app_args.cuts, app_args.use_mpi, stderr);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, decomp, stderr);

  int my_rank = 0;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    gkyl_comm_get_rank(comm, &my_rank);
#endif

  // Electron species.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { -ctx.vpar_max_elc, 0.0 },
    .upper = { ctx.vpar_max_elc, ctx.mu_max_elc },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,
    .no_by = true, 

    .projection = {
      .proj_id = GKYL_PROJ_BIMAXWELLIAN,
      .ctx_density = &ctx,
      .density = eval_density_elc,
      .ctx_upar = &ctx,
      .upar= eval_upar_elc,
      .ctx_temppar = &ctx,
      .temppar = eval_temp_par_elc,
      .ctx_tempperp = &ctx,
      .tempperp = eval_temp_perp_elc,
      .correct_all_moms = true,
    },
    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.TeIC, // Temperature used to calculate coulomb logarithm
      .self_nu = evalNuElcInit,
      .ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
      .correct_all_moms = true,
      .use_last_converged = true,
      .iter_eps = 1e-12,
      .max_iter = 10,
    },
    
    .num_diag_moments = 6,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "BiMaxwellianMoments" },
  };

  // Ion species.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .lower = { -ctx.vpar_max_ion, 0.0 },
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,
    .no_by = true, 

    .projection = {
      .proj_id = GKYL_PROJ_BIMAXWELLIAN,
      .ctx_density = &ctx,
      .density = eval_density_ion,
      .ctx_upar = &ctx,
      .upar= eval_upar_ion,
      .ctx_temppar = &ctx,
      .temppar = eval_temp_par_ion,
      .ctx_tempperp = &ctx,
      .tempperp = eval_temp_perp_ion,
      .correct_all_moms = true,
    },
    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.TiIC, // Temperature used to calculate coulomb logarithm
      .self_nu = evalNuIonInit,
      .ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
      .correct_all_moms = true,
      .use_last_converged = true,
      .iter_eps = 1e-12,
      .max_iter = 10,
    },

    .num_diag_moments = 6,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "BiMaxwellianMoments" },
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .fem_parbc = GKYL_FEM_PARPROJ_NONE,
  };

  // GK app.
  struct gkyl_gk app_inp = {
    .name = "gk_bgk_cross_relax_1x2v_p1",

    .cdim = 1, .vdim = 2,
    .lower = { -2.0 }, 
    .upper = { 2.0 }, 
    .cells = { cells_x[0] },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .mapc2p = mapc2p,
      .c2p_ctx = &ctx,
      .bmag_func = bmag_func,
      .bmag_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, ion },
    .skip_field = true, // Do not evolve the field.
    .field = field,

    .use_gpu = app_args.use_gpu,

    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp->ranges[my_rank],
      .comm = comm
    }
  };


  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

  // Initial and final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
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

  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
  write_data(&trig_write, app, t_curr, false);

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

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, t_curr > t_end);
    write_data(&trig_write, app, t_curr, t_curr > t_end);

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
        calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, true);
        write_data(&trig_write, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

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

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  gkyl_gyrokinetic_comms_release(decomp, comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}