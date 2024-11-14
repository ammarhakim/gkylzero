#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct vp_sheath_ctx {
  int cdim, vdim; // Dimensionality.

  double epsilon0; // Permittivity of free space.
  double mass_elc; // Electron mass.
  double mass_ion; // Proton mass.
  double charge_elc; // Electron charge.
  double charge_ion; // Proton charge.

  double n0; // Reference density.
  double Te0; // Reference electron temperature.
  double Ti0; // Reference ion temperature.

  double lambda_D; // Debye length.
  double nu_ee; // Electron-electron collision frequency.
  double nu_ii; // Ion-ion collision frequency.

  double x_min, x_max; // Extents of the x grid.
  double vx_min_elc, vx_max_elc; // Extents of the electron vx grid.
  double vx_min_ion, vx_max_ion; // Extents of the ion vx grid.
  int Nx; // Number of cells along x.
  int Nvx; // Number of cells along vx.
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order; // Polynomial order of the basis.
  double Lx; // Length of the domain along x.

  double L_src; // Source length.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int int_diag_calc_num; // Number of times to compute integrated diagnostics.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

void
eval_dens_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  double n0 = app->n0;
  fout[0] = n0;
}

void
eval_vdrift_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  fout[0] = 0.0;
}

void
eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  double temp = app->Te0;
  fout[0] = temp;
}

void
eval_dens_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  double n0 = app->n0;
  fout[0] = n0;
}

void
eval_vdrift_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  fout[0] = 0.0;
}

void
eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  double temp = app->Ti0;
  fout[0] = temp;
}

void
eval_source_dens_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  double L_src = app->L_src;

  if (fabs(x) < L_src) {
    fout[0] = (L_src - fabs(x))/L_src;
  } else {
    fout[0] = 0.0;
  }
}

void
eval_source_vdrift_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  fout[0] = 0.0;
}

void
eval_source_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  double temp = app->Ti0;

  fout[0] = temp;
}

void
eval_source_dens_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  double L_src = app->L_src;

  if (fabs(x) < L_src) {
    fout[0] = (L_src - fabs(x))/L_src;
  } else {
    fout[0] = 0.0;
  }
}

void
eval_source_vdrift_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  fout[0] = 0.0;
}

void
eval_source_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  double temp = app->Ti0;

  fout[0] = temp;
}

void
eval_nu_elc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  double nu = app->nu_ee;
  double lambda_D = app->lambda_D;

  // Set collision frequency.
  fout[0] = nu/(1 + exp(fabs(x)/(6.0*lambda_D) - 8.0/1.5));
}

void
eval_nu_ion(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], v = xn[1];
  struct vp_sheath_ctx *app = ctx;
  double nu = app->nu_ii;
  double lambda_D = app->lambda_D;

  // Set collision frequency.
  fout[0] = nu/(1 + exp(fabs(x)/(6.0*lambda_D) - 8.0/1.5));
}

struct vp_sheath_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 1; // Dimensionality.

  double epsilon0 = 1.0; // Permittivity of free space.
  double eV = 1.0; // Elementary charge.
  double mass_elc = 1.0; // Electron mass.
  double mass_ion = 1836.153; // Proton mass.
  double charge_elc = -eV; // Electron charge.
  double charge_ion =  eV; // Proton charge.

  double n0 = 1.0; // Reference density.
  double Te0 = 1.0*eV; // Reference electron temperature.
  double Ti0 = Te0; // Reference ion temperature.

  // Thermal speeds.
  double vte = sqrt(Te0/mass_elc);
  double vti = sqrt(Ti0/mass_ion);

  // Plasma frequency and electron Debye length.
  double omega_pe = sqrt((n0 * pow(charge_elc,2.0))/(epsilon0*mass_elc));
  double lambda_D = sqrt((epsilon0 * Te0)/(n0 * pow(charge_elc,2.0)));

  // Collision frequencies.
  double nu_ee = vte/(50.0*lambda_D);
  double nu_ii = vti/(50.0*lambda_D);  

  // Grid extents.
  double Lx = 256.0*lambda_D; // Length of the x domain.
  double x_min = -Lx/2.0, x_max = Lx/2.0;
  double vx_min_elc = -6.0*vte, vx_max_elc = 6.0*vte;
  double vx_min_ion = -6.0*vti, vx_max_ion = 6.0*vti;
  // Number of cells and order of the polynomial basis.
  int Nx = 256;
  int Nvx = 64;
  int poly_order = 2;

  double L_src = 100.0*lambda_D; // Length of the source.

  double t_end = 20.0/omega_pe; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int int_diag_calc_num = num_frames*100; // Number of times to compute integrated diagnostics.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct vp_sheath_ctx ctx = {
    .cdim              = cdim              ,
    .vdim              = vdim              ,
                                           
    .epsilon0          = epsilon0          ,
    .mass_elc          = mass_elc          ,
    .mass_ion          = mass_ion          ,
    .charge_elc        = charge_elc        ,
    .charge_ion        = charge_ion        ,
                                           
    .n0                = n0                ,
    .Te0               = Te0               ,
    .Ti0               = Ti0               ,
                                           
    .lambda_D          = lambda_D          ,
    .nu_ee             = nu_ee             ,
    .nu_ii             = nu_ii             ,  

    .Lx                = Lx                ,
    .x_min             = x_min             ,
    .x_max             = x_max             ,
    .vx_min_elc        = vx_min_elc        ,
    .vx_max_elc        = vx_max_elc        ,
    .vx_min_ion        = vx_min_ion        ,
    .vx_max_ion        = vx_max_ion        ,
    .Nx                = Nx                ,
    .Nvx               = Nvx               ,
    .cells             = {Nx, Nvx}         ,
    .poly_order        = poly_order        ,

    .L_src             = L_src             , 
                                           
    .t_end             = t_end             ,
    .num_frames        = num_frames        ,
    .int_diag_calc_num = int_diag_calc_num ,
    .dt_failure_tol    = dt_failure_tol    ,
    .num_failures_max  = num_failures_max  , 
  };
  return ctx;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  }
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_write)
{
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_vlasov_app_write(app, t_curr, frame);

    gkyl_vlasov_app_calc_mom(app);
    gkyl_vlasov_app_write_mom(app, t_curr, frame);

    gkyl_vlasov_app_calc_field_energy(app, t_curr);
    gkyl_vlasov_app_write_field_energy(app);

    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
    gkyl_vlasov_app_write_integrated_mom(app);
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

  struct vp_sheath_ctx ctx = create_ctx(); // Simulation context.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  struct gkyl_comm *comm = gkyl_vlasov_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { ctx.vx_min_elc},
    .upper = { ctx.vx_max_elc}, 
    .cells = { cells_v[0] },

    .num_init = 1,
    .projection[0] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
      .density = eval_dens_elc,
      .V_drift = eval_vdrift_elc,
      .temp = eval_temp_elc,
      .ctx_density = &ctx,
      .ctx_V_drift = &ctx,
      .ctx_temp = &ctx,
    },

    .source = {
      .source_id = GKYL_BFLUX_SOURCE,
      .source_length = ctx.L_src,
      .source_species = "ion",
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = eval_source_dens_elc,
        .V_drift = eval_source_vdrift_elc,
        .temp = eval_source_temp_elc,
        .ctx_density = &ctx,
        .ctx_V_drift = &ctx,
        .ctx_temp = &ctx,
      },
    },

    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .self_nu = eval_nu_elc,
      .ctx = &ctx,
      .fixed_temp_relax = true,
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ABSORB, },
      .upper = { .type = GKYL_SPECIES_ABSORB, },
    },

    .num_diag_moments = 1,
    .diag_moments = { "LTEMoments" },
  };

  // Ions.
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .lower = { ctx.vx_min_ion},
    .upper = { ctx.vx_max_ion}, 
    .cells = { cells_v[0] },

    .num_init = 1,
    .projection[0] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
      .density = eval_dens_ion,
      .V_drift = eval_vdrift_ion,
      .temp = eval_temp_ion,
      .ctx_density = &ctx,
      .ctx_V_drift = &ctx,
      .ctx_temp = &ctx,
    },

    .source = {
      .source_id = GKYL_BFLUX_SOURCE,
      .source_length = ctx.L_src,
      .source_species = "ion",
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_VLASOV_LTE,
        .density = eval_source_dens_ion,
        .V_drift = eval_source_vdrift_ion,
        .temp = eval_source_temp_ion,
        .ctx_density = &ctx,
        .ctx_V_drift = &ctx,
        .ctx_temp = &ctx,
      },
    },

    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .self_nu = eval_nu_ion,
      .ctx = &ctx,
      .fixed_temp_relax = true,
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ABSORB, },
      .upper = { .type = GKYL_SPECIES_ABSORB, },
    },

    .num_diag_moments = 1,
    .diag_moments = { "LTEMoments" },
  };

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0,
    .poisson_bcs = {
      .lo_type = { GKYL_POISSON_DIRICHLET },
      .up_type = { GKYL_POISSON_DIRICHLET },
      .lo_value = { 0.0 }, .up_value = { 0.0 }
    },
  };

  // VP app
  struct gkyl_vm vm = {
    .name = "vp_sheath_ss_bgk_1x1v_p2",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .lower = { ctx.x_min },
    .upper = { ctx.x_max },
    .cells = { cells_x[0] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,
    .is_electrostatic = true,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm
    }
  };

  // Create app object.
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // Initial & final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_vlasov_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_vlasov_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
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
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
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

      gkyl_vlasov_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_vlasov_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_vlasov_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_vlasov_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_vlasov_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
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

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  freeresources:
  // Free resources after simulation completion.
  gkyl_vlasov_app_release(app);
  gkyl_vlasov_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
