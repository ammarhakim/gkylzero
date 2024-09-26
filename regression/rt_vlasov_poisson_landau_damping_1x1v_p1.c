#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct vp_langmuir_ctx {
  int cdim, vdim; // Dimensionality.
  double epsilon_0; // Permittivity of vacuum.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double n_elc; // Number density.
  double temp_elc; // Electron temperature.

  double knumber; // Perturbation wavenumber.
  double amplitude; // Perturbation amplitude

  double Lx; // Domain length
  double x_min; // Minimum x coordinate.
  double x_max; // Maximum x coordinate.
  double vx_min; // Minimum vx coordinate.
  double vx_max; // Maximum vx coordinate.
  int Nx; // Number of cells in x.
  int Nvx; // Number of cells in vx.
  int poly_order; // Polynomial order of the basis.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

static inline double
maxwellian(const double *v, double n, const double *u, double vt, int vdim)
{
  // v: velocity coordinate vector.
  // n: density.
  // u: mean flow speed (vector).
  // vt: thermal speed.
  double vsq = 0;
  for (int d=0; d<vdim; d++)
    vsq += pow(v[d] - u[d],2);
  return n/sqrt(pow(2.0*M_PI*pow(vt,2),vdim))*exp(-vsq/(2.0*pow(vt,2)));
}

void
eval_distf_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  const double *v = &xn[1];

  struct vp_langmuir_ctx *app = ctx;
  int vdim = app->vdim;
  double n = app->n_elc, vt = sqrt(app->temp_elc/app->mass_elc);
  double udrift[vdim];
  for (int d=0; d<vdim; d++)
    udrift[d] = 0.0;
  double k = app->knumber, delta = app->amplitude;

  fout[0] = (1.0 + delta*cos(k*x)) * maxwellian(v, n, udrift, vt, vdim);
}

struct vp_langmuir_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 1; // Dimensionality.

  double epsilon_0 = 1.0; // Permittivity of vacuum.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.

  double n_elc = 1.0; // Number density.
  double temp_elc = 1.0; // Electron temperature.

  // Thermal speed, plasma frequency and Debye length.
  double vt_elc = sqrt(temp_elc/mass_elc);
  double omega_pe = sqrt(pow(charge_elc,2)*n_elc/(epsilon_0*mass_elc));
  double lambda_D = vt_elc/omega_pe;

  double knumber = 0.5/lambda_D; // Perturbation wavenumber.
  double amplitude = 1e-4; // Perturbation amplitude

  double Lx = 2.*M_PI/knumber; // Domain length
  double x_min = -Lx/2.0; // Minimum x coordinate.
  double x_max =  Lx/2.0; // Maximum x coordinate.
  double vx_min = -6.0*vt_elc; // Minimum vx coordinate.
  double vx_max =  6.0*vt_elc; // Maximum vx coordinate.
  int Nx = 64; // Number of cells in x.
  int Nvx = 64; // Number of cells in vx.
  int poly_order = 1; // Polynomial order of the basis.

  double t_end = 15.0/omega_pe; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct vp_langmuir_ctx ctx = {
    .cdim = cdim, .vdim = vdim,
    .epsilon_0         = epsilon_0        ,
    .mass_elc          = mass_elc         , 
    .charge_elc        = charge_elc       , 
    .n_elc             = n_elc            , 
    .temp_elc          = temp_elc         , 
    .knumber           = knumber          , 
    .amplitude         = amplitude        , 
    .Lx                = Lx               , 
    .x_min             = x_min            , 
    .x_max             = x_max            , 
    .vx_min            = vx_min           , 
    .vx_max            = vx_max           , 
    .Nx                = Nx               , 
    .Nvx               = Nvx              , 
    .poly_order        = poly_order       ,
    .t_end             = t_end            , 
    .num_frames        = num_frames       , 
    .int_diag_calc_num = int_diag_calc_num, 
    .dt_failure_tol    = dt_failure_tol   , 
    .num_failures_max  = num_failures_max , 
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

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct vp_langmuir_ctx ctx = create_ctx(); // Context for init functions.

  // Electrons.
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { ctx.vx_min},
    .upper = { ctx.vx_max}, 
    .cells = { ctx.Nvx },

    .num_init = 1,
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = eval_distf_elc,
      .ctx_func = &ctx,
    },

    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2" },
  };

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon_0,
    .poisson_bcs = {
      .lo_type = { GKYL_POISSON_PERIODIC },
      .up_type = { GKYL_POISSON_PERIODIC },
    },
  };

  // VP app.
  struct gkyl_vm vm = {
    .name = "vp_landau_damping_1x1v_p1",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .lower = { ctx.x_min },
    .upper = { ctx.x_max },
    .cells = { ctx.Nx },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .cfl_frac = 0.6,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 1,
    .species = { elc },
    .field = field,
    .is_electrostatic = true,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
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
  gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  freeresources:
  // Free resources after simulation completion.
  gkyl_vlasov_app_release(app);
  
  return 0;
}
