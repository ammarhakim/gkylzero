#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_gr_minkowski.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct quadrants_2d_ctx
{
  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.

  double rho_ul; // Upper left fluid mass density.
  double u_ul; // Upper left fluid x-velocity.
  double v_ul; // Upper left fluid y-velocity.
  double p_ul; // Upper left fluid pressure.

  double rho_ur; // Upper right fluid mass density.
  double u_ur; // Upper right fluid x-velocity.
  double v_ur; // Upper right fluid y-velocity.
  double p_ur; // Upper left fluid pressure.
  
  double rho_ll; // Lower left fluid mass density.
  double u_ll; // Lower left fluid x-velocity.
  double v_ll; // Lower left fluid y-velocity.
  double p_ll; // Lower left fluid pressure.

  double rho_lr; // Lower right fluid mass density.
  double u_lr; // Lower right fluid x-velocity.
  double v_lr; // Lower right fluid y-velocity.
  double p_lr; // Lower right fluid pressure.

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime;

  // Simulation parameters.
  int Nx; // Cell count (x-direction).
  int Ny; // Cell count (y-direction).
  double Lx; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  double loc; // Fluid boundaries (both x and y coordinates).
};

struct quadrants_2d_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.

  double rho_ul = 0.1; // Upper-left fluid mass density.
  double u_ul = 0.99; // Upper-left fluid x-velocity.
  double v_ul = 0.0; // Upper-left fluid y-velocity.
  double p_ul = 1.0; // Upper-left fluid pressure.

  double rho_ur = 0.1; // Upper-right fluid mass density.
  double u_ur = 0.0; // Upper-right fluid x-velocity.
  double v_ur = 0.0; // Upper-right fluid y-velocity.
  double p_ur = 0.01; // Upper-right fluid pressure.
  
  double rho_ll = 0.5; // Lower-left fluid mass density.
  double u_ll = 0.0; // Lower-left fluid x-velocity.
  double v_ll = 0.0; // Lower-left fluid y-velocity.
  double p_ll = 1.0; // Lower-left fluid pressure.

  double rho_lr = 0.1; // Lower-right fluid mass density.
  double u_lr = 0.0; // Lower-right fluid x-velocity.
  double v_lr = 0.99; // Lower-right fluid y-velocity.
  double p_lr = 1.0; // Lower-right fluid pressure.

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  // Simulation parameters.
  int Nx = 800; // Cell count (x-direction).
  int Ny = 800; // Cell count (y-direction).
  double Lx = 1.0; // Domain size (x-direction).
  double Ly = 1.0; // Domain size (y-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 0.4; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double loc = 0.5; // Fluid boundaries (both x and y coordinates).

  struct quadrants_2d_ctx ctx = {
    .gas_gamma = gas_gamma,
    .rho_ul = rho_ul,
    .u_ul = u_ul,
    .v_ul = v_ul,
    .p_ul = p_ul,
    .rho_ur = rho_ur,
    .u_ur = u_ur,
    .v_ur = v_ur,
    .p_ur = p_ur,
    .rho_ll = rho_ll,
    .u_ll = u_ll,
    .v_ll = v_ll,
    .p_ll = p_ll,
    .rho_lr = rho_lr,
    .u_lr = u_lr,
    .v_lr = v_lr,
    .p_lr = p_lr,
    .spacetime = spacetime,
    .Nx = Nx,
    .Ny = Ny,
    .Lx = Lx,
    .Ly = Ly,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .loc = loc,
  };

  return ctx;
}

void
evalGREulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct quadrants_2d_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;

  double rho_ul = app->rho_ul;
  double u_ul = app->u_ul;
  double v_ul = app->v_ul;
  double p_ul = app->p_ul;

  double rho_ur = app->rho_ur;
  double u_ur = app->u_ur;
  double v_ur = app->v_ur;
  double p_ur = app->p_ur;

  double rho_ll = app->rho_ll;
  double u_ll = app->u_ll;
  double v_ll = app->v_ll;
  double p_ll = app->p_ll;

  double rho_lr = app->rho_lr;
  double u_lr = app->u_lr;
  double v_lr = app->v_lr;
  double p_lr = app->p_lr;

  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  double loc = app->loc;

  double rho = 0.0;
  double u = 0.0;
  double v = 0.0;
  double p = 0.0;

  if (y > loc) {
    if (x < loc) {
      rho = rho_ul; // Fluid mass density (upper-left).
      u = u_ul; // Fluid x-velocity (upper-left).
      v = v_ul; // Fluid y-velocity (upper-left).
      p = p_ul; // Fluid pressure (upper-left).
    }
    else {
      rho = rho_ur; // Fluid mass density (upper-right).
      u = u_ur; // Fluid x-velocity (upper-right).
      v = v_ur; // Fluid y-velocity (upper-right).
      p = p_ur; // Fluid pressure (upper-right).
    }
  }
  else {
    if (x < loc) {
      rho = rho_ll; // Fluid mass density (lower-left).
      u = u_ll; // Fluid x-velocity (lower-left).
      v = v_ll; // Fluid y-velocity (lower-left).
      p = p_ll; // Fluid pressure (lower-left).
    }
    else {
      rho = rho_lr; // Fluid mass density (lower-right).
      u = u_lr; // Fluid x-velocity (lower-right).
      v = v_lr; // Fluid y-velocity (lower-right).
      p = p_lr; // Fluid pressure (lower-right).
    }
  }

  double spatial_det, lapse;
  double *shift = malloc(sizeof(double) * 3);
  bool in_excision_region;

  double **spatial_metric = malloc(sizeof(double*) * 3);
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = malloc(sizeof(double) * 3);
  }

  double **inv_spatial_metric = malloc(sizeof(double*) * 3);
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = malloc(sizeof(double) * 3);
  }

  spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
  
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
  spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);

  double *vel = malloc(sizeof(double) * 3);
  double v_sq = 0.0;
  vel[0] = u; vel[1] = v; vel[2] = 0.0;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      v_sq += spatial_metric[i][j] * vel[i] * vel[j];
    }
  }

  double W = 1.0 / (sqrt(1.0 - v_sq));
  if (v_sq > 1.0 - pow(10.0, -8.0)) {
    W = 1.0 / sqrt(1.0 - pow(10.0, -8.0));
  }

  double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));
  
  // Set fluid mass density.
  fout[0] = sqrt(spatial_det) * rho * W;
  // Set fluid momentum density.
  fout[1] = sqrt(spatial_det) * rho * h * (W * W) * u;
  fout[2] = sqrt(spatial_det) * rho * h * (W * W) * v;
  fout[3] = 0.0;
  // Set fluid total energy density.
  fout[4] = sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W));

  // Set spatial metric determinant.
  fout[5] = spatial_det;
  // Set lapse gauge variable.
  fout[6] = lapse;
  // Set shift gauge variables.
  fout[7] = shift[0]; fout[8] = shift[1]; fout[9] = shift[2];

  // Set spatial metric tensor.
  fout[10] = spatial_metric[0][0]; fout[11] = spatial_metric[0][1]; fout[12] = spatial_metric[0][2];
  fout[13] = spatial_metric[1][0]; fout[14] = spatial_metric[1][1]; fout[15] = spatial_metric[1][2];
  fout[16] = spatial_metric[2][0]; fout[17] = spatial_metric[2][1]; fout[18] = spatial_metric[2][2];

  // Set inverse spatial metric tensor.
  fout[19] = inv_spatial_metric[0][0]; fout[20] = inv_spatial_metric[0][1]; fout[21] = inv_spatial_metric[0][2];
  fout[22] = inv_spatial_metric[1][0]; fout[23] = inv_spatial_metric[1][1]; fout[24] = inv_spatial_metric[1][2];
  fout[25] = inv_spatial_metric[2][0]; fout[26] = inv_spatial_metric[2][1]; fout[27] = inv_spatial_metric[2][2];

  // Set excision boundary conditions.
  if (in_excision_region) {
    for (int i = 0; i < 28; i++) {
      fout[i] = 0.0;
    }

    fout[28] = -1.0;
  }
  else {
    fout[28] = 1.0;
  }
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_moment_app_write(app, t_curr, frame);
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

  struct quadrants_2d_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  // Fluid equations.
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_inew(
    &(struct gkyl_wv_gr_euler_inp) {
        .gas_gamma = ctx.gas_gamma,
        .spacetime = ctx.spacetime,
        .rp_type = WV_GR_EULER_RP_LAX,
        .use_gpu = app_args.use_gpu,
    }
  );

  struct gkyl_moment_species fluid = {
    .name = "gr_euler",
    .equation = gr_euler,
    .evolve = true,
    .init = evalGREulerInit,
    .ctx = &ctx,
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
    .name = "gr_quadrants_2d",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { ctx.Lx, ctx.Ly },
    .cells = { NX, NY },

    .scheme_type = GKYL_MOMENT_WAVE_PROP,
    .mp_recon = app_args.mp_recon,

    .cfl_frac = ctx.cfl_frac,

    .num_species = 1,
    .species = { fluid },

    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp->ranges[my_rank],
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
  write_data(&io_trig, app, t_curr, false);

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
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

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

  // Free resources after simulation completion.
  gkyl_wv_eqn_release(gr_euler);
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