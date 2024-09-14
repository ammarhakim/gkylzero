#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>


static inline double sq(double x) { return x*x; }

struct annulus_sodshock_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  double rhol; // Left/inner density.
  double ul; // Left/inner velocity (r-direction).
  double ul_theta; // Left/inner velocity (theta-direction).
  double templ; // Left/inner temperature.

  double rhor; // Right/outer density.
  double ur; // Right/outer velocity (r-direction).
  double ur_theta; // Right/outer velocity (theta-direction).
  double tempr; // Right/outer temperature.

  // Simulation parameters.
  int Nr; // Cell count (configuration space: r-direction).
  int Ntheta; // Cell count (configuration space: theta-direction).
  int Nv; // Cell count (velocity space: all directions).
  double Rmin; // Domain size (configuration space: r-direction, min value).
  double Rmax; // Domain size (configuration space: r-direction, max value).
  double theta_min; // Domain size (configuration space: theta-direction, min value).
  double theta_max; // Domain size (configuration space: theta-direction, max value).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  double charge; // charge
  double mass; // mass
  double vt; // thermal velocity
  double midplane; // midplane of the theta coord.
};

struct annulus_sodshock_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  double rhol = 1.0; // Left/inner density.
  double ul = 0.0; // Left/inner velocity (r-direction).
  double ul_theta = 0.0; // Left/inner velocity (theta-direction).
  double templ = 1.0; // Left/inner temperature.

  double rhor = 0.125; // Right/outer density.
  double ur = 0.0; // Right/outer velocity (r-direction).
  double ur_theta = 0.0; // Right/outer velocity (theta-direction).
  double tempr = sqrt(0.1 / 0.125); // Right/outer temperature.

  // Simulation parameters.
  int Nr = 64; // Cell count (configuration space: x-direction).
  int Ntheta = 1; // Cell count (configuration space: y-direction).
  int Nv = 12; // Cell count (velocity space: all directions).
  double Rmin = 0.5; // Domain size (configuration space: r-direction, min value).
  double Rmax = 1.5; // Domain size (configuration space: r-direction, max value).
  double theta_min = 0.0; // Domain size (configuration space: theta-direction, min value).
  double theta_max = 2*pi; // Domain size (configuration space: theta-direction, max value).
  int poly_order = 2; // Polynomial order.
  double cfl_frac = 0.9; // CFL coefficient.

  double t_end = 0.1; // Final simulation time.
  int num_frames = 5; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double charge = 1.0; // charge
  double mass = 1.0; // mass
  double vt = 1.0; // thermal velocity
  double midplane = 1.0; // Midplane location in R for the jump in quantities

  struct annulus_sodshock_ctx ctx = {
    .pi = pi,
    .rhol = rhol,
    .ul = ul,
    .ul_theta = ul_theta,
    .templ = templ,
    .rhor = rhor,
    .ur = ur,
    .ur_theta = ur_theta,
    .tempr = tempr,
    .Nr = Nr,
    .Ntheta = Ntheta,
    .Nv = Nv,
    .Rmin = Rmin,
    .Rmax = Rmax,
    .theta_min = theta_min,
    .theta_max = theta_max,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .charge = charge,
    .mass = mass,
    .vt = vt,
    .midplane = midplane,
  };

  return ctx;
}


void
evalNu(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double theta = xn[0], v = xn[1];
  fout[0] = 15000.0;
}

void 
h_ij_inv(double t, const double* xn, double* fout, void* ctx)
{
  // Inverse metric tensor, must be symmetric!
  // [h^{xx},h^{xy},h^{yy}]
  struct annulus_sodshock_ctx *app = (struct annulus_sodshock_ctx *)ctx;
  double q_R = xn[0], q_theta = xn[1];
  const double q[2] = {q_R, q_theta};

  // [h^{thetatheta},h^{thetaphi},h^{phiphi}]
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 1.0 / pow(q[0], 2);
}

void 
det_h(double t, const double* xn, double* fout, void* ctx)
{
  // determinant of the metric tensor: J = det(h_{ij})
  struct annulus_sodshock_ctx *app = (struct annulus_sodshock_ctx *)ctx;
  double q_R = xn[0], q_theta = xn[1];
  const double q[2] = {q_R, q_theta};
  fout[0] = q_R;
}

void 
hamil(double t, const double* xn, double* fout, void* ctx)
{
  // Canonical coordinates:
  double q_R = xn[0], q_theta = xn[1], p_R_dot = xn[2], p_theta_dot = xn[3];
  const double q[2] = {q_R, q_theta};
  const double w[2] = {p_R_dot, p_theta_dot};
  struct annulus_sodshock_ctx *app = (struct annulus_sodshock_ctx *)ctx;
  double *h_inv = malloc(3 * sizeof(double));
  h_ij_inv(t, xn, h_inv, ctx); 
  fout[0] = 0.5 * h_inv[0] * w[0] * w[0] + 
            0.5 * (2.0* h_inv[1] * w[1] * w[0]) + 
            0.5 * h_inv[2] * w[1] * w[1];
  free(h_inv);
}

void
evalDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double R = xn[0], theta = xn[1];
  struct annulus_sodshock_ctx *app = ctx;
  double pi = app -> pi;
  double rhol = app -> rhol;
  double rhor = app -> rhor;
  double midplane = app -> midplane;

  double rho = 0.0;

  if (R < midplane) {
    rho = rhol; // Density (left/inner).
  }
  else {
    rho = rhor; // Density (right/outer).
  }

  // Set the density.
  double det_h_val;
  det_h(t, xn, &det_h_val, ctx); 
  fout[0] = rho*det_h_val;
}

void
evalVDriftInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double R = xn[0], theta = xn[1];
  struct annulus_sodshock_ctx *app = ctx;
  double pi = app -> pi;
  double ul_R = app -> ul;
  double ur_R = app -> ur;
  double ul_theta = app -> ul_theta;
  double ur_theta = app -> ur_theta;
  double midplane = app -> midplane;

  double u_R = 0.0;
  double u_theta = 0.0;


  if (R < midplane) {
    u_R = ul_R; // r-velocity (left/inner).
    u_theta = ul_theta; // theta-velocity (left/inner).
  }
  else {
    u_R = ur_R; // r-velocity (right/outer).
    u_theta = ur_theta; // theta-velocity (right-outer).
  }

  // Set the velocity.
  fout[0] = u_R; 
  fout[1] = u_theta;
}

void
evalTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double R = xn[0], theta = xn[1];
  struct annulus_sodshock_ctx *app = ctx;
  double pi = app -> pi;
  double templ = app -> templ;
  double tempr = app -> tempr;
  double midplane = app -> midplane;

  double rho = 0.0;
  double temperature = 0.0;

  if (R < midplane) {
    temperature = templ; // Temperature (left/inner).
  }
  else {
    temperature = tempr; // Temperature (right/outer).
  }

  // Set temperature
  fout[0] = temperature;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_vlasov_app_write(app, tcurr, iot->curr-1);
    gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, iot->curr-1);
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
  struct annulus_sodshock_ctx ctx = create_ctx(); // context for init functions

  int NR = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nr);
  int NTHETA = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ntheta);

  // electrons
  struct gkyl_vlasov_species neut = {
    .name = "neut",
    .model_id = GKYL_MODEL_CANONICAL_PB,
    .charge = ctx.charge, .mass = ctx.mass,
    .lower = { -8.0*ctx.vt, -8.0*ctx.vt},
    .upper = { 8.0*ctx.vt, 8.0*ctx.vt},
    .cells = { ctx.Nv, ctx.Nv },
    .hamil = hamil,
    .h_ij_inv = h_ij_inv,
    .det_h = det_h,
    .hamil_ctx = &ctx,
    .h_ij_inv_ctx = &ctx,
    .det_h_ctx = &ctx,

    // Reflective boundary condition
    .bcx = {GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT},

    .num_init = 1, 
    .projection[0] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
      .density = evalDensityInit,
      .ctx_density = &ctx,
      .V_drift = evalVDriftInit,
      .ctx_V_drift = &ctx,
      .temp = evalTempInit,
      .ctx_temp = &ctx,
      .correct_all_moms = true, 
    },

    // .collisions =  {
    //   .collision_id = GKYL_BGK_COLLISIONS,

    //   .ctx = &ctx,
    //   .self_nu = evalNu,
    //   .has_implicit_coll_scheme = true,
    //   .correct_all_moms = true, 
    // },

    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "LTEMoments" },
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "can_pb_bgk_surf_annulus_sodshock_p2",

    .cdim = 2, .vdim = 2,
    .lower = { ctx.Rmin, ctx.theta_min },
    .upper = { ctx.Rmax, ctx.theta_max },
    .cells = { NR, NTHETA },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 1,
    .periodic_dirs = {1},

    .num_species = 1,
    .species = { neut },
    .skip_field = true,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.t_end;
  double dt = tend-tcurr;
  int nframe = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);

  // Write the inital timestep
  write_data(&io_trig, app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    gkyl_vlasov_app_write_integrated_energy_f(app);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    write_data(&io_trig, app, tcurr);
    
    step += 1;
  }
  gkyl_vlasov_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  // simulation complete, free app
  gkyl_vlasov_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of forward-Euler calls %ld\n", stat.nfeuler);
  printf("Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    printf("Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    printf("Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  printf("Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  printf("Species RHS calc took %g secs\n", stat.species_rhs_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
