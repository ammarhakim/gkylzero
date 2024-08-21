#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>


static inline double sq(double x) { return x*x; }

struct can_pb_ctx {
  double charge; // charge
  double mass; // mass
  double vt; // thermal velocity
  double Lx; // size of the box (q0)
  double Ly; // size of the box (q1)
  double R; // Radius of the sphere
};

void
evalNu(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 15000.0;
}

void 
h_ij_inv(double t, const double* xn, double* fout, void* ctx)
{
  // [h^{xx},h^{xy},h^{yy}]
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 1.0;
}

void 
det_h(double t, const double* xn, double* fout, void* ctx)
{
  fout[0] = 1.0;
}

void 
hamil(double t, const double* xn, double* fout, void* ctx)
{
  // Canonical coordinates:
  double q_x = xn[0], q_y = xn[1], p_x_dot = xn[2], p_y_dot = xn[3];
  const double q[2] = {q_x, q_y};
  const double w[2] = {p_x_dot, p_y_dot};
  struct can_pb_ctx *app = (struct can_pb_ctx *)ctx;
  double R = app->R;
  double *h_inv = malloc(3 * sizeof(double));
  h_ij_inv(t, xn, h_inv, ctx); 
  fout[0] = 0.5 * h_inv[0] * w[0] * w[0] + 
            0.5 * (2.0* h_inv[1] * w[1] * w[0]) + 
            0.5 * h_inv[2] * w[1] * w[1];
  free(h_inv);
}

void
evalDistFuncSquare(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct free_stream_ctx *app = ctx;
  double x = xn[0], y = xn[1], vx = xn[2], vy = xn[3];
  if ((vx>-1.0 && vx<1.0) && (vy>-1.0 && vy<1.0)) {
    fout[0] = 0.5;
  } else {
    fout[0] = 0.0;
  }
}

void
evalDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Inputs Jv*n
  struct can_pb_ctx *app = ctx;
  double x = xn[0];
  double y = xn[1];
  double det_h_val;
  det_h(t, xn, &det_h_val, ctx); 
  fout[0] = 1.0*det_h_val;
}

void
evalVDriftInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = 0.0;
  fout[1] = 0.0;
}

void
evalTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = 1.0;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_vlasov_app_write(app, tcurr, iot->curr-1);
    gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, iot->curr-1);
  }
}


struct can_pb_ctx
create_ctx(void)
{
  struct can_pb_ctx ctx = {
    .mass = 1.0,
    .charge = 1.0,
    .vt = 1.0,
    .Lx = 1.0,
    .Ly = 1.0,
    .R = 1.0,
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct can_pb_ctx ctx = create_ctx(); // context for init functions

  // electrons
  struct gkyl_vlasov_species neut = {
    .name = "neut",
    .model_id = GKYL_MODEL_CANONICAL_PB,
    .charge = ctx.charge, .mass = ctx.mass,
    .lower = { -6.0*ctx.vt, -6.0*ctx.vt},
    .upper = { 6.0*ctx.vt, 6.0*ctx.vt}, 
    .cells = { 32, 32 },
    .hamil = hamil,
    .h_ij_inv = h_ij_inv,
    .det_h = det_h,
    .hamil_ctx = &ctx,
    .h_ij_inv_ctx = &ctx,
    .det_h_ctx = &ctx,

    .num_init = 1, 
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFuncSquare,
      .ctx_func = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNu,
      .correct_all_moms = true, 
    },

    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "LTEMoments" },
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "can_pb_ex_bgk_surf_flat_sq_ic_p2",

    .cdim = 2, .vdim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { ctx.Lx, ctx.Ly },
    .cells = { 2, 2 },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },

    .num_species = 1,
    .species = { neut },
    .skip_field = true,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 0.1;
  double dt = tend-tcurr;
  int nframe = 1;
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  
  write_data(&io_trig, app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
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
