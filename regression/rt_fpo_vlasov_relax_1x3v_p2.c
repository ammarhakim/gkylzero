#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct free_stream_ctx {
  double charge; // charge
  double mass; // mass
  double vt; // thermal velocity
  double Lx; // size of the box
  double gamma; // collision frequency factor in FPO
  double t_end; // end time of simulation
  int num_frame; // number of frames to write out
};

static inline double sq(double x) { return x*x; }

static inline double
bump_maxwellian(double n, double vx, double vy, double vz, double ux, double uy, double uz, double vt, double bA, double bUx, double bUy, double bUz, double bS, double bVt)
{
  double v2 = (vx - ux)*(vx - ux) + (vy - uy)*(vy - uy) + (vz - uz)*(vz - uz);
  double bv2 = (vx - bUx)*(vx - bUx) + (vy - bUy)*(vy - bUy) + (vz - bUz)*(vz - bUz);
  return n/pow(sqrt(2*M_PI*vt*vt), 3)*exp(-v2/(2*vt*vt)) + n/pow(sqrt(2*M_PI*bVt*bVt), 3)*exp(-bv2/(2*bVt*bVt))*(bA*bA)/(bv2 + bS*bS);
}

void
evalDistFuncSquare(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct free_stream_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  double square_fac = 4.0*app->vt;
  if(vx>-square_fac && vx<square_fac && vy>-square_fac && vy<square_fac && vz>-square_fac && vz<square_fac) {
    fout[0] = 0.5;
  } else {
    fout[0] = 0.0;
  }
}

void
evalDistFuncBump(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct free_stream_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  fout[0] = bump_maxwellian(1.0, vx, vy, vz, 0.0, 0.0, 0.0, app->vt, 
    sqrt(0.1), 4*sqrt(0.25/3), 0.0, 0.0, 0.12, 3.0*app->vt);
}

void
evalGamma(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct free_stream_ctx *app = ctx;  
  fout[0] = app->gamma;
}

struct free_stream_ctx
create_ctx(void)
{
  struct free_stream_ctx ctx = {
    .mass = 1.0,
    .charge = 0.0,
    .vt = 1.0,
    .Lx = 1.0,
    .gamma = 1.0, 
    .t_end = 1.0, 
    .num_frame = 1, 
  };
  return ctx;
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
  struct free_stream_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 2);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 16); 

  struct gkyl_vlasov_species square = {
    .name = "square",
    .charge = ctx.charge, .mass = ctx.mass,
    .lower = { -8.0*ctx.vt, -8.0*ctx.vt, -8.0*ctx.vt },
    .upper = { 8.0*ctx.vt, 8.0*ctx.vt, 8.0*ctx.vt }, 
    .cells = { NV, NV, NV },

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFuncSquare,
      .ctx_func = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_FPO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalGamma,
    },
    
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2" },
  };

  struct gkyl_vlasov_species bump = {
    .name = "bump",
    .charge = ctx.charge, .mass = ctx.mass,
    .lower = { -5.0*ctx.vt, -5.0*ctx.vt, -5.0*ctx.vt },
    .upper = { 5.0*ctx.vt, 5.0*ctx.vt, 5.0*ctx.vt }, 
    .cells = { NV, NV, NV },

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFuncBump,
      .ctx_func = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_FPO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalGamma,
    },
    
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2" },
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "fpo_vlasov_relax_1x3v_p2",

    .cdim = 1, .vdim = 3,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { NX },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { square, bump },
    .skip_field = true,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.t_end;
  double dt = tend-tcurr;
  int nframe = ctx.num_frame;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_vlasov_app_calc_integrated_mom(app, tcurr); 

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    gkyl_vlasov_app_calc_integrated_mom(app, tcurr);   
    
    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_vlasov_app_stat_write(app);
  gkyl_vlasov_app_write_integrated_mom(app);

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
  printf("Species collisions took %g secs\n", stat.species_coll_mom_tm);  
  printf("Species collisions took %g secs\n", stat.species_coll_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
