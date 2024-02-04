#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <rt_arg_parse.h>

struct gk_lbo_ctx {
  double B0; // reference magnetic field
  double n0; // reference density
  double u0;
  double vt;
  double nu;
  double ab;
  double ub;
  double sb;
  double vtb;

  double finalTime;
  int numFrames; // number of output frames
};


void
eval_tophat(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lbo_ctx *app = ctx;
  double n0 = app->n0, u0 = app->u0, vt = app->vt;
  double x = xn[0], v = xn[1], mu = xn[2];
  double v0 = sqrt(3)*vt;
  if (fabs(v) < v0)
    fout[0] = n0/2.0/v0;
  else
    fout[0] = 0.0;
}

void
eval_bump(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lbo_ctx *app = ctx;
  double n0 = app->n0, u0 = app->u0, vt = app->vt;
  double x = xn[0], v = xn[1], mu = xn[2];
  double B0 = app->B0;
  double ab = app->ab, ub = app->ub, sb = app->sb, vtb = app->vtb;
  double vsq = ((v-u0)/(sqrt(2.0)*vt))*((v-u0)/(sqrt(2.0)*vt)) + mu*B0;
  double vsqb = ((v-u0)/(sqrt(2.0)*vtb))*((v-u0)/(sqrt(2.0)*vtb)) + mu*B0;
  fout[0] = (n0/sqrt(2.0*M_PI*vt))*exp(-vsq) + (n0/sqrt(2.0*M_PI*vtb))*exp(-vsqb)*(ab*ab)/((v-ub)*(v-ub)+sb*sb);
}

void
evalNu(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lbo_ctx *app = ctx;
  fout[0] = app->nu;
}


void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void
bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lbo_ctx *gc = ctx;
  fout[0] = gc->B0;
}

struct gk_lbo_ctx
create_ctx(void)
{
  double finalTime = 100.0; 
  double numFrames = 1;

  double n0 = 1.0;
  double u0 = 0.0;
  double vt = 1.0/3.0;
  double nu = 0.01;
  double B0 = 1.0;
  // bump params
  double ab   = sqrt(0.1);                      // Amplitude of bump.
  double ub   = 4*sqrt( ((3*vt/2)*(3*vt/2))/3);         // Location of bump.
  double sb   = 0.12;                                // Softening factor to avoid divergence.
  double vtb  = 1.0;                                 // Thermal speed of Maxwellian in bump.

  struct gk_lbo_ctx ctx = {
    .n0 = n0,
    .u0 = u0,
    .B0 = B0,
    .vt = vt,
    .nu = nu,
    .ab = ab,
    .ub = ub,
    .sb = sb,
    .vtb = vtb,
    .finalTime = finalTime, 
    .numFrames = numFrames, 
  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_gyrokinetic_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_gyrokinetic_app_write(app, tcurr, iot->curr-1);
    gkyl_gyrokinetic_app_calc_mom(app); gkyl_gyrokinetic_app_write_mom(app, tcurr, iot->curr-1);
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

  struct gk_lbo_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 2);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 32);
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], 16);

  // electrons
  struct gkyl_gyrokinetic_species square = {
    .name = "square",
    .charge = 1.0, .mass = 1.0,
    .lower = { -8*ctx.vt, 0.0},
    .upper = { 8*ctx.vt, 12.0*ctx.vt*ctx.vt/2/ctx.B0}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0,

    .ctx_dist= &ctx,
    .init_dist = eval_tophat,

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .normNu = true,
      .bmag_mid = ctx.B0,
      .norm_nu_facs = {0.0, ctx.nu*sqrt(0.11404*0.11404*0.11404)/1.01036},
      .self_nu = evalNu,
      .num_cross_collisions = 1,
      .collide_with = { "bump" },
    },


    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // ions
  struct gkyl_gyrokinetic_species bump = {
    .name = "bump",
    .charge = 1.0, .mass = 1.0,
    .lower = { -8*ctx.vt, 0.0},
    .upper = { 8*ctx.vt, 12.0*ctx.vt*ctx.vt/2/ctx.B0}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0,


    .ctx_dist= &ctx,
    .init_dist = eval_bump,

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .normNu = true,
      .bmag_mid = ctx.B0,
      .norm_nu_facs = {0.0, ctx.nu*sqrt(0.39677*0.39677*0.39677)/1.10187},
      .self_nu = evalNu,
      .num_cross_collisions = 1,
      .collide_with = { "square" },
    },

    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ADIABATIC,
    .electron_mass = 1,
    .electron_charge = 1,
    .electron_temp = ctx.vt,
    .bmag_fac = ctx.B0, 
    .fem_parbc = GKYL_FEM_PARPROJ_NONE, 
  };

  // GK app
  struct gkyl_gk gk = {
    .name = "gk_lborelax_varnu_1x2v_p1",

    .cdim = 1, .vdim = 2,
    .lower = { 0.0},
    .upper = { 1.0},
    .cells = { NX},
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = &ctx,
      .bmag_func = bmag_func, // mapping of computational to physical space
      .bmag_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = {0},
    .skip_field = true,

    .num_species = 2,
    .species = { square, bump },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.finalTime;
  double dt = tend-tcurr;
  int nframe = ctx.numFrames;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_gyrokinetic_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    if (step % 100 == 0) {
      gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);
    }
    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }
  gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);
  gkyl_gyrokinetic_app_write_field_energy(app);
  gkyl_gyrokinetic_app_stat_write(app);
  
  // fetch simulation statistics
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  // simulation complete, free app
  gkyl_gyrokinetic_app_release(app);
  
  return 0;
}
