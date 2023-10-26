#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm.h>
#include <rt_arg_parse.h>

struct pkpm_ot_ctx {
  double epsilon0;
  double mu0;
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te_Ti; // electron to ion temperature ratio
  double n0;
  double vAe;
  double B0;
  double beta;
  double vtElc;
  double vtIon;
  double nuElc;
  double nuIon;
  double delta_u0;
  double delta_B0;
  double L;
  double Lz;
  double tend;
  double min_dt;
  bool use_gpu;
};

static inline double
maxwellian(double n, double v, double vth)
{
  double v2 = v*v;
  return n/sqrt(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

void
evalDistFuncElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ot_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], z = xn[2], vx = xn[3];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double Lx = app->L;
  double Ly = app->L;
  double u0x = app->delta_u0;
  double u0y = app->delta_u0;
  double B0x = app->delta_B0;
  double B0y = app->delta_B0;
  
  double fv = maxwellian(app->n0, vx, app->vtElc);
    
  fout[0] = fv;
  fout[1] = app->vtElc*app->vtElc*fv;
}
void
evalDistFuncIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ot_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], z = xn[2], vx = xn[3];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double Lx = app->L;
  double Ly = app->L;
  double u0x = app->delta_u0;
  double u0y = app->delta_u0;
  double B0x = app->delta_B0;
  double B0y = app->delta_B0;

  double fv = maxwellian(app->n0, vx, app->vtIon);
    
  fout[0] = fv;
  fout[1] = app->vtIon*app->vtIon*fv;
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ot_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], z = xn[2];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double mu0 = app->mu0;
  double _2pi = 2.0*M_PI;
  double _4pi = 4.0*M_PI;
  double Lx = app->L;
  double Ly = app->L;
  double Lz = app->Lz;
  double u0x = app->delta_u0;
  double u0y = app->delta_u0;
  double B0x = app->delta_B0;
  double B0y = app->delta_B0;

  double ne = app->n0;

  double Jx = (_2pi/Lz)*0.5*B0y*(cos(_2pi*x/Lx - _2pi*z/Lz) + cos(_2pi*x/Lx + _2pi*z/Lz) ) / mu0;
  Jx = Jx + (_2pi/Lz)*0.5*B0y*(cos(_4pi*x/Lx - _2pi*z/Lz) - cos(_4pi*x/Lx + _2pi*z/Lz)) / mu0;
  double Jy = (_2pi/Lz)*B0x*cos(_2pi*y/Ly - _2pi*z/Lz) / mu0;
  double Jz = (_2pi/Lx)*0.5*B0y*(cos(_2pi*x/Lx - _2pi*z/Lz) - cos(_2pi*x/Lx + _2pi*z/Lz) ) / mu0;
  Jz = Jz + (_4pi/Lx)*0.5*B0y*(cos(_4pi*x/Lx - _2pi*z/Lz) + cos(_4pi*x/Lx + _2pi*z/Lz) ) / mu0;
  Jz = Jz + (_2pi/Ly)*B0x*cos(_2pi*y/Ly - _2pi*z/Lz) / mu0;

  double vdrift_x = -u0x*sin(_2pi*y/Ly - _2pi*z/Lz) - Jx / (qi*ne);
  double vdrift_y = 0.5*u0y*(sin(_2pi*x/Lx + _2pi*z/Lz) + sin(_2pi*x/Lx - _2pi*z/Lz)) - Jy / (qi*ne);
  vdrift_y = vdrift_y + 0.5*u0y*(sin(_4pi*x/Lx - _2pi*z/Lz) - sin(_4pi*x/Lx + _2pi*z/Lz));
  double vdrift_z = -Jz / (qi*ne);

  fout[0] = me*ne*vdrift_x;
  fout[1] = me*ne*vdrift_y;
  fout[2] = me*ne*vdrift_z;
}

void
evalFluidIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ot_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], z = xn[2];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double mu0 = app->mu0;
  double _2pi = 2.0*M_PI;
  double _4pi = 4.0*M_PI;
  double Lx = app->L;
  double Ly = app->L;
  double Lz = app->Lz;
  double u0x = app->delta_u0;
  double u0y = app->delta_u0;
  double B0x = app->delta_B0;
  double B0y = app->delta_B0;

  double ni = app->n0;

  double vdrift_x = -u0x*sin(_2pi*y/Ly - _2pi*z/Lz);
  double vdrift_y = 0.5*u0y*(sin(_2pi*x/Lx + _2pi*z/Lz) + sin(_2pi*x/Lx - _2pi*z/Lz));
  vdrift_y = vdrift_y + 0.5*u0y*(sin(_4pi*x/Lx - _2pi*z/Lz) - sin(_4pi*x/Lx + _2pi*z/Lz));
  double vdrift_z = 0.0;

  fout[0] = mi*ni*vdrift_x;
  fout[1] = mi*ni*vdrift_y;
  fout[2] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ot_ctx *app = ctx;

  double x = xn[0], y = xn[1], z = xn[2];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double mu0 = app->mu0;
  double _2pi = 2.0*M_PI;
  double _4pi = 4.0*M_PI;
  double Lx = app->L;
  double Ly = app->L;
  double Lz = app->Lz;
  double u0x = app->delta_u0;
  double u0y = app->delta_u0;
  double B0x = app->delta_B0;
  double B0y = app->delta_B0;
  double B0 = app->B0;

  double ne = app->n0;
  double ni = app->n0;

  double Jx = (_2pi/Lz)*0.5*B0y*(cos(_2pi*x/Lx - _2pi*z/Lz) + cos(_2pi*x/Lx + _2pi*z/Lz) ) / mu0;
  Jx = Jx + (_2pi/Lz)*0.5*B0y*(cos(_4pi*x/Lx - _2pi*z/Lz) - cos(_4pi*x/Lx + _2pi*z/Lz)) / mu0;
  double Jy = (_2pi/Lz)*B0x*cos(_2pi*y/Ly - _2pi*z/Lz) / mu0;
  double Jz = (_2pi/Lx)*0.5*B0y*(cos(_2pi*x/Lx - _2pi*z/Lz) - cos(_2pi*x/Lx + _2pi*z/Lz) ) / mu0;
  Jz = Jz + (_4pi/Lx)*0.5*B0y*(cos(_4pi*x/Lx - _2pi*z/Lz) + cos(_4pi*x/Lx + _2pi*z/Lz) ) / mu0;
  Jz = Jz + (_2pi/Ly)*B0x*cos(_2pi*y/Ly - _2pi*z/Lz) / mu0;

  double Bx = -B0x*sin(_2pi*y/Ly - _2pi*z/Lz);
  double By = 0.5*B0y*(sin(_2pi*x/Lx - _2pi*z/Lz) - sin(_2pi*x/Lx + _2pi*z/Lz));
  By = By + 0.5*B0y*(sin(_4pi*x/Lx - _2pi*z/Lz) + sin(_4pi*x/Lx + _2pi*z/Lz));
  double Bz = B0;

  // Assumes qi = abs(qe)
  double u_xe = -u0x*sin(_2pi*y/Ly - _2pi*z/Lz) - Jx / (qi*ne);
  double u_ye = 0.5*u0y*(sin(_2pi*x/Lx + _2pi*z/Lz) + sin(_2pi*x/Lx - _2pi*z/Lz)) - Jy / (qi*ne);
  u_ye = u_ye + 0.5*u0y*(sin(_4pi*x/Lx - _2pi*z/Lz) - sin(_4pi*x/Lx + _2pi*z/Lz));
  double u_ze = -Jz / (qi*ne);

  // E = - v_e x B ~  (J - u) x B
  double Ex = - (u_ye*Bz - u_ze*By);
  double Ey = - (u_ze*Bx - u_xe*Bz);
  double Ez = - (u_xe*By - u_ye*Bx);
  
  fout[0] = Ex; fout[1] = Ey, fout[2] = Ez;
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ot_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ot_ctx *app = ctx;
  fout[0] = app->nuIon;
}

struct pkpm_ot_ctx
create_ctx(void)
{
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // pemiability of free space

  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge
  double massIon = 25.0; // ion mass
  double chargeIon = 1.0; // ion charge

  double Te_Ti = 1.0; // ratio of electron to ion temperature
  double n0 = 1.0; // initial number density

  double vAe = 0.5;
  double B0 = vAe*sqrt(mu0*n0*massElc);
  double beta = 0.08;
  double vtElc = vAe*sqrt(beta/2.0);

  // ion velocities
  double vAi = vAe/sqrt(massIon);
  double vtIon = vtElc/sqrt(massIon); //Ti/Te = 1.0

  // ion cyclotron frequency and gyroradius
  double omegaCi = chargeIon*B0/massIon;
  double di = vAi/omegaCi;

  // collision frequencies
  double nuElc = 0.01*omegaCi;
  double nuIon = 0.01*omegaCi/sqrt(massIon);

  // OT initial conditions
  double elong = 5.0;
  double u0x = vAi/elong;
  double u0y = vAi/elong;
  double B0x = B0/elong;
  double B0y = B0/elong;

  // domain size and simulation time
  double L = 20.48*di;
  double Lz = elong*L;
  double tend = 100.0/omegaCi;
  
  struct pkpm_ot_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .chargeElc = chargeElc,
    .massElc = massElc,
    .chargeIon = chargeIon,
    .massIon = massIon,
    .Te_Ti = Te_Ti,
    .n0 = n0,
    .vAe = vAe,
    .B0 = B0,
    .beta = beta,
    .vtElc = vtElc,
    .vtIon = vtIon,
    .nuElc = nuElc,
    .nuIon = nuIon,
    .delta_u0 = u0x,
    .delta_B0 = B0x,
    .L = L,
    .Lz = Lz,
    .tend = tend,
    .min_dt = 1.0e-2, 
  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_pkpm_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) 
    gkyl_pkpm_app_write(app, tcurr, iot->curr-1);
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 16);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 16);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
     
  struct pkpm_ot_ctx ctx = create_ctx(); // context for init functions

  // electrons
  struct gkyl_pkpm_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -6.0 * ctx.vtElc},
    .upper = { 6.0 * ctx.vtElc}, 
    .cells = { VX },

    .ctx_dist = &ctx,
    .ctx_fluid = &ctx,
    .init_dist = evalDistFuncElc,
    .init_fluid = evalFluidElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
    },    

    .diffusion = {.D = 1.0e-3, .order=4},
  };
  
  // ions
  struct gkyl_pkpm_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -6.0 * ctx.vtIon},
    .upper = { 6.0 * ctx.vtIon}, 
    .cells = { VX },

    .ctx_dist = &ctx,
    .ctx_fluid = &ctx,
    .init_dist = evalDistFuncIon,
    .init_fluid = evalFluidIon,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
    },    

    .diffusion = {.D = 1.0e-3, .order=4},
  };

  // field
  struct gkyl_pkpm_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc
  };

  // pkpm app
  struct gkyl_pkpm pkpm = {
    .name = "pkpm_3d_ot_p1",

    .cdim = 3, .vdim = 1,
    .lower = { 0.0, 0.0, 0.0 },
    .upper = { ctx.L, ctx.L, ctx.Lz },
    .cells = { NX, NX, 5*NX },
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    
    .num_periodic_dir = 3,
    .periodic_dirs = { 0, 1, 2 },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&pkpm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.tend;
  double dt = tend-tcurr;
  int nframe = 100;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_pkpm_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_pkpm_app_calc_field_energy(app, tcurr);
  gkyl_pkpm_app_calc_integrated_L2_f(app, tcurr);
  gkyl_pkpm_app_calc_integrated_mom(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_pkpm_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_pkpm_update(app, dt);
    gkyl_pkpm_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    if (step % 100 == 0) {
      gkyl_pkpm_app_calc_field_energy(app, tcurr);
      gkyl_pkpm_app_calc_integrated_L2_f(app, tcurr);
      gkyl_pkpm_app_calc_integrated_mom(app, tcurr);
    }
    if (!status.success) {
      gkyl_pkpm_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    if (status.dt_actual < ctx.min_dt) {
      gkyl_pkpm_app_cout(app, stdout, "** Time step crashing! Aborting simulation and writing out last output ....\n");
      gkyl_pkpm_app_write(app, tcurr, 1000);
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }
  gkyl_pkpm_app_calc_field_energy(app, tcurr);
  gkyl_pkpm_app_calc_integrated_L2_f(app, tcurr);
  gkyl_pkpm_app_calc_integrated_mom(app, tcurr);
  gkyl_pkpm_app_write_field_energy(app);
  gkyl_pkpm_app_write_integrated_L2_f(app);
  gkyl_pkpm_app_write_integrated_mom(app);
  gkyl_pkpm_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_pkpm_stat stat = gkyl_pkpm_app_stat(app);

  gkyl_pkpm_app_cout(app, stdout, "\n");
  gkyl_pkpm_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_pkpm_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_pkpm_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_pkpm_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_pkpm_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_pkpm_app_cout(app, stdout, "Fluid Species RHS calc took %g secs\n", stat.fluid_species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species PKPM Vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_pkpm_app_cout(app, stdout, "EM Variables (bvar) calculation took %g secs\n", stat.field_em_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm);
  gkyl_pkpm_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);

  gkyl_pkpm_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_pkpm_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  // simulation complete, free app
  gkyl_pkpm_app_release(app);
  
  return 0;
}
