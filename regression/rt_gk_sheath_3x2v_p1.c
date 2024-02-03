#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_gyrokinetic.h>
#include <rt_arg_parse.h>

struct gk_app_ctx {
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double c_s; // sound speed
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double B0; // reference magnetic field
  double n0; // reference density
  double Lx; // Box size in x.
  double Ly; // Box size in y.
  double Lz; // Box size in z.
  double n_src; // Source density.
  double T_src; // Source temperature.
  double xmu_src; // Source location in x.
  double xsigma_src; // Source spread in x.
  double R0; // Reference major radius.
  double R; // Major radius.
  double a0; // Reference minor radius.
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double finalTime; // end time
  int numFrames; // number of output frames
};

// Source profiles.
void eval_source_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double n_src = app->n_src;
  double Ls = app->Lz/4.;
  double xmu_src = app->xmu_src;
  double xsigma_src = app->xsigma_src;

  double floor_src = 0.1;

  if (fabs(z) < Ls) {
    fout[0] = GKYL_MAX2(exp(-pow(x-xmu_src,2)/(pow(2*xsigma_src,2))), floor_src);
  } else {
    fout[0] = 1.e-40;
  }
  fout[0] = n_src*fout[0];
}

void eval_source_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void eval_source_temp(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double T_src = app->T_src;
  double xmu_src = app->xmu_src;
  double xsigma_src = app->xsigma_src;

  if (x < xmu_src + 3.*xsigma_src) {
    fout[0] = T_src;
  } else {
    fout[0] = (3./8.)*T_src;
  }
}

// Initial conditions.
void eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double n0 = app->n0;
  double Ls = app->Lz/4.;
  double mi = app->massIon;
  double xnCenterZ[] = {xn[0], xn[1], 0.};
  double n_src[1];
  eval_source_density(t, xnCenterZ, n_src, ctx);
  double T_src[1];
  eval_source_temp(t, xnCenterZ, T_src, ctx);

  double c_ss = sqrt((5./3.)*T_src[0]/mi);
  double nPeak = 4.*sqrt(5)/3./c_ss*Ls/2.*n_src[0];
  if (fabs(z) <= Ls) {
    fout[0] = nPeak*(1.+sqrt(1.-pow(z/Ls,2)))/2;
  } else {
    fout[0] = nPeak/2.;
  }
}

void eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double Te = app->Te;
  double xmu_src = app->xmu_src;
  double xsigma_src = app->xsigma_src;

  fout[0] = Te;
  if (x < xmu_src + 3.*xsigma_src) {
    fout[0] = (5./4.)*Te;
  } else {
    fout[0] = 0.5*Te;
  }
}

void eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double Ti = app->Ti;
  double xmu_src = app->xmu_src;
  double xsigma_src = app->xsigma_src;

  if (x < xmu_src + 3.*xsigma_src) {
    fout[0] = (5./4.)*Ti;
  } else {
    fout[0] = 0.5*Ti;
  }
}

void eval_nuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void eval_nuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuIon;
}

void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double R0 = app->R0;
  double a0 = app->a0;

  double x = xc[0], y = xc[1], z = xc[2];

  double R = x;
  double phi = z/(R0+a0);
  double X = R*cos(phi);
  double Y = R*sin(phi);
  double Z = y;

  xp[0] = X;  xp[1] = Y;  xp[2] = Z;
}

void bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double R0 = app->R0;
  double R = app->R;

  double x = xc[0];

  fout[0] = app->B0*R/x;
}

struct gk_app_ctx
create_ctx(void)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS; // Proton mass.
  double me = GKYL_ELECTRON_MASS; // Electron mass.

  double mi = 2.014*mp; // ion mass
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Reference density and temperature.
  double Te = 40.0*eV;
  double Ti = 40.0*eV;
  double n0 = 7.0e18;

  // Geometry and magnetic field.
  double B_axis = 0.5;
  double R0     = 0.85;
  double a0     = 0.15;
  double R      = R0 + a0;
  double B0     = B_axis*(R0/R);

  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double c_s = sqrt(Te/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Collision parameters.
  double nuFrac = 0.1;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq
  
  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).
  double Lx = 50*rho_s;
  double Ly = 100*rho_s;
  double Lz = 4.;

  // Source parameters.
  double n_src = 1.4690539*3.612270e+23;
  double T_src = 2.*Te;
  double xmu_src = R;
  double xsigma_src = 0.005;

  double vpar_max_elc = 4.0*vtElc;
  double mu_max_elc = (3./2.)*0.5*me*pow(4.0*vtElc,2)/(2.0*B0);

  double vpar_max_ion = 4.0*vtIon;
  double mu_max_ion = (3./2.)*0.5*mi*pow(4.0*vtIon,2)/(2.0*B0);

  double finalTime = .5e-6; 
  double numFrames = 1;

  struct gk_app_ctx ctx = {
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .Te = Te, 
    .Ti = Ti, 
    .c_s = c_s, 
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .B0 = B0, 
    .n0 = n0, 
    .Lx = Lx,
    .Ly = Ly,
    .Lz = Lz,
    .n_src = n_src,
    .T_src = T_src,
    .xmu_src    = xmu_src,
    .xsigma_src = xsigma_src,
    .R0 = R0,
    .a0 = a0,
    .R = R,
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
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

  struct gk_app_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 4);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[0], 1);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[0], 8);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 6);
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], 4);

  // electrons
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0,

    .ctx_density = &ctx,
    .init_density = eval_density,
    .ctx_upar = &ctx,
    .init_upar= eval_upar,
    .ctx_temp = &ctx,
    .init_temp = eval_temp_elc,
    .is_maxwellian = true,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = eval_nuElc,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },

    .source = {
      .source_id = GKYL_MAXWELLIAN_SOURCE,
      .write_source = true,
      .ctx_density = &ctx,
      .density_profile = eval_source_density,
      .ctx_upar = &ctx,
      .upar_profile = eval_source_upar,
      .ctx_temp = &ctx,
      .temp_profile = eval_source_temp,
    },
    
    .bcx = { GKYL_SPECIES_ZERO_FLUX, GKYL_SPECIES_ZERO_FLUX },
    .bcz = { GKYL_SPECIES_GK_SHEATH, GKYL_SPECIES_GK_SHEATH },

    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0,

    .ctx_density = &ctx,
    .init_density = eval_density,
    .ctx_upar = &ctx,
    .init_upar = eval_upar,
    .ctx_temp = &ctx,
    .init_temp = eval_temp_ion,
    .is_maxwellian = true,

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = eval_nuIon,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },

    .source = {
      .source_id = GKYL_MAXWELLIAN_SOURCE,
      .write_source = true,
      .ctx_density = &ctx,
      .density_profile = eval_source_density,
      .ctx_upar = &ctx,
      .upar_profile = eval_source_upar,
      .ctx_temp = &ctx,
      .temp_profile = eval_source_temp,
    },

    .bcx = { GKYL_SPECIES_ZERO_FLUX, GKYL_SPECIES_ZERO_FLUX },
    .bcz = { GKYL_SPECIES_GK_SHEATH, GKYL_SPECIES_GK_SHEATH },
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .bmag_fac = ctx.B0,
    .fem_parbc = GKYL_FEM_PARPROJ_NONE,
    .poisson_bcs = {.lo_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_PERIODIC},
                    .up_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_PERIODIC},
                    .lo_value = {0.0, 0.0}, .up_value = {0.0, 0.0}},
  };

  // GK app
  struct gkyl_gk gk = {
    .name = "gk_sheath_3x2v_p1",

    .cdim = 3, .vdim = 2,
    .lower = { ctx.R-ctx.Lx/2.0, -ctx.Ly/2.0, -ctx.Lz/2.0 },
    .upper = { ctx.R+ctx.Lx/2.0,  ctx.Ly/2.0,  ctx.Lz/2.0 },
    .cells = { NX, NY, NZ },
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
    .periodic_dirs = { 1 },

    .num_species = 2,
    .species = { elc, ion },
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
  gkyl_gyrokinetic_app_calc_integrated_mom(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    if (step % 10 == 0) {
      gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);
      gkyl_gyrokinetic_app_calc_integrated_mom(app, tcurr);
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
  gkyl_gyrokinetic_app_calc_integrated_mom(app, tcurr);
  gkyl_gyrokinetic_app_write_field_energy(app);
  gkyl_gyrokinetic_app_write_integrated_mom(app);
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
