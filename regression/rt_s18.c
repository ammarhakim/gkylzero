#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <rt_arg_parse.h>
#include <gkyl_tok_geo.h>

struct gk_step_ctx {
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double vtIon;
  double vtElc;
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuElcIon; // electron-ion collision frequency
  double nuIonElc; // ion-electron collision frequency
  double nuFrac;
  double B0; // reference magnetic field
  double n0; // reference density
  double nsource;
  // Source parameters
  double Ti_source; // Source electron temperature
  double Te_source;
  double cx;
  double cz;
  // Simulation parameters
  double Ly; // Box size in y
  double Lz; // Box size in z
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double finalTime; // end time
  int numFrames; // number of output frames
};

struct gkyl_tok_geo_efit_inp inp = {
  // psiRZ and related inputs
  .filepath = "./data/eqdsk/step.geqdsk",
  .rzpoly_order = 2,
  .fluxpoly_order = 1,
  .plate_spec = false,
  .quad_param = {  .eps = 1e-10 }
};


struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.2,
    .rright= 6.2,
    .rleft= 2.0,
    .rmin = 1.1,
    .rmax = 6.2,
    .zmin = -8.3,
    .zmax = 8.3,
    .write_node_coord_array = true,
    .node_file_nm = "step_outboard_fixed_z_nodes.gkyl"
  };

void
eval_lab_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double n0 = app->n0;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double vte_sq = app->Te/app->massElc;
  double n = n0*exp(-(x-xcenter)*(x-xcenter)/2/cx/cx) * exp(-z*z/2/cz/cz);
  if (n/n0 < 1e-5)
    n = n0*1e-5;
  fout[0] = n;
  fout[1] = 0.0;
  fout[2] = 3.0*n*vte_sq;
}

void
eval_lab_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double n0 = app->n0;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double vti_sq = app->Ti/app->massIon;
  double n = n0*exp(-(x-xcenter)*(x-xcenter)/2/cx/cx) * exp(-z*z/2/cz/cz);
  if (n/n0 < 1e-5)
    n = n0*1e-5;
  fout[0] = n;
  fout[1] = 0.0;
  fout[2] = 3.0*n*vti_sq;
}

void
eval_lab_source_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double nsource = app->nsource;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double vte_sq = app->Te_source/app->massElc;
  double n = nsource*exp(-(x-xcenter)*(x-xcenter)/2/cx/cx) * exp(-z*z/2/cz/cz);
  if (n/nsource < 1e-5)
    n = nsource*1e-5;
  fout[0] = n;
  fout[1] = 0.0;
  fout[2] = 3.0*n*vte_sq;
}

void
eval_lab_source_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double nsource = app->nsource;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double vti_sq = app->Ti_source/app->massIon;
  double n = nsource*exp(-(x-xcenter)*(x-xcenter)/2/cx/cx) * exp(-z*z/2/cz/cz);
  if (n/nsource < 1e-5)
    n = nsource*1e-5;
  fout[0] = n;
  fout[1] = 0.0;
  fout[2] = 3.0*n*vti_sq;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  fout[0] = app->nuIon;
}

void
bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  fout[0] = app->B0;
}

double plasma_frequency(double n, double m)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  return sqrt(n*eV*eV/m/eps0);
}
double coulomb_log(double ns, double nr, double ms, double mr, double Ts, double Tr, double qs, double qr)
{

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double hbar = GKYL_PLANCKS_CONSTANT_H/2/M_PI;
  double vts = sqrt(Ts/ms);
  double vtr = sqrt(Tr/mr);
  double wps = plasma_frequency(ns,ms);
  double wpr = plasma_frequency(nr,mr);
  double inner1 = wps*wps/(Ts/ms + 3*Ts/ms) + wpr*wpr/(Tr/mr + 3*Ts/ms);
  double u = 3*(vts*vts + vtr*vtr);
  double msr = ms*mr/(ms+mr);
  double inner2 = fmax(fabs(qs*qr)/(4*M_PI*eps0*msr*u*u), hbar/(2*sqrt(eV)*msr*u));
  double inner = (1/inner1)*(1/inner2/inner2) + 1;
  return 0.5*log(inner);
}

double norm_nu_func(double nuFrac, double ns, double nr, double ms, double mr, double qs, double qr, double Ts, double Tr)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double clog = coulomb_log(ns,nr,ms,mr,Ts, Tr, qs, qr);
  double vts = sqrt(Ts/ms);
  double vtr = sqrt(Tr/mr);
  return nuFrac/ms*(1/mr+1/ms)*qs*qs*qr*qr*clog/(6*pow(M_PI,1.5)*eps0*eps0);
}

struct gk_step_ctx
create_ctx(void)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double tempfac = 16.0;
  double Te = 100*2.8*eV/tempfac;
  double Ti = 150*2.8*eV;
  double B0 = 2.51; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19/2.8; // Particle density in 1/m^3

  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(tempfac*Te/me); // use the high temperature for the grid.

  // Source parameters.
  double nsource = 3.9e23/2.8; // peak source rate in particles/m^3/s 
  double Ti_source = 285*eV*2.8;
  double Te_source = Ti_source/tempfac;
  double cx = 0.0065612*2*4;
  double cz = 0.4916200*1.4;

  // Collision parameters.
  double nuFrac = 0.25;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  double nuElcIon = nuElc*sqrt(2);
  double nuIonElc = nuElcIon/(mi/me);


  // Simulation box size (m).
  double Ly = 0.02;
  double Lz = M_PI*2 - 2e-10;

  double vpar_max_elc = 8.0*vtElc;
  double mu_max_elc = 12*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 8.0*vtIon;
  double mu_max_ion = 12*mi*vtIon*vtIon/(2.0*B0);

  double finalTime = 1.0e-3; 
  double numFrames = 200;

  struct gk_step_ctx ctx = {
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .Te = Te, 
    .Ti = Ti, 
    .vtIon = vtIon,
    .vtElc = vtElc,
    .nuElc = nuElc, 
    .nuFrac = nuFrac,
    .nuIon = nuIon, 
    .nuElcIon = nuElcIon,
    .nuIonElc= nuIonElc,
    .B0 = B0, 
    .n0 = n0, 
    .Ti_source = Ti_source, 
    .Te_source = Te_source, 
    .nsource = nsource,
    .cx = cx,
    .cz = cz,
    .Ly = Ly,  
    .Lz = Lz, 
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

  struct gk_step_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 18);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], 16);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 32);
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], 8);

  // electrons
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = { ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_LAB, 
      .ctx_lab_moms = &ctx,
      .lab_moms = eval_lab_elc,    
    },

    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB},
    .bcy = { GKYL_SPECIES_GK_SHEATH, GKYL_SPECIES_GK_SHEATH },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massElc, ctx.massElc, ctx.chargeElc, ctx.chargeElc, ctx.Te, ctx.Te),
      .cross_nu_fac = {norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massElc, ctx.massIon, ctx.chargeElc, ctx.chargeIon, ctx.Te, ctx.Ti)},
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .projection = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_LAB, 
        .ctx_lab_moms = &ctx,
        .lab_moms = eval_lab_source_elc,    
      }, 
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_LAB, 
      .ctx_lab_moms = &ctx,
      .lab_moms = eval_lab_ion,    
    },

    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ABSORB},
    .bcy = { GKYL_SPECIES_GK_SHEATH, GKYL_SPECIES_GK_SHEATH },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massIon, ctx.massIon, ctx.chargeIon, ctx.chargeIon, ctx.Ti, ctx.Ti),
      .cross_nu_fac = {norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massIon, ctx.massElc, ctx.chargeIon, ctx.chargeElc, ctx.Ti, ctx.Te)},
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .projection = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_LAB, 
        .ctx_lab_moms = &ctx,
        .lab_moms = eval_lab_source_ion,    
      }, 
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .bmag_fac = ctx.B0, 
    .fem_parbc = GKYL_FEM_PARPROJ_NONE, 
    .poisson_bcs = {.lo_type = {GKYL_POISSON_DIRICHLET}, 
                    .up_type = {GKYL_POISSON_DIRICHLET}, 
                    .lo_value = {0.0}, .up_value = {0.0}}, 
  };

  // GK app
  struct gkyl_gk gk = {
    .name = "s18",

    .cdim = 2, .vdim = 2,
    .lower = { 0.934, -ctx.Lz/2.0 },
    .upper = { 1.4688, ctx.Lz/2.0 },
    .cells = { NX, NZ },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      //.geometry_id = GKYL_GEOMETRY_FROMFILE,
      .world = {0.0},
      .geometry_id = GKYL_TOKAMAK,
      .tok_efit_info = &inp,
      .tok_grid_info = &ginp,
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {  },

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
