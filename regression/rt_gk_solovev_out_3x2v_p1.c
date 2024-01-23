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

struct gk_solovev_ctx {
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
  // Source parameters
  double lambda_source;
  double x_source;
              
              
  // Simulation parameters
  double Lx; // Box size in x
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
  .filepath = "./efit_data/solovev.geqdsk",
  .rzpoly_order = 2,
  .fluxpoly_order = 1,
  .plate_spec = false,
  .quad_param = {  .eps = 1e-10 }
};

struct gkyl_tok_geo_grid_inp ginp = {
  .ftype = GKYL_SOL_DN_OUT,
  .rclose = 3.0, // any number larger than ~2 will do
  .zmin = -1.5,
  .zmax = 1.5,
  .write_node_coord_array = true,
  .node_file_nm = "solovev_nodes.gkyl"
}; 


void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double Lz = app->Lz;
  double Ls = Lz/4;
  double floor = 0.1;


  // find source density at z = 0
  double source_density = 0;
  double source_floor = 1e-10;
  if (x > x_source)
     source_floor = 1e-2; // higher floor to left of source peak
  source_density = fmax(exp(-(x-x_source)*(x-x_source)/((2*lambda_source)*(2*lambda_source))), source_floor);

  // find source temp at z = 0
  double source_temp = 0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  if (x > x_source - 3*lambda_source)
     source_temp = 80*eV;
  else
     source_temp = 30*eV;

  // now compute initial desity
  double effective_source = 3.800419e+23*fmax(source_density, floor);
  double c_ss = sqrt(5.0/3.0*source_temp/app->massIon);
  double n_peak = 4*sqrt(5)/3/c_ss*Ls*effective_source/2;
  double perturb = 0;
  if (fabs(z) <= Ls)
     fout[0]  = n_peak*(1+sqrt(1-(z/Ls)*(z/Ls)))/2*(1+perturb);
  else
     fout[0] = n_peak/2*(1+perturb);
}

void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double eV = GKYL_ELEMENTARY_CHARGE;
  if (x > x_source - 3*lambda_source)
     fout[0] = 50*eV;
  else
     fout[0] = 20*eV;
}

void
eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double eV = GKYL_ELEMENTARY_CHARGE;
  if (x > x_source - 3*lambda_source)
     fout[0] = 50*eV;
  else
     fout[0] = 20*eV;
}

void
eval_density_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double Lz = app->Lz;
  double source_floor = 1e-10;
  if (x > x_source)
     source_floor = 1e-2; // higher floor to left of source peak
  if (fabs(z) < Lz/4)
     fout[0] = 3.800419e+23*fmax(exp(-(x-x_source)*(x-x_source)/((2*lambda_source)*(2*lambda_source))), source_floor);
  else
     fout[0] = 3.800419e+23*1e-40;
}

void
eval_upar_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double eV = GKYL_ELEMENTARY_CHARGE;
  if (x > x_source - 3*lambda_source)
     fout[0] = 80*eV;
  else
     fout[0] = 30*eV;
}

void
eval_temp_ion_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];
  double lambda_source = app->lambda_source;
  double x_source = app->x_source;
  double eV = GKYL_ELEMENTARY_CHARGE;
  if (x > x_source - 3*lambda_source)
     fout[0] = 80*eV;
  else
     fout[0] = 30*eV;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *app = ctx;
  fout[0] = app->nuIon;
}

void
bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_solovev_ctx *gc = ctx;
  fout[0] = gc->B0;
}

struct gk_solovev_ctx
create_ctx(void)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double Te = 40.0*eV;
  double Ti = 40.0*eV;
  double B0 = 0.55; // Magnetic field magnitude in Tesla
  double n0 = 7.0e18; // Particle density in 1/m^3

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

  // Simulation box size (psi, alpha, theta).
  double q0 = 2.0;
  double r0 = 0.40705706492831;
  double Lx = 0.02;
  double Ly = 50 * rho_s * q0/r0 ; // should be 0.107890816895
  double Lz = 3.14*2;

  // Source parameters.
  double x_source = -0.07;
  double lambda_source = 0.03*Lx;

  // Velocity Grid
  double vpar_max_elc = 4.0*vtElc;
  double mu_max_elc = 0.75*me*(4.0*vtElc)*(4.0*vtElc)/(2.0*B0);

  double vpar_max_ion = 4.0*vtIon;
  double mu_max_ion = 0.75*mi*(4.0*vtIon)*(4.0*vtIon)/(2.0*B0);

  double finalTime = 1.0e-6; 
  double numFrames = 1;

  struct gk_solovev_ctx ctx = {
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
    .lambda_source = lambda_source,
    .x_source = x_source,
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

  struct gk_solovev_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 36);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 1);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], 12);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 16);
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
      .proj_id = GKYL_PROJ_MAXWELLIAN, 
      .ctx_density = &ctx,
      .density = eval_density,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_elc,      
    },
    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .projection = {
        .proj_id = GKYL_PROJ_MAXWELLIAN, 
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_elc_source,      
      }, 
    },
    .bcx = { GKYL_SPECIES_ZERO_FLUX, GKYL_SPECIES_ZERO_FLUX },
    .bcz = { GKYL_SPECIES_GK_SHEATH, GKYL_SPECIES_GK_SHEATH },
        
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
      .proj_id = GKYL_PROJ_MAXWELLIAN, 
      .ctx_density = &ctx,
      .density = eval_density,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ion,      
    },
    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .projection = {
        .proj_id = GKYL_PROJ_MAXWELLIAN, 
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_ion_source,      
      }, 
    },
    .bcx = { GKYL_SPECIES_ZERO_FLUX, GKYL_SPECIES_ZERO_FLUX },
    .bcz = { GKYL_SPECIES_GK_SHEATH, GKYL_SPECIES_GK_SHEATH },
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
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
    .name = "gk_solovev_out_3x2v_p1",

    .cdim = 3, .vdim = 2,
    .lower = { -0.08, -ctx.Ly/2.0, -ctx.Lz/2.0 },
    .upper = { -0.06, ctx.Ly/2.0, ctx.Lz/2.0 },
    .cells = { NX, NY, NZ },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_TOKAMAK,
      .tok_efit_info = &inp,
      .tok_grid_info = &ginp,
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
