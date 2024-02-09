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
  double massAr; // Argon mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double TAr; // Argon temperature
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double B0; // reference magnetic field
  double n0; // reference density
  double nsource;
  // Source parameters
  double T_source; // Source electron temperature
  double cx;
  double cz;
  // Simulation parameters
  double Lz; // Box size in z
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double vpar_max_Ar; // Velocity space extents in vparallel for Li ions
  double mu_max_Ar; // Velocity space extents in mu for Li ions    
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
    .zmin = -5.14213,
    .zmax = 5.14226,
    .write_node_coord_array = true,
    .node_file_nm = "step_outboard_fixed_z_nodes.gkyl"
  };

void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double n0 = app->n0;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = n0*exp(-(x-xcenter)*(x-xcenter)/(2.0*cx*cx)) * exp(-z*z/(2.0*cz*cz));
  if (n/n0 < 1.0e-5)
    n = n0*1.0e-5;
  fout[0] = n;
}

void
eval_density_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double n0 = app->n0;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = 0.99*n0*exp(-(x-xcenter)*(x-xcenter)/(2.0*cx*cx)) * exp(-z*z/(2.0*cz*cz));
  if (n/n0 < 1.0e-5)
    n = n0*1.0e-5;
  fout[0] = n;
}

void
eval_density_ar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double n0 = app->n0;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = 0.0;
  fout[0] = n;
}

void
eval_density_ar0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double n0 = 0.01*app->n0;
  double xmax = 1.509820; // taken from gkyl output file. CHECK THIS!
  double lambda_ar_x = .1;
  double lambda_ar_z = 0.25;
  double zmin = -3.14; 
  double xcenter = 1.2014;
  // edit this profile appropriately
  double n; 
  if (z <= 0) {
    n = n0*exp((x-xmax)/lambda_ar_x) * exp(-(z-zmin)/lambda_ar_z);
  }
  else {
    n = n0*exp((x-xmax)/lambda_ar_x) * exp((z+zmin)/lambda_ar_z);
  }
  if (n/n0 < 1e-5)
    n = n0*1e-5;
  fout[0] = n;
}

void
eval_vm_max(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1], vx = xn[2], vy = xn[3], vz = xn[4];
  double n0 = 0.01*app->n0;
  double vt2Ar = app->TAr/app->massAr;
  double xmax = 1.509820; // taken from gkyl output file. CHECK THIS!
  double lambda_ar_x = .1;
  double lambda_ar_z = 0.25;
  double zmin = -3.14; 
  double xcenter = 1.2014;
  // edit this profile appropriately
  double n; 
  if (z <= 0.0) {
    n = n0*exp((x-xmax)/lambda_ar_x) * exp(-(z-zmin)/lambda_ar_z);
  }
  else {
    n = n0*exp((x-xmax)/lambda_ar_x) * exp((z+zmin)/lambda_ar_z);
  }
  if (n/n0 < 1.0e-5)
    n = n0*1.0e-5;
  
  fout[0] = n/pow(2*M_PI*vt2Ar,1.5)*exp(-(vx*vx + vy*vy* + vz*vz)/(2.0*vt2Ar));
}

void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; 
}

void
eval_udrift(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; 
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->Te;
  fout[0] = T;
}

void
eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->Ti;
  fout[0] = T;
}

void
eval_temp_ar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->TAr;
  fout[0] = T;
}

void
eval_density_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double nsource = app->nsource;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = nsource*exp(-(x-xcenter)*(x-xcenter)/(2.0*cx*cx)) * exp(-z*z/(2.0*cz*cz));
  if (n/nsource < 1.0e-5)
    n = nsource*1.0e-5;
  fout[0] = n;
}

void
eval_upar_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double T = app->T_source;
  fout[0] = T;
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

struct gk_step_ctx
create_ctx(void)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // D ion mass
  double me = GKYL_ELECTRON_MASS;
  double mAr = 39.95*GKYL_PROTON_MASS; // Ar ion mass
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double Te = 100*2.8*eV;
  double Ti = 150*2.8*eV;
  double TAr = 30*eV; 
  double B0 = 2.51; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19/2.8; // Particle density in 1/m^3

  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double vtAr = sqrt(TAr/mAr);

  // Source parameters.
  double nsource = 3.9e23/2.8; // peak source rate in particles/m^3/s 
  double T_source = 285*eV*2.8;
  double cx = 0.0065612*9;
  double cz = 0.4916200;

  // Collision parameters.
  double nuFrac = 0.25;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).
  double Lz = 3.14*2;

  double vpar_max_elc = 8.0*vtElc;
  double mu_max_elc = 12*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 8.0*vtIon;
  double mu_max_ion = 12*mi*vtIon*vtIon/(2.0*B0);

  double vpar_max_Ar = 8.0*vtAr;
  double mu_max_Ar = 12*mi*vtAr*vtAr/(2.0*B0);

  double finalTime = 1.0e-8; 
  double numFrames = 1;

  struct gk_step_ctx ctx = {
    .chargeElc = qe, 
    .massElc = me,
    .chargeIon = qi, 
    .massIon = mi,
    .massAr = mAr,
    .Te = Te, 
    .Ti = Ti,
    .TAr = TAr,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .B0 = B0, 
    .n0 = n0, 
    .T_source = T_source, 
    .nsource = nsource,
    .cx = cx,
    .cz = cz,
    .Lz = Lz, 
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion,
    .vpar_max_Ar = vpar_max_Ar, 
    .mu_max_Ar = mu_max_Ar,
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 8);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], 16);
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
      .upar = eval_upar,
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
        .temp = eval_temp_source,      
      }, 
    },
    .react_neut = { 
      .num_react = 2, 
      .react_type = { 
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ELC, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar1", // ion is always the higher charge state  
          .donor_nm = "Ar0", // interacts with elc to give up charge 
          .charge_state = 0, // corresponds to lower charge state (donor) 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB, 
          .type_self = GKYL_SELF_ELC, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar1", 
          .recvr_nm = "Ar0", 
          .charge_state = 0, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
      },
    },  
    //.react= { 
    //  .num_react = 2, 
    //  .react_type = { 
    //	  { .react_id = GKYL_REACT_IZ, 
    //      .type_self = GKYL_SELF_ELC, 
    //      .ion_id = GKYL_ION_AR, 
    //      .elc_nm = "elc", 
    //      .ion_nm = "Ar2", // ion is always the higher charge state  
    //      .donor_nm = "Ar1", // interacts with elc to give up charge 
    //      .charge_state = 1, // corresponds to lower charge state (donor) 
    //      .ion_mass = ctx.massAr, 
    //      .elc_mass = ctx.massElc, 
    //    }, 
    //    { .react_id = GKYL_REACT_RECOMB, 
    //      .type_self = GKYL_SELF_ELC, 
    //      .ion_id = GKYL_ION_AR, 
    //      .elc_nm = "elc", 
    //      .ion_nm = "Ar2", 
    //      .recvr_nm = "Ar1", 
    //      .charge_state = 1, 
    //      .ion_mass = ctx.massAr, 
    //      .elc_mass = ctx.massElc, 
    //    }, 
    //  },
    //},  
    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 
    .bcx = { GKYL_SPECIES_ZERO_FLUX, GKYL_SPECIES_ZERO_FLUX },
    .bcy = { GKYL_SPECIES_GK_SHEATH, GKYL_SPECIES_GK_SHEATH },
    
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
    .polarization_density = 0.99*ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN, 
      .ctx_density = &ctx,
      .density = eval_density,
      .ctx_upar = &ctx,
      .upar = eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
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
        .upar = eval_upar_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_source,      
      }, 
    },
    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 
    .bcx = { GKYL_SPECIES_ZERO_FLUX, GKYL_SPECIES_ZERO_FLUX },
    .bcy = { GKYL_SPECIES_GK_SHEATH, GKYL_SPECIES_GK_SHEATH },
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // Ar1+ ions
  struct gkyl_gyrokinetic_species Ar1 = {
    .name = "Ar1",
    .charge = ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = 0.01*ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN, 
      .ctx_density = &ctx,
      .density = eval_density_ar,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },
    .react_neut = {
      .num_react = 2, 
      .react_type = { 
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ION, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar1", 
          .donor_nm = "Ar0", 
          .charge_state = 0, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB, 
          .type_self = GKYL_SELF_ION, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar1", 
          .recvr_nm = "Ar0", 
          .charge_state = 0, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
      }, 
    }, 
    .bcx = { GKYL_SPECIES_ZERO_FLUX, GKYL_SPECIES_ZERO_FLUX },
    .bcy = { GKYL_SPECIES_GK_SHEATH, GKYL_SPECIES_GK_SHEATH },
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // neutral Ar
  struct gkyl_gyrokinetic_neut_species Ar0 = {
    .name = "Ar0", .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, -ctx.vpar_max_Ar, -ctx.vpar_max_Ar},
    .upper = { ctx.vpar_max_Ar, ctx.vpar_max_Ar, ctx.vpar_max_Ar },
    .cells = { NV, NV, NV},
    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .ctx_func = &ctx,
      .func = eval_vm_max,
    },

    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ZERO_FLUX },
    .bcy = { GKYL_SPECIES_ZERO_FLUX, GKYL_SPECIES_ZERO_FLUX },
    
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2"}, //, "M2par", "M2perp" },
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
    .name = "gk_step_out_ar_2x2v_p1",

    .cdim = 2, .vdim = 2,
    .lower = { 0.934, -ctx.Lz/2.0 },
    .upper = { 1.5098198350000001, ctx.Lz/2.0 },
    .cells = { NX, NZ },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .world = {0.0},
      .geometry_id = GKYL_TOKAMAK,
      .tok_efit_info = &inp,
      .tok_grid_info = &ginp,
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {  },

    .num_species = 3,
    .species = { elc, ion, Ar1},
    .num_neut_species = 1,
    .neut_species = {Ar0},
    
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
