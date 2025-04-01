#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_math.h>

#include <rt_arg_parse.h>

// Define the context of the simulation. This is basically all the globals
struct gk_app_ctx {
  int cdim, vdim; // Dimensionality.
  
  // Geometry and magnetic field.
  double a_shift;   // Parameter in Shafranov shift.
  double Z_axis;    // Magnetic axis height [m].
  double R_axis;    // Magnetic axis major radius [m].

  double R0;        // Major radius of the simulation box [m].
  double a_mid;     // Minor radius at outboard midplane [m].
  double r0;        // Minor radius of the simulation box [m].
  double B0;        // Magnetic field magnitude in the simulation box [T].
  double kappa;     // Elongation (=1 for no elongation).
  double delta;     // Triangularity (=0 for no triangularity).

  double x_LCFS;    // Radial location of the last closed flux surface.

  // Plasma parameters.
  double me;  double qe;
  double mi;  double qi;
  double n0;  double Te0;  double Ti0; 

  // Collisions.
  double nuFrac;  double nuElc;  double nuIon;

  // Source parameters.
  double n_srcOMP;        // Amplitude of the OMP source
  double x_srcOMP;        // Radial location of the OMP source.
  double Te_srcOMP;       // Te for the OMP source.
  double Ti_srcOMP;       // Ti for the OMP source.
  double sigma_srcOMP;    // Radial spread of the OMP source.
  double n_srcGB;         // Amplitude of the grad-B source
  double x_srcGB;         // Radial location of the grad-B source.
  double sigma_srcGB;     // Radial spread of the grad-B source.
  double bfac_srcGB;      // Field aligned spread of the grad-B source. 
  double Te_srcGB;        // Te for the grad-B source.
  double Ti_srcGB;        // Ti for the grad-B source.
  double floor_src;       // Source floor.

  // Grid parameters.
  double Lx;        // Domain size in radial direction.
  double Lz;        // Domain size along magnetic field.
  double x_min;  double x_max;
  double z_min;  double z_max;
  int Nx;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;
  double vpar_max_elc;  double mu_max_elc;
  double vpar_max_ion;  double mu_max_ion;

  double t_end;   int num_frames;
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};


// Common source profiles.
void density_srcOMP(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_app_ctx *app = ctx;
  double n_srcOMP = app->n_srcOMP;
  double x_srcOMP = app->x_srcOMP;
  double sigma_srcOMP = app->sigma_srcOMP;
  double floor_src = app->floor_src;

  fout[0] = n_srcOMP*(exp(-(pow(x-x_srcOMP,2))/(2.*pow(sigma_srcOMP,2)))+floor_src);
}
void zero_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

// Electron source profiles.
void density_elc_srcGB(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_app_ctx *app = ctx;
  double n_srcGB = app->n_srcGB;
  double x_srcGB = app->x_srcGB;
  double sigma_srcGB = app->sigma_srcGB;
  double bfac_srcGB = app->bfac_srcGB;
  double floor_src = app->floor_src;

  fout[0] = n_srcGB*exp(-pow(x-x_srcGB,2)/(2.*pow(sigma_srcGB,2)))
           *GKYL_MAX2(sin(z)*exp(-pow(fabs(z),1.5)/(2*pow(bfac_srcGB,2))),0.);
}
void temp_elc_srcOMP(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_app_ctx *app = ctx;
  double x_srcOMP = app->x_srcOMP;
  double sigma_srcOMP = app->sigma_srcOMP;
  double Te_srcOMP = app->Te_srcOMP;

  if (x < x_srcOMP + 3*sigma_srcOMP) {
    fout[0] = Te_srcOMP;
  } else {
    fout[0] = Te_srcOMP*3./8.;
  }
}
void temp_elc_srcGB(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_app_ctx *app = ctx;
  double Te_srcGB = app->Te_srcGB;

  fout[0] = Te_srcGB;
}

// Ion source profiles.
void density_ion_srcGB(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_app_ctx *app = ctx;
  double n_srcGB = app->n_srcGB;
  double x_srcGB = app->x_srcGB;
  double sigma_srcGB = app->sigma_srcGB;
  double bfac_srcGB = app->bfac_srcGB;
  double floor_src = app->floor_src;

  fout[0] = -n_srcGB*exp(-pow(x-x_srcGB,2)/(2.*pow(sigma_srcGB,2)))
           *GKYL_MAX2(sin(z)*exp(-pow(fabs(z),1.5)/(2*pow(bfac_srcGB,2))),0.);
}
void temp_ion_srcOMP(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_app_ctx *app = ctx;
  double x_srcOMP = app->x_srcOMP;
  double sigma_srcOMP = app->sigma_srcOMP;
  double Ti_srcOMP = app->Ti_srcOMP;

  if (x < x_srcOMP + 3*sigma_srcOMP) {
    fout[0] = Ti_srcOMP;
  } else {
    fout[0] = Ti_srcOMP*3./8.;
  }
}
void temp_ion_srcGB(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_app_ctx *app = ctx;
  double Ti_srcGB = app->Ti_srcGB;

  fout[0] = Ti_srcGB;
}

// Ion initial conditions
void density_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_app_ctx *app = ctx;
  double n0 = app->n0;

  fout[0] = n0*1.1*(0.5*(1.+tanh(2.*(2.-25.*x)))+0.01);
}
void temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_app_ctx *app = ctx;
  double Ti0 = app->Ti0;

  fout[0] = Ti0*((1./3.)*(2.+tanh(2.*(2.-25.*x)))+0.01);
}

// Electron initial conditions
void temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_app_ctx *app = ctx;
  double Te0 = app->Te0;

  fout[0] = Te0*((1./3.)*(2.+tanh(2.*(2.-25.*x)))+0.01);
}

// Collision frequencies.
void nuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuElc;
}
void nuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuIon;
}

struct gk_app_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Geometry and magnetic field.
  double a_shift   = 0.6;                // Parameter in Shafranov shift.
  double Z_axis    = 0.00232616113;      // Magnetic axis height [m].
  double R_axis    = 1.72068012;         // Magnetic axis major radius [m].
  double B_axis    = 2.0;                // Magnetic field at the magnetic axis [T].
  double R_LCFSmid = 2.2801477223421736; // Major radius of the LCFS at the outboard midplane [m].
  double Rmid_min  = R_LCFSmid - 0.1;    // Minimum midplane major radius of simulation box [m].
  double Rmid_max  = R_LCFSmid + 0.05;   // Maximum midplane major radius of simulation box [m].
  double R0        = 0.5*(Rmid_min+Rmid_max);  // Major radius of the simulation box [m].

  double a_mid     = R_LCFSmid-R_axis;   // Minor radius at outboard midplane [m].
  // Redefine a_mid with Shafranov shift, to ensure LCFS radial location.
  a_mid = R_axis/a_shift - sqrt(R_axis*(R_axis - 2*a_shift*R_LCFSmid + 2*a_shift*R_axis))/a_shift;

  double r0        = R0-R_axis;          // Minor radius of the simulation box [m].
  double B0        = B_axis*(R_axis/R0); // Magnetic field magnitude in the simulation box [T].
  double kappa     = 1.488;              // Elongation (=1 for no elongation).
  double delta     = 0.1;                // Triangularity (=0 for no triangularity).

  // Grid parameters
  int Nx = 4;
  int Nz = 16;
  int poly_order = 1;
  double x_LCFS = 0.036;     // Keep in mind (x_max - x_LCFS)/dx must be an integer
  double x_max = 0.048;
  double x_min = x_max - 2.0*(x_max-x_LCFS);
  double Lx = x_max - x_min;  // Domain size along x.

  double Lz        = 2.*(M_PI-1e-14);      // Domain size along magnetic field.
  double z_min     = -Lz/2.;
  double z_max     =  Lz/2.;



  // Plasma parameters. Chosen based on the value of a cubic sline
  // between the last TS data inside the LCFS and the probe data in
  // in the far SOL, near R=0.475 m.
  double AMU = 2.01410177811;
  double mi  = mp*AMU;   // Deuterium ions.
  double Te0 = 200*eV;
  double Ti0 = 200*eV;
  double n0  = 2.0e19;   // [1/m^3]

  double nuFrac = 1.0;
  // Electron-electron collision freq.
  double logLambdaElc = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuElc = nuFrac * logLambdaElc * pow(eV, 4) * n0 /
    (6*sqrt(2.) * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(me) * pow(Te0,3./2.));
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4) * n0 /
    (12 * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(mi) * pow(Ti0,3./2.));

  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi); // Thermal speeds.

  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Source parameters
  double n_srcOMP = 3.e19;
  double x_srcOMP = x_min;
  double Te_srcOMP = 2*Te0;
  double Ti_srcOMP = 2*Ti0;
  double sigma_srcOMP = 0.03*Lx;
  double n_srcGB = 0.04*1.3553306021427228e+23;
  double x_srcGB = x_min;
  double sigma_srcGB = 10*rho_s;
  double bfac_srcGB = 1.2;
  double Te_srcGB = 350*eV;
  double Ti_srcGB = 350*eV;
  double floor_src = 1e-2;

  // Velocity  Grid parameters
  int Nvpar = 8;
  int Nmu = 4;
  double vpar_max_elc = 4.*vte;
  double mu_max_elc = me*pow(4*vte,2)/(2*B0);
  double vpar_max_ion = 4.*vti;
  double mu_max_ion = mi*pow(4*vti,2)/(2*B0);

  double t_end = 1.e-6;
  int num_frames = 10;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .a_shift = a_shift,
    .R_axis = R_axis,
    .R0     = R0    ,
    .a_mid  = a_mid ,
    .r0     = r0    ,
    .B0     = B0    ,
    .kappa  = kappa ,
    .delta  = delta ,
    .Lx     = Lx    ,
    .Lz     = Lz    ,
    .x_min = x_min,  .x_max = x_max,
    .z_min = z_min,  .z_max = z_max,

    .x_LCFS = x_LCFS,
  
    .me = me,  .qe = qe,
    .mi = mi,  .qi = qi,
    .n0 = n0,  .Te0 = Te0,  .Ti0 = Ti0,
  
    .nuFrac = nuFrac,  .nuElc = nuElc,  .nuIon = nuIon,
  
    .n_srcOMP     = n_srcOMP    ,
    .x_srcOMP     = x_srcOMP    ,
    .Te_srcOMP    = Te_srcOMP   ,
    .Ti_srcOMP    = Ti_srcOMP   ,
    .sigma_srcOMP = sigma_srcOMP,
    .n_srcGB      = n_srcGB     ,
    .x_srcGB      = x_srcGB     ,
    .sigma_srcGB  = sigma_srcGB ,
    .bfac_srcGB   = bfac_srcGB  ,
    .Te_srcGB     = Te_srcGB    ,
    .Ti_srcGB     = Ti_srcGB    ,
    .floor_src    = floor_src   ,
  
    .Nx = Nx,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .vpar_max_elc = vpar_max_elc,  .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,  .mu_max_ion = mu_max_ion,

    .t_end = t_end,  .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
  }
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_gyrokinetic_app_write(app, t_curr, frame);

    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_write_field_energy(app);

    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
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

  struct gk_app_ctx ctx = create_ctx(); // context for init functions

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // electrons
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe, .mass = ctx.me,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
      .density = density_init,
      .upar = zero_func,
      .temp = temp_elc,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = nuElc,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 2,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = density_srcOMP,
        .upar = zero_func,
        .temp = temp_elc_srcOMP,
      },
      .projection[1] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = density_elc_srcGB,
        .upar = zero_func,
        .temp = temp_elc_srcGB,
      },
    },

    .diffusion = {
      .num_diff_dir = 1,
      .diff_dirs = { 0 },
      .D = { 0.5 },
      .order = 2,
    },

    .bcx = {
      .lower={.type = GKYL_SPECIES_ABSORB,},
      .upper={.type = GKYL_SPECIES_ABSORB,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_IWL,
              .aux_parameter = ctx.x_LCFS,},
      .upper={.type = GKYL_SPECIES_GK_IWL,
              .aux_parameter = ctx.x_LCFS,},
    },

    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
      .density = density_init,
      .upar = zero_func,
      .temp = temp_ion,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = nuIon,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 2,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = density_srcOMP,
        .upar = zero_func,
        .temp = temp_ion_srcOMP,
      },
      .projection[1] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = density_ion_srcGB,
        .upar = zero_func,
        .temp = temp_ion_srcGB,
      },
    },

    .diffusion = {
      .num_diff_dir = 1,
      .diff_dirs = { 0 },
      .D = { 0.5 },
      .order = 2,
    },

    .bcx = {
      .lower={.type = GKYL_SPECIES_ABSORB,},
      .upper={.type = GKYL_SPECIES_ABSORB,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_IWL,
              .aux_parameter = ctx.x_LCFS,},
      .upper={.type = GKYL_SPECIES_GK_IWL,
              .aux_parameter = ctx.x_LCFS,},
    },

    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES_IWL,
    .xLCFS = ctx.x_LCFS,
    .poisson_bcs = {.lo_type = {GKYL_POISSON_DIRICHLET},
                    .up_type = {GKYL_POISSON_DIRICHLET},
                    .lo_value = {0.0}, .up_value = {0.0}},
  };

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .filepath = "./data/eqdsk/g171646.02800", // equilibrium to use
    .rz_poly_order = 2,                       // polynomial order for psi(R,Z) used for field line tracing
    .flux_poly_order = 1,                     // polynomial order for fpol(psi)
    .reflect = true,                          // Reflect lower half of psi(R,Z) for up-down symmetry
  };

  struct gkyl_tok_geo_grid_inp grid_inp = {
    .ftype = GKYL_IWL ,     // type of geometry
    .rclose = 2.5,
    .rright = 2.5,                // Closest R to outboard SOL
    .rleft = 0.0,                 // closest R to inboard SOL
    .rmin = 0.5,                  // smallest R in machine
    .rmax = 2.5,                  // largest R in machine
    .zmin = -0.83,                 // Z of upper plate
    .zmax = 0.83,                  // Z of lower plate
    .use_cubics = false,          // Whether to use cubic representation of psi(R,Z) for field line tracing
  };

  // GK app
  struct gkyl_gk gk = {
    .name = "gk_d3d_iwl_numerical_2x2v_p1",

    .cfl_frac_omegaH = 1.0,
    .cfl_frac = 1.0,

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .lower = { ctx.x_min, ctx.z_min },
    .upper = { ctx.x_max, ctx.z_max },
    .cells = { cells_x[0], cells_x[1] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {
      .world = {0.0},
      .geometry_id = GKYL_TOKAMAK,
      .efit_info = efit_inp,
      .tok_grid_info = grid_inp,
    },

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };


  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);

  // Initial and final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, t_curr);
  }  

  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
  write_data(&trig_write, app, t_curr, false);

  // Initial time-step.
  double dt = t_end-t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  printf("Starting main loop ...\n");
  long step = 1, num_steps = app_args.num_steps;
  while ((t_curr < t_end) && (step <= num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", t_curr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, t_curr > t_end);
    write_data(&trig_write, app, t_curr, t_curr > t_end);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, true);
        write_data(&trig_write, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }
  printf(" ... finished\n");

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

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.n_io);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  freeresources:
  // simulation complete, free app
  gkyl_gyrokinetic_app_release(app);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
