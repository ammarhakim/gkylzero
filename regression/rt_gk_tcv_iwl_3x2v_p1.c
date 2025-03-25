#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <gkyl_math.h>
#include <rt_arg_parse.h>

// Define the context of the simulation. This stores global parameters.
struct gk_app_ctx {
    int cdim, vdim;
    // Geometry and magnetic field parameters
    double a_shift, Z_axis, R_axis, R0, a_mid, x_inner, r0, B0, kappa, delta, q0, Bref, x_LCFS;
    // Plasma parameters
    double me, qe, mi, qi, n0, Te0, Ti0;
    // Collision parameters
    double nuFrac, nuElc, nuIon;
    // Source parameters
    double num_sources;
    double n_srcOMP, x_srcOMP, Te_srcOMP, Ti_srcOMP, sigma_srcOMP, power_src;
    double n_srcRECY, x_srcRECY, Te_srcRECY, Ti_srcRECY, sigmax_srcRECY, sigmaz_srcRECY, floor_src;
    bool adaptive_source;
    char adapt_mom_type [26]; // Type of adaptive moment for the source
    double adapt_mom_rate_target;
    // Grid parameters
    double Lx, Ly, Lz;
    double x_min, x_max, y_min, y_max, z_min, z_max;
    int num_cell_x, num_cell_y, num_cell_z, num_cell_vpar, num_cell_mu;
    int cells[GKYL_MAX_DIM], poly_order;
    double vpar_max_elc, mu_max_elc, vpar_max_ion, mu_max_ion;
    // Simulation control parameters
    double final_time, write_phase_freq;
    int num_frames, int_diag_calc_num, num_failures_max;
    double dt_failure_tol;
};

// Function prototypes (defined at the end of the file)
static double r_x(double x, double a_mid, double x_inner);
static double qprofile(double r, double R_axis);
static void zero_func(double t, const double *xn, double *fout, void *ctx);
static void density_srcOMP(double t, const double *xn, double *fout, void *ctx);
static void temp_elc_srcOMP(double t, const double *xn, double *fout, void *ctx);
static void temp_ion_srcOMP(double t, const double *xn, double *fout, void *ctx);
static void density_srcRECY(double t, const double *xn, double *fout, void *ctx);
static void temp_elc_srcRECY(double t, const double *xn, double *fout, void *ctx);
static void temp_ion_srcRECY(double t, const double *xn, double *fout, void *ctx);
static void density_init(double t, const double *xn, double *fout, void *ctx);
static void temp_elc(double t, const double *xn, double *fout, void *ctx);
static void temp_ion(double t, const double *xn, double *fout, void *ctx);
static void nuElc(double t, const double *xn, double *fout, void *ctx);
static void nuIon(double t, const double *xn, double *fout, void *ctx);
static void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx);
static void mapc2p_vel_elc(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx);
static void mapc2p_vel_ion(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx);
static void bmag_func(double t, const double *xn, double *bmag, void *ctx);
static void bc_shift_func_lo(double t, const double *xn, double *fout, void *ctx);
static void bc_shift_func_up(double t, const double *xn, double *fout, void *ctx);
static void write_data(struct gkyl_tm_trigger *iot_conf, struct gkyl_tm_trigger *iot_phase,
                       gkyl_gyrokinetic_app *app, double t_curr, bool force_write);
static void calc_integrated_diagnostics(struct gkyl_tm_trigger *iot, gkyl_gyrokinetic_app *app, double t_curr, bool force_write);

struct gk_app_ctx create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.
  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Geometry and magnetic field.
  double a_shift   = 0.5;               // Parameter in Shafranov shift.
  double Z_axis    = 0.1414361745;       // Magnetic axis height [m].
  double R_axis    = 0.8867856264;       // Magnetic axis major radius [m].
  double B_axis    = 1.4;                // Magnetic field at the magnetic axis [T].
  double R_LCFSmid = 1.0870056099999; // Major radius of the LCFS at the outboard midplane [m
  double x_inner   = 0.04;               // Radial extent inside LCFS    
  double x_outer   = 0.08;               // Radial extent outside LCFS
  double Rmid_min  = R_LCFSmid - x_inner;      // Minimum midplane major radius of simulation box [m].
  double Rmid_max  = R_LCFSmid + x_outer;      // Maximum midplane major radius of simulation box [m].
  double R0        = 0.5*(Rmid_min+Rmid_max);  // Major radius of the simulation box [m].
  double a_mid     = R_LCFSmid-R_axis;   // Minor radius at outboard midplane [m].
  // Redefine a_mid with Shafranov shift, to ensure LCFS radial location.
  a_mid = R_axis/a_shift - sqrt(R_axis*(R_axis - 2*a_shift*R_LCFSmid + 2*a_shift*R_axis))/a_shift;
  double r0        = R0-R_axis;           // Minor radius of the simulation box [m].
  double B0        = B_axis*(R_axis/R0);  // Magnetic field magnitude in the simulation box [T].
  double kappa     = 1.4;                // Elongation (=1 for no elongation).
  double delta     = -0.38;                // Triangularity (=0 for no triangularity).

  // Plasma parameters. Chosen based on the value of a cubic sline
  // between the last TS data inside the LCFS and the probe data in
  // in the far SOL, near R=0.475 m.
  double AMU = 2.01410177811;
  double mi  = mp*AMU;   // Deuterium ions.
  double Te0 = 100*eV;
  double Ti0 = 100*eV;
  double n0  = 2.0e19;   // [1/m^3]
  double Bref = 1.129;   // Reference magnetic field [T].
  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi); // Thermal speeds.
  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Configuration domain parameters 
  double Lx        = Rmid_max-Rmid_min;   // Domain size along x.
  double Lz        = 2.*M_PI-1e-10;       // Domain size along magnetic field.
  double Ly        = 150*rho_s;                 // from rt_gk_d3d_iwl_3x2v_p1
  double x_min     = 0.;
  double x_max     = Lx;
  double y_min     = -Ly/2.;
  double y_max     =  Ly/2.;
  double z_min     = -Lz/2.;
  double z_max     =  Lz/2.;
  double x_LCFS    = R_LCFSmid - Rmid_min; // Radial location of the last closed flux surface.

  // Magnetic safety factor in the center of domain.
  double q0        = qprofile(r_x(0.5*(x_min+x_max),a_mid,x_inner),R_axis);   

  // Collision frequencies
  double nuFrac = 0.1;
  // Electron-electron collision freq.
  double logLambdaElc = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuElc = nuFrac * logLambdaElc * pow(eV, 4) * n0 /
    (6*sqrt(2.) * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(me) * pow(Te0,3./2.));
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4) * n0 /
    (12 * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(mi) * pow(Ti0,3./2.));

  // Source parameters
  double num_sources = 1;
  bool adaptive_source = true;
  char adapt_mom_type [26];
  strcpy(adapt_mom_type, "M2");
  double adapt_mom_rate_target = 0.5e6; // Input power of the source in W (we put 600kw because the recy source is ~100)
  // OMP source parameters
  double n_srcOMP = 0.75*2.4e23;
  double x_srcOMP = x_min;
  double Te_srcOMP = 3*Te0;
  double Ti_srcOMP = 3*Ti0;
  double sigma_srcOMP = 0.03*Lx;
  double floor_src = 1e-2;
  double power_src = 0.5e6;
  // Recycling source parameters
  double n_srcRECY = 0.99*n_srcOMP;
  double x_srcRECY = 0.5*x_LCFS;
  double Te_srcRECY = 20*eV;
  double Ti_srcRECY = 20*eV;
  double sigmax_srcRECY = 0.25*x_LCFS;
  double sigmaz_srcRECY = 0.25*Lz;
  // Grid parameters
  int num_cell_x = 16;
  int num_cell_y = 10;
  int num_cell_z = 12;
  int num_cell_vpar = 8;
  int num_cell_mu = 4;
  int poly_order = 1;
  // Velocity box dimensions
  double vpar_max_elc = 5.*vte;
  double mu_max_elc   = 1*me*pow(4*vte,2)/(2*B0);
  double vpar_max_ion = 5.*vti;
  double mu_max_ion   = 1*mi*pow(4*vti,2)/(2*B0);
  double final_time = 5.e-8;
  int num_frames = 1;
  double write_phase_freq = 1;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-3; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .a_shift = a_shift,
    .R_axis = R_axis,
    .R0     = R0    ,
    .a_mid  = a_mid ,
    .x_inner = x_inner,
    .r0     = r0    ,
    .B0     = B0    ,
    .kappa  = kappa ,
    .delta  = delta ,
    .q0     = q0    ,
    .Lx     = Lx    ,
    .Ly     = Ly    ,
    .Lz     = Lz    ,
    .x_min = x_min,  .x_max = x_max,
    .y_min = y_min,  .y_max = y_max,
    .z_min = z_min,  .z_max = z_max,
    .Bref = Bref,
    .x_LCFS = x_LCFS,
    .me = me,  .qe = qe,
    .mi = mi,  .qi = qi,
    .n0 = n0,  .Te0 = Te0,  .Ti0 = Ti0,
    .nuFrac = nuFrac,  .nuElc = nuElc,  .nuIon = nuIon,
    .num_sources = num_sources,
    .adaptive_source = adaptive_source,
    .adapt_mom_rate_target = adapt_mom_rate_target,
    .n_srcOMP     = n_srcOMP    ,
    .x_srcOMP     = x_srcOMP    ,
    .Te_srcOMP    = Te_srcOMP   ,
    .Ti_srcOMP    = Ti_srcOMP   ,
    .sigma_srcOMP = sigma_srcOMP,
    .n_srcRECY     = n_srcRECY    ,
    .x_srcRECY     = x_srcRECY    ,
    .Te_srcRECY    = Te_srcRECY   ,
    .Ti_srcRECY    = Ti_srcRECY   ,
    .sigmax_srcRECY = sigmax_srcRECY,
    .sigmaz_srcRECY = sigmaz_srcRECY,
    .floor_src    = floor_src   ,
    .num_cell_x     = num_cell_x,
    .num_cell_y     = num_cell_y,
    .num_cell_z     = num_cell_z,
    .num_cell_vpar  = num_cell_vpar,
    .num_cell_mu    = num_cell_mu,
    .cells = {num_cell_x, num_cell_y, num_cell_z, num_cell_vpar, num_cell_mu},
    .poly_order   = poly_order,
    .vpar_max_elc = vpar_max_elc,  .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,  .mu_max_ion = mu_max_ion,
    .write_phase_freq = write_phase_freq,
    .final_time = final_time,  .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  strcpy(ctx.adapt_mom_type, adapt_mom_type);
  return ctx;
}

int 
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
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
  int my_rank = 0;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    gkyl_comm_get_rank(comm, &my_rank);
#endif
  
  // electrons
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe, .mass = ctx.me,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

/*
    .flr = {
      .type = GKYL_GK_FLR_PADE_CONST,
      .Tperp = 2*ctx.Te0 //200eV
    },
*/

/*
    .init_from_file = {
       .type = GKYL_IC_IMPORT_F,
       .file_name = "restart_24x12_1.2ms-elc.gkyl"
    },
*/  
    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

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
      .num_sources = ctx.num_sources,
      .evolve = true,
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
        .density = density_srcRECY,
        .upar = zero_func,
        .temp = temp_elc_srcRECY,
      },
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = {"HamiltonianMoments"},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {"HamiltonianMoments"},
      }
    },
    .bcx = {
      .lower={.type = GKYL_SPECIES_ABSORB,},
      .upper={.type = GKYL_SPECIES_ABSORB,},
    },
    .bcz = {
      .lower={.type = GKYL_SPECIES_GK_IWL,
              .aux_parameter = ctx.x_LCFS,
              .aux_profile = bc_shift_func_lo, // taken from rt d3d 
              .aux_ctx = &ctx,
      },
      .upper={.type = GKYL_SPECIES_GK_IWL,
              .aux_parameter = ctx.x_LCFS,
              .aux_profile = bc_shift_func_up, // taken from rt d3d 
              .aux_ctx = &ctx,
      },
    },
    .num_diag_moments = 9,
    .diag_moments = {"HamiltonianMoments", "BiMaxwellianMoments", "M0", "M1", "M2par", "M2perp", "M2", "M3par", "M3perp"},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { "HamiltonianMoments" },
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { "HamiltonianMoments" },
    },
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,
/*
    .flr = {
      .type = GKYL_GK_FLR_PADE_CONST,
      .Tperp = 2*ctx.Ti0 //200eV
    },
*/
/*
    .init_from_file = {
       .type = GKYL_IC_IMPORT_F,
       .file_name = "restart_24x12_1.2ms-ion.gkyl"
    },
*/
    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

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
      .num_sources = ctx.num_sources,
      .evolve = true,
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
        .density = density_srcRECY,
        .upar = zero_func,
        .temp = temp_ion_srcRECY,
      },
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = {"HamiltonianMoments"},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {"HamiltonianMoments"},
      }
    },
    .bcx = {
      .lower={.type = GKYL_SPECIES_ABSORB,},
      .upper={.type = GKYL_SPECIES_ABSORB,},
    },
    .bcz = {
      .lower={.type = GKYL_SPECIES_GK_IWL,
              .aux_parameter = ctx.x_LCFS,
              .aux_profile = bc_shift_func_lo, // taken from rt d3d 
              .aux_ctx = &ctx,
      },
      .upper={.type = GKYL_SPECIES_GK_IWL,
              .aux_parameter = ctx.x_LCFS,
              .aux_profile = bc_shift_func_up, // taken from rt d3d 
              .aux_ctx = &ctx,
      },
    },
    .num_diag_moments = 9,
    .diag_moments = {"HamiltonianMoments", "BiMaxwellianMoments", "M0", "M1", "M2par", "M2perp", "M2", "M3par", "M3perp"},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { "HamiltonianMoments" },
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { "HamiltonianMoments" },
    },
  };

  struct gkyl_poisson_bias_plane target_corner_bc = {
    .dir = 0, // Direction perpendicular to the plane.
    .loc = ctx.x_LCFS, // Location of the plane in the 'dir' dimension.
    .val = 0.0, // Biasing value.
  };
  
  struct gkyl_poisson_bias_plane_list bias_plane_list = {
    .num_bias_plane = 1,
    .bp = &target_corner_bc,
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES_IWL,
    // .fem_parbc = GKYL_FEM_PARPROJ_DIRICHLET,
    .polarization_bmag = ctx.Bref,
    .xLCFS = ctx.x_LCFS,
    .poisson_bcs = {.lo_type = {GKYL_POISSON_DIRICHLET},
                    .up_type = {GKYL_POISSON_DIRICHLET},
                    .lo_value = {0.0}, .up_value = {0.0},
                   },
    .bias_plane_list = &bias_plane_list,
    .time_rate_diagnostics = true,
  };

  // Geometry
  struct gkyl_gyrokinetic_geometry geometry = {
    .geometry_id = GKYL_MAPC2P,
    .world = {0.},
    .mapc2p = mapc2p, // mapping of computational to physical space
    .c2p_ctx = &ctx,
    .bmag_func = bmag_func, // magnetic field magnitude
    .bmag_ctx = &ctx    
  };

  // Parallelism
  struct gkyl_app_parallelism_inp parallelism = {
    .comm = comm,
    .cuts = {app_args.cuts[0], app_args.cuts[1], app_args.cuts[2]},
    .use_gpu = app_args.use_gpu,
  };

  // Adaptive source parameters
  struct gkyl_gyrokinetic_source_adaptive adapt_src_params = {
    .mom_rate_target = ctx.adapt_mom_rate_target,
    .num_boundaries = 1,
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };
  strcpy(adapt_src_params.mom_type, ctx.adapt_mom_type);

  // GK app
  struct gkyl_gk *gk = gkyl_malloc(sizeof *gk);
  memset(gk, 0, sizeof(*gk));

  strcpy(gk->name, "rt_gk_tcv_iwl_3x2v_p1");

  gk->cfl_frac_omegaH = 1.0e9;
  gk->cfl_frac = 1.0;
  gk->cdim = 3;
  gk->vdim = 2;
  gk->lower[0] = ctx.x_min;
  gk->lower[1] = ctx.y_min;
  gk->lower[2] = ctx.z_min;
  gk->upper[0] = ctx.x_max;
  gk->upper[1] = ctx.y_max;
  gk->upper[2] = ctx.z_max;
  gk->cells[0] = cells_x[0];
  gk->cells[1] = cells_x[1];
  gk->cells[2] = cells_x[2];
  gk->poly_order = ctx.poly_order;
  gk->basis_type = app_args.basis_type;
  gk->geometry = geometry;
  gk->num_periodic_dir = 1;
  gk->periodic_dirs[0] = 1;
  gk->num_species = 2;
  gk->species[0] = elc;
  gk->species[1] = ion;
  gk->field = field;
  gk->parallelism = parallelism;
  gk->adaptive_source = ctx.adaptive_source;
  gk->adaptive_src_params = adapt_src_params;

  // create app object
  if (my_rank == 0) printf("Creating app object ...\n");
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(gk);
  
  // Initial and final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.final_time;
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
  struct gkyl_tm_trigger trig_write_conf = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_write_phase = { .dt = t_end/(ctx.write_phase_freq*num_frames), .tcurr = t_curr, .curr = frame_curr};
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  if (!app_args.is_restart) {
    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
    write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false);
  }

  // Initial time-step.
  double dt = t_end-t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  // MAIN TIME LOOP
  if (my_rank == 0) printf("Starting main loop ...\n");
  long step = 1, num_steps = app_args.num_steps;
  while ((t_curr < t_end) && (step <= num_steps)) {
    // if (step == 1 || step % 100 == 0)
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g mus ...", t_curr*1.0e6);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g mus\n", status.dt_actual*1.0e6);

    if (!status.success) {
      if (my_rank == 0) gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, t_curr > t_end);
    write_data(&trig_write_conf, &trig_write_phase, app, t_curr, t_curr > t_end);

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
        write_data(&trig_write_conf, &trig_write_phase, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }
  if (my_rank == 0) printf(" ... finished\n");
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

// End of main function
// ---------------------------------------------------------------

// Function definitions
// ---------------------------------------------------------------

// Geometry related functions 
double r_x(double x, double a_mid, double x_inner)
{
  return x+a_mid-x_inner;
}

double qprofile(double r, double R_axis) 
{
  // Magnetic safety factor as a function of minor radius r.
  double a[4] = {484.0615913225881, -1378.25993228584, 1309.3099150729233, 
    -414.13270311478726};
  return a[0]*pow(r+R_axis,3.0) + a[1]*pow(r+R_axis,2.0) + a[2]*(r+R_axis) + a[3];
}

double R_rtheta(double r, double theta, void *ctx)
{
  // Major radius as a function of minor radius r and poloidal angle theta.
  struct gk_app_ctx *app = ctx;
  double a_shift = app->a_shift;
  double R_axis = app->R_axis;
  double delta = app->delta;
  return R_axis - a_shift*r*r/(2.*R_axis) + r*cos(theta + asin(delta)*sin(theta));
}

double Z_rtheta(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Z_axis = app->Z_axis;
  double kappa = app->kappa;
  return Z_axis + kappa*r*sin(theta);
}

double dRdr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double a_shift = app->a_shift;
  double R_axis = app->R_axis;
  double delta = app->delta;
  return - a_shift*r/(R_axis) + cos(theta + asin(delta)*sin(theta));
}

double dRdtheta(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double delta = app->delta;
  return -r*sin(theta + asin(delta)*sin(theta))*(1.+asin(delta)*cos(theta));
}

double dZdr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double kappa = app->kappa;
  return kappa*sin(theta);
}

double dZdtheta(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double kappa = app->kappa;
  return kappa*r*cos(theta);
}

double Jr(double r, double theta, void *ctx)
{
  return R_rtheta(r,theta,ctx)*( dRdr(r,theta,ctx) * dZdtheta(r,theta,ctx)
		                -dRdtheta(r,theta,ctx) * dZdr(r,theta,ctx) );
}

struct integrand_ctx {
  struct gk_app_ctx *app_ctx;
  double r;
};

double integrand(double t, void *int_ctx)
{
  struct integrand_ctx *inctx = int_ctx;
  double r = inctx->r;
  struct gk_app_ctx *app = inctx->app_ctx;
  return Jr(r,t,app) / pow(R_rtheta(r,t,app),2);
}

double dPsidr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  integral = gkyl_dbl_exp(integrand, &tmp_ctx, 0., 2.*M_PI, 7, 1e-10);

  double B0 = app->B0;
  double a_mid = app->a_mid;
  double R_axis = app->R_axis;
  return ( B0*R_axis/(2.*M_PI*qprofile(r,R_axis)))*integral.res;
}

double alpha(double r, double theta, double phi, void *ctx)
{
  double twrap = theta;
  while (twrap < -M_PI) twrap = twrap+2.*M_PI;
  while (M_PI < twrap) twrap = twrap-2.*M_PI;

  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  if (0. < twrap) {
    integral = gkyl_dbl_exp(integrand, &tmp_ctx, 0., twrap, 7, 1e-10);
  } else {
    integral = gkyl_dbl_exp(integrand, &tmp_ctx, twrap, 0., 7, 1e-10);
    integral.res = -integral.res;
  }

  double B0 = app->B0;
  double R_axis = app->R_axis;

  return phi - B0*R_axis*integral.res/dPsidr(r,theta,ctx);
}

double Bphi(double R, void *ctx)
{
  // Toroidal magnetic field.
  struct gk_app_ctx *app = ctx;
  double B0 = app->B0;
  double R0 = app->R0;
  return B0*R0/R;
}

double gradr(double r, double theta, void *ctx)
{
  return (R_rtheta(r,theta,ctx)/Jr(r,theta,ctx))*sqrt(pow(dRdtheta(r,theta,ctx),2) + pow(dZdtheta(r,theta,ctx),2));
}

void zero_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

// Source functions
void density_srcOMP(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double n_srcOMP = app->n_srcOMP;
  double x_srcOMP = app->x_srcOMP;
  double sigma_srcOMP = app->sigma_srcOMP;
  double floor_src = app->floor_src;

  fout[0] = n_srcOMP*(exp(-(pow(x-x_srcOMP,2))/(2.*pow(sigma_srcOMP,2)))+floor_src);
}
void temp_elc_srcOMP(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

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
void temp_ion_srcOMP(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

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
// Recycling source
void density_srcRECY(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double n_src = app->n_srcRECY;
  double x_src = app->x_srcRECY;
  double sigmax_src = app->sigmax_srcRECY;
  double sigmaz_src = app->sigmaz_srcRECY;
  double floor_src = app->floor_src;
  
  double z_envelope = exp(-(pow(z + M_PI,2)) / (2.0 * pow(sigmaz_src,2))) + exp(-(pow(z - M_PI,2)) / (2.0 * pow(sigmaz_src,2)));
  fout[0] = z_envelope * n_src * (exp(-(pow(x - x_src,2)) / (2.0 * pow(sigmax_src,2))) + floor_src);
}
void temp_elc_srcRECY(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double Te_src = app->Te_srcRECY;

  fout[0] = Te_src*3./8.;
}
void temp_ion_srcRECY(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double Ti_src = app->Ti_srcRECY;

  fout[0] = Ti_src;
}

// Density initial condition
void density_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];
  struct gk_app_ctx *app = ctx;
  double n0 = app->n0;
  fout[0] = n0*(0.5*(1.+tanh(3.*(.1-10.*x)))+0.01);
}

// Electron temperature initial conditions
void temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];
  struct gk_app_ctx *app = ctx;
  double Te0 = app->Te0;
  fout[0] = 6.*Te0*(0.5*(1.+tanh(3.*(-.1-10.*x)))+0.01);
}

// Ion temperature initial conditions
void temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];
  struct gk_app_ctx *app = ctx;
  double Ti0 = app->Ti0;
  fout[0] = 6.*Ti0*(0.5*(1.+tanh(3.*(-.1-10.*x)))+0.01);
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

// Geometry evaluation functions for the gk app
void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];
  struct gk_app_ctx *app = ctx;
  double r0 = app->r0;
  double q0 = app->q0;
  double a_mid = app->a_mid;
  double x_inner = app->x_inner;
  double r = r_x(x,a_mid,x_inner);
  // Map to cylindrical (R, Z, phi) coordinates.
  double R   = R_rtheta(r, z, ctx);
  double Z   = Z_rtheta(r, z, ctx);
  double phi = -q0/r0*y - alpha(r, z, 0, ctx);
  // Map to Cartesian (X, Y, Z) coordinates.
  double X = R*cos(phi);
  double Y = R*sin(phi);
  xp[0] = X; xp[1] = Y; xp[2] = Z;
}

// Taken from rt gk d3d 3x2c, is this the non uniform v grid mapping?
void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;
  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_elc*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_elc*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_elc*2.0*pow(cvpar,2);
  // Quadratic map in mu.
  vp[1] = mu_max_elc*pow(cmu,2);
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;
  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_ion*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_ion*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_ion*2.0*pow(cvpar,2);
  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}

void bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];
  struct gk_app_ctx *app = ctx;
  double a_mid = app->a_mid;
  double x_inner = app->x_inner;
  double r = r_x(x,a_mid,x_inner);
  double Bt = Bphi(R_rtheta(r,z,ctx),ctx);
  double Bp = dPsidr(r,z,ctx)/R_rtheta(r,z,ctx)*gradr(r,z,ctx);
  fout[0] = sqrt(Bt*Bt + Bp*Bp);
}

// Taken from rt gk d3d 3x2v
void bc_shift_func_lo(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0];
  struct gk_app_ctx *app = ctx;
  double r0 = app->r0;
  double q0 = app->q0;
  double a_mid = app->a_mid;
  double Lz = app->Lz;
  double x_inner = app->x_inner;
  double r = r_x(x,a_mid,x_inner);
  double R_axis = app->R_axis;
  fout[0] = -r0/q0*qprofile(r,R_axis)*Lz;
}

void bc_shift_func_up(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  bc_shift_func_lo(t, xc, fout, ctx);
  fout[0] *= -1;
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
write_data(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

    gkyl_gyrokinetic_app_write_conf(app, t_curr, frame);

    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_write_field_energy(app);

    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

    gkyl_gyrokinetic_app_write_phase(app, t_curr, frame);
  }
}
