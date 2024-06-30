#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

struct gk_step_ctx {
  int cdim, vdim; // Dimensionality.
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
  double nuFrac; // Factor to multiply collision frequencies
  double B0; // reference magnetic field
  double n0; // reference density
  double nsource;
  // Source parameters
  double T_source; // Source electron temperature
  double cx;
  double cz;
  // Simulation parameters
  int Nx; // Cell count (configuration space: x-direction).
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double lower_x, upper_x; // Limits of the domain along x.
  double Lx; // Domain size (configuration space: x-direction).
  double Lz; // Domain size (configuration space: z-direction).
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double t_end; // end time
  int num_frames; // number of output frames
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
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
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double n0 = app->n0;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = n0*exp(-(x-xcenter)*(x-xcenter)/2/cx/cx) * exp(-z*z/2/cz/cz);
  if (n/n0 < 1e-5)
    n = n0*1e-5;
  fout[0] = n;
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
eval_density_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double nsource = app->nsource;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = nsource*exp(-(x-xcenter)*(x-xcenter)/2/cx/cx) * exp(-z*z/2/cz/cz);
  if (n/nsource < 1e-5)
    n = nsource*1e-5;
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

double plasma_frequency(double n, double m)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  return sqrt(n*eV*eV/m/eps0);
}

struct gk_step_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double Te = 100*2.8*eV;
  double Ti = 150*2.8*eV;
  double B0 = 2.51; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19/2.8; // Particle density in 1/m^3
                             
  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);

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
  double lower_x = 0.934;
  double upper_x = 1.4688;
  double Lx = upper_x - lower_x;
  double Lz = 3.14*2;

  double vpar_max_elc = 4.0*vtElc;
  double mu_max_elc = 18*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 4.0*vtIon;
  double mu_max_ion = 18*mi*vtIon*vtIon/(2.0*B0);

  // Number of cells.
  int Nx = 4;
  int Nz = 16;
  int Nvpar = 16;
  int Nmu = 8;

  double t_end = 5.0e-7; 
  double num_frames = 1;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_step_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .Te = Te, 
    .Ti = Ti, 
    .vtIon = vtIon,
    .vtElc = vtElc,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuFrac = nuFrac,
    .B0 = B0, 
    .n0 = n0, 
    .T_source = T_source, 
    .nsource = nsource,
    .cx = cx,
    .cz = cz,
    .lower_x = lower_x,
    .upper_x = upper_x,
    .Lx = Lx,
    .Lz = Lz, 
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    .Nx = Nx,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .t_end = t_end, 
    .num_frames = num_frames, 
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_multib_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    //gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_multib_app_calc_integrated_mom(app, t_curr);
  }
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_multib_app* app, double t_curr, bool force_write)
{
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_gyrokinetic_multib_app_write(app, t_curr, frame);

    gkyl_gyrokinetic_multib_app_calc_mom(app);
    gkyl_gyrokinetic_multib_app_write_mom(app, t_curr, frame);
    //gkyl_gyrokinetic_multib_app_write_source_mom(app, t_curr, frame);

    //gkyl_gyrokinetic_multib_app_calc_field_energy(app, t_curr);
    //gkyl_gyrokinetic_multib_app_write_field_energy(app);

    gkyl_gyrokinetic_multib_app_calc_integrated_mom(app, t_curr);
    gkyl_gyrokinetic_multib_app_write_integrated_mom(app);
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

  struct gk_step_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Electrons.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_elc,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 1,
      .collide_with = { "ion"},
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_source,
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

    .bcx = {
      .lower={.type = GKYL_SPECIES_ZERO_FLUX,},
      .upper={.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // Ions.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
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
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm
      .self_nu = evalNuIon,
      .num_cross_collisions = 1,
      .collide_with = { "elc"},
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_source,
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

    .bcx = {
      .lower={.type = GKYL_SPECIES_ZERO_FLUX,},
      .upper={.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .fem_parbc = GKYL_FEM_PARPROJ_NONE, 
    .poisson_bcs = {.lo_type = {GKYL_POISSON_DIRICHLET}, 
                    .up_type = {GKYL_POISSON_DIRICHLET}, 
                    .lo_value = {0.0}, .up_value = {0.0}}, 
  };

  /* Block layout

     +------+
     |2(bup)|
     |      |
     +------+
     |1(bmid)|
     |      |
     +------+
     |0(blo)|
     |      |
     +------+
    
  */  




  struct gkyl_gk blo = {
    .lower = { ctx.lower_x, -ctx.Lz/2.0 },
    .upper = { ctx.upper_x,  -ctx.Lz/4.0 },
    .cells = { cells_x[0], cells_x[1]/4 },


    .geometry = {
      .geometry_id = GKYL_GEOMETRY_FROMFILE,
      //.world = {0.0},
      //.geometry_id = GKYL_TOKAMAK,
      //.tok_efit_info = &inp,
      //.tok_grid_info = &ginp,
    },

    .block_connections = {
      .connections[0] = { // x-direction
        { .bid = 0, .dir = 0, .edge = GKYL_BLOCK_EDGE_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_BLOCK_EDGE_PHYSICAL },
      },
      .connections[1] = { // z-direction
        { .bid = 0, .dir = 1, .edge = GKYL_BLOCK_EDGE_PHYSICAL }, // physical boundary
        { .bid = 1, .dir = 1, .edge = GKYL_BLOCK_EDGE_LOWER_POSITIVE },
      }
    }, 

    .cuts = { 1, 1 },
  };

  struct gkyl_gk bmid = {
    .lower = { ctx.lower_x, -ctx.Lz/4.0 },
    .upper = { ctx.upper_x,  ctx.Lz/4.0 },
    .cells = { cells_x[0], cells_x[1]/2 },

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_FROMFILE,
      //.world = {0.0},
      //.geometry_id = GKYL_TOKAMAK,
      //.tok_efit_info = &inp,
      //.tok_grid_info = &ginp,
    },

    .block_connections = {
      .connections[0] = { // x-direction
        { .bid = 1, .dir = 0, .edge = GKYL_BLOCK_EDGE_PHYSICAL }, // physical boundary
        { .bid = 1, .dir = 0, .edge = GKYL_BLOCK_EDGE_PHYSICAL },
      },
      .connections[1] = { // z-direction
        { .bid = 0, .dir = 1, .edge = GKYL_BLOCK_EDGE_UPPER_POSITIVE}, 
        { .bid = 2, .dir = 1, .edge = GKYL_BLOCK_EDGE_LOWER_POSITIVE },
      }
    }, 

    .cuts = { 1, 1 },
  };

  struct gkyl_gk bup= {
    .lower = { ctx.lower_x, ctx.Lz/4.0 },
    .upper = { ctx.upper_x,  ctx.Lz/2.0 },
    .cells = { cells_x[0], cells_x[1]/4 },

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_FROMFILE,
      //.world = {0.0},
      //.geometry_id = GKYL_TOKAMAK,
      //.tok_efit_info = &inp,
      //.tok_grid_info = &ginp,
    },

    .block_connections = {
      .connections[0] = { // x-direction
        { .bid = 2, .dir = 0, .edge = GKYL_BLOCK_EDGE_PHYSICAL }, // physical boundary
        { .bid = 2, .dir = 0, .edge = GKYL_BLOCK_EDGE_PHYSICAL },
      },
      .connections[1] = { // z-direction
        { .bid = 1, .dir = 1, .edge = GKYL_BLOCK_EDGE_UPPER_POSITIVE}, // physical boundary
        { .bid = 2, .dir = 1, .edge = GKYL_BLOCK_EDGE_PHYSICAL},
      }
    }, 

    .cuts = { 1, 1 },
  };

  // GK app
  struct gkyl_gk_multib app_inp = {
    .name = "gk_mb_step_out_2x2v_p1",

    .cdim = 2, .vdim = 2,
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .num_blocks = 3,
    .blocks = { &blo, &bmid, &bup },

    .num_periodic_dir = 0,
    .periodic_dirs = {  },

    .num_species = 2,
    .species = { elc, ion},

    .field = field,

    .use_gpu = app_args.use_gpu,
    .use_mpi = app_args.use_mpi,



  };

  // Create app object.
  gkyl_gyrokinetic_multib_app *app = gkyl_gyrokinetic_multib_app_new(&app_inp);

  // Initial and final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_multib_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_multib_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_gyrokinetic_multib_app_apply_ic(app, t_curr);
  }  

  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
  write_data(&trig_write, app, t_curr, false);

  double dt = t_end-t_curr; // Initial time step.
  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1, num_steps = app_args.num_steps;
  while ((t_curr < t_end) && (step <= num_steps)) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Taking time-step at t = %g ...", t_curr);
    struct gkyl_update_status status = gkyl_gyrokinetic_multib_update(app, dt);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_multib_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
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

      gkyl_gyrokinetic_multib_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_multib_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_multib_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_multib_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_multib_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
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

  gkyl_gyrokinetic_multib_app_stat_write(app);
  
  // Fetch simulation statistics.
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_multib_app_stat(app);

  gkyl_gyrokinetic_multib_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);

  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_multib_app_release(app);


#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif  

  return 0;
}
