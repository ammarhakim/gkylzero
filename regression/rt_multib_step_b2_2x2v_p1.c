#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_null_comm.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

struct gkyl_block_geom*
create_block_geom(void)
{
  // Only do b2 in the block layout below.
  struct gkyl_block_geom *bgeom = gkyl_block_geom_new(2, 1);

  /* Block layout and coordinates

   x  
   ^  
   |
   4  +------------------+------------------+------------------+
   |  |b1                |b2                |b3                |
   |  |lower outer SOL   |middle outer sol  |upper outer sol   |
   |  |                  |                  |                  |
   3  +------------------+------------------+------------------+
   |  |b0               x|o b10            %|$ b4              |
   |  |lower outer PF   x|o outer core     %|$ upper outer PF  |
   |  |                 x|o                %|$                 |
   |  +------------------+------------------+------------------+
   2  +------------------+------------------+------------------+
   |  |b9               x|o b11            %|$ b5              |
   |  |lower inner PF   x|o inner core     %|$ upper inner PF  |
   |  |                 x|o                %|$                 |
   1  +------------------+------------------+------------------+
   |  |b8                |b7                |b6                |
   |  |lower inner SOL   |middle inner SOL  |upper inner SOL   |
   |  |                  |                  |                  |
   0  +------------------+------------------+------------------+

      0 -----------1------------2------------3 -> z

      Edges that touch coincide are physically connected unless
      otherwise indicated by a special symbol. Edges with a special
      symbol such as o,x,%, or % are instead connected to the other
      edge with the same symbol. Edges that do not coincide with
      another edge are a physical boundary.
  */  

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .filepath = "./data/eqdsk/step.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psisep = efit->psisep;
  gkyl_efit_release(efit);
  double psi_up_core = 1.8;
  double psi_up_pf = 1.8;
  double psi_lo_outer_sol = 0.934;
  double psi_lo_inner_sol = 1.45;

  int npsi_outer_sol = 4;
  int npsi_core = 4;
  int npsi_inner_sol = 4;
  int npsi_lower_pf = 4;
  int npsi_upper_pf = 4;

  double ntheta_lower  = 8;
  double ntheta_middle = 8;
  double ntheta_upper  = 8;

  double Lz = (M_PI-1e-14)*2.0;
  double theta_lo = -Lz/2.0, theta_up = Lz/2.0;

  // block 2. Middle outer SOL.
  gkyl_block_geom_set_block(bgeom, 0, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_middle},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_SOL_DN_OUT_MID,
          .rclose = 6.2,                // closest R to region of interest
          .rright = 6.2,                // Closest R to outboard SOL
          .rleft = 2.0,                 // closest R to inboard SOL
          .rmin = 1.1,                  // smallest R in machine
          .rmax = 6.2,                  // largest R in machine
          .use_cubics = false,          // Whether to use cubic representation of psi(R,Z) for field line tracing
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}
      }
    }
  );

  return bgeom;
}

struct gk_step_ctx {
  int cdim, vdim; // Dimensionality.
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double massAr; // Argon mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double TAr; // Argon temperature
  double vtIon;
  double vtElc;
  double vtAr;
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuFrac; // Factor to multiply collision frequencies
  double B0; // reference magnetic field
  double n0; // reference density
  double n0Ar; // Argon reference density
  double nsource;
  // Source parameters
  double T_source; // Source electron temperature
  double cx;
  double cz;
  // Space-time parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double vpar_max_Ar; // Velocity space extents in vparallel for Ar
  double mu_max_Ar; // Velocity space extents in mu for Ar
  double t_end; // end time
  int num_frames; // number of output frames
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct gk_step_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double mAr = 39.95*GKYL_PROTON_MASS; // Ar ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double Te = 100*2.8*eV;
  double Ti = 150*2.8*eV;
  double TAr = 40.0*eV;
  double B0 = 2.51; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19/2.8; // Particle density in 1/m^3
  double n0Ar = n0*0.0001/3.0; // Particle density in 1/m^3
                             
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

  double vpar_max_elc = 4.0*vtElc;
  double mu_max_elc = 18*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 4.0*vtIon;
  double mu_max_ion = 18*mi*vtIon*vtIon/(2.0*B0);

  double vpar_max_Ar = 4.0*vtAr;
  double mu_max_Ar = 18.*mAr*vtAr*vtAr/(2.0*B0);

  // Number of cells.
  int Nx = 4;
  int Nz = 8;
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
    .massAr = mAr,
    .Te = Te, 
    .Ti = Ti, 
    .TAr = TAr, 
    .vtIon = vtIon,
    .vtElc = vtElc,
    .vtAr = vtAr,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuFrac = nuFrac,
    .B0 = B0, 
    .n0 = n0, 
    .n0Ar = n0Ar,
    .T_source = T_source, 
    .nsource = nsource,
    .cx = cx,
    .cz = cz,
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    .vpar_max_Ar = vpar_max_Ar, 
    .mu_max_Ar = mu_max_Ar, 
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
init_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
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
init_density_Ar0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0Ar;
  double cz = app->cz/1.4/2.0;
  double zcenter = 3.14;
  double n = 0.0;
  if (z>0)
    n = n0 * exp(-(z-zcenter)*(z-zcenter)/(2.0*cz*cz));
  else
    n = n0 * exp(-(z+zcenter)*(z+zcenter)/(2.0*cz*cz));
  if (n < 1.0e8)
    n = 1.0e8;

  fout[0] = n;
}

void
init_density_Ar1(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0e5;
}

void
init_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
init_udrift_Ar0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; 
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
init_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->Te;
  fout[0] = T;
}

void
init_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->Ti;
  fout[0] = T;
}

void
init_temp_Ar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->TAr;
  fout[0] = T;
}

void
init_density_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
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
init_upar_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
init_temp_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double T = app->T_source;
  fout[0] = T;
}

void
init_nu_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  fout[0] = input->nuElc;
}

void
init_nu_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  fout[0] = input->nuIon;
}


void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_multib_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
//    gkyl_gyrokinetic_multib_app_calc_field_energy(app, t_curr);
//    gkyl_gyrokinetic_multib_app_calc_integrated_mom(app, t_curr);
  }
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_multib_app* app, double t_curr, bool force_write)
{
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_gyrokinetic_multib_app_write(app, t_curr, frame);

//    gkyl_gyrokinetic_multib_app_calc_field_energy(app, t_curr);
//    gkyl_gyrokinetic_multib_app_write_field_energy(app);
//
//    gkyl_gyrokinetic_multib_app_calc_integrated_mom(app, t_curr);
//    gkyl_gyrokinetic_multib_app_write_integrated_mom(app);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  
  // Construct block geometry.
  struct gkyl_block_geom *bgeom = create_block_geom();
  int nblocks = gkyl_block_geom_num_blocks(bgeom);

  struct gk_step_ctx ctx = create_ctx(); // Context for init functions.
  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Elc Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb elc_blocks[1];
  elc_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = init_density,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = init_density_source,
        .ctx_upar = &ctx,
        .upar = init_upar_source,
        .ctx_temp = &ctx,
        .temp = init_temp_source,      
      }, 
    },

  };


  struct gkyl_block_physical_bcs elc_phys_bcs[] = {
    // block 2 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    { .bidx = 0, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
  };

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { cells_v[0], cells_v[1] },
    .no_by = true,
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = init_nu_elc,
      .num_cross_collisions = 2,
      .collide_with = { "ion", "Ar1" },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 

    .radiation = {
      .radiation_id = GKYL_GK_RADIATION, 
      .num_cross_collisions = 1, 
      .collide_with = { "Ar1" },
      .z = 18,
      .charge_state = 1,
      .num_of_densities = 1, // Must be 1 for now
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

    .duplicate_across_blocks = true,
    .blocks = elc_blocks,
    .num_physical_bcs = 4,
    .bcs = elc_phys_bcs,
  };


  // Ion Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb ion_blocks[1];
  ion_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = init_density,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = init_density_source,
        .ctx_upar = &ctx,
        .upar = init_upar_source,
        .ctx_temp = &ctx,
        .temp = init_temp_source,      
      }, 
    },

  };

  struct gkyl_block_physical_bcs ion_phys_bcs[] = {
    // block 2 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    { .bidx = 0, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
  };

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },
    .no_by = true,
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = init_nu_ion,
      .num_cross_collisions = 2,
      .collide_with = { "elc", "Ar1" },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 
  
    .duplicate_across_blocks = true,
    .blocks = ion_blocks,
    .num_physical_bcs = 4,
    .bcs = ion_phys_bcs,
  };

  // Ar+1 Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb Ar1_blocks[1];
  Ar1_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .polarization_density = ctx.n0Ar,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = init_density_Ar1,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_Ar,
    },

  };

  struct gkyl_block_physical_bcs Ar1_phys_bcs[] = {
    // block 2 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    { .bidx = 0, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
  };

  struct gkyl_gyrokinetic_multib_species Ar1 = {
    .name = "Ar1",
    .charge = ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = {  ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { cells_v[0], cells_v[1] },
    .no_by = true,
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0Ar, // Density used to calculate coulomb logarithm
      .T_ref = ctx.TAr, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = init_nu_ion,
      .num_cross_collisions = 2,
      .collide_with = { "elc", "ion" },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
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

  
    .duplicate_across_blocks = true,
    .blocks = Ar1_blocks,
    .num_physical_bcs = 4,
    .bcs = Ar1_phys_bcs,
  };

  // Neutral Ar0 Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_neut_species_pb Ar0_blocks[1];
  Ar0_blocks[0] = (struct gkyl_gyrokinetic_multib_neut_species_pb) {

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = init_density_Ar0,
      .ctx_udrift = &ctx,
      .udrift = init_udrift_Ar0,
      .ctx_temp = &ctx,
      .temp = init_temp_Ar,
    },

  };

  struct gkyl_block_physical_bcs Ar0_phys_bcs[] = {
    // block 2 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
  };

  struct gkyl_gyrokinetic_multib_neut_species Ar0 = {
    .name = "Ar0",
    .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, -ctx.vpar_max_Ar, -ctx.vpar_max_Ar},
    .upper = {  ctx.vpar_max_Ar,  ctx.vpar_max_Ar,  ctx.vpar_max_Ar },
    .cells = { cells_v[0], cells_v[0], cells_v[0] },
    .is_static = true,
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2"},

    .duplicate_across_blocks = true,
    .blocks = Ar0_blocks,
    .num_physical_bcs = 4,
    .bcs = Ar0_phys_bcs,
  };

  // Field object
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    // No block specific field info for this simulation
  };

  struct gkyl_block_physical_bcs field_phys_bcs[] = {
    // block 2 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_POISSON_DIRICHLET},
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_POISSON_DIRICHLET},
  };

  struct gkyl_gyrokinetic_multib_field field = {
    .duplicate_across_blocks = true,
    .blocks = field_blocks, 
    .num_physical_bcs = 2,
    .bcs = field_phys_bcs,
  };

  struct gkyl_gyrokinetic_multib app_inp = {
    .name = "multib_step_b2_2x2v_p1",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .use_gpu = app_args.use_gpu,

    .block_geom = bgeom,
    .cfl_frac = 0.9,
    
    .enforce_positivity = false,

    .num_species = 3,
    .species = { elc, ion, Ar1},

    .num_neut_species = 1,
    .neut_species = { Ar0 },

    .field = field,

    .comm = comm
  };

  // Create app object.
  struct gkyl_gyrokinetic_multib_app *app = gkyl_gyrokinetic_multib_app_new(&app_inp);

  // Initial and final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
//  if (app_args.is_restart) {
//    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);
//
//    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
//      gkyl_gyrokinetic_multib_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
//        gkyl_array_rio_status_msg(status.io_status));
//      goto freeresources;
//    }
//
//    frame_curr = status.frame;
//    t_curr = status.stime;
//
//    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
//    gkyl_gyrokinetic_multib_app_cout(app, stdout, " at time = %g\n", t_curr);
//  }
//  else {
    gkyl_gyrokinetic_multib_app_apply_ic(app, t_curr);
//  }

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
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
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

//  gkyl_gyrokinetic_app_stat_write(app);
//
//  // Fetch simulation statistics.
//  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  freeresources:
  // Free resources after simulation completion.
  gkyl_block_geom_release(bgeom);
  gkyl_gyrokinetic_multib_app_release(app);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  
  return 0;
}
