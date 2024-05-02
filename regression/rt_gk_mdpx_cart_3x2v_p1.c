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

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

struct gk_mdpx_ctx {
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massAr; // Argon mass
  double Te; // electron temperature
  double TAr; // Argon temperature
  double vtElc;
  double vtAr;
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuFrac; // Factor to multiply collision frequencies
  double B0; // reference magnetic field
  double n0; // reference density
  double r_prof; // Radial extent of density profile initial conditions.
  double L_prof; // Length of density profile initial conditions.  
  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Ny; // Cell count (configuration space: y-direction).
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  double Lx; // Domain size (configuration space: x-direction).
  double Ly; // Domain size (configuration space: y-direction).
  double Lz; // Domain size (configuration space: z-direction).
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_Ar; // Velocity space extents in vparallel for Ar
  double mu_max_Ar; // Velocity space extents in mu for Ar
  double t_end; // end time
  int num_frames; // number of output frames
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mdpx_ctx *app = ctx;
  double x = xn[0], y = xn[1], z = xn[2];

  double r_prof = app->r_prof;
  double L_prof = app->L_prof;
  double n0 = app->n0;

  double r = sqrt(x * x + y * y);
  pcg32_random_t rng = gkyl_pcg32_init(0);
  double perturb = 2.0e-3 * (1.0 - 0.5 * gkyl_pcg32_rand_double(&rng));
  // Set number density
  fout[0] = n0 * (1.0 + perturb) * (1.0 - tanh((r - r_prof) / L_prof));
}

void
eval_density_arneut(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0e21;
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
  struct gk_mdpx_ctx *app = ctx;
  double T = app->Te;
  fout[0] = T;
}

void
eval_temp_ar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mdpx_ctx *app = ctx;
  double T = app->TAr;
  fout[0] = T;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mdpx_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mdpx_ctx *app = ctx;
  fout[0] = app->nuIon;
}

static inline void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  // Set physical coordinates (X, Y, Z) from computational coordinates (x, y, z).
  xp[0] = zc[0]; xp[1] = zc[1]; xp[2] = zc[2];
}

void
bmag_func(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_mdpx_ctx *app = ctx;

  double B0 = app->B0;

  // Set magnetic field strength.
  fout[0] = B0;
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

struct gk_mdpx_ctx
create_ctx(void)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mAr = 39.95*GKYL_PROTON_MASS; // Ar ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double Te = 4.0*eV;
  double TAr = eV/40.0;
  double B0 = 1.0; // Magnetic field magnitude in Tesla
  double n0 = 1.0e14; // Particle density in 1/m^3
                             
  // Derived parameters.
  double vtElc = sqrt(Te/me);
  double vtAr = sqrt(TAr/mAr);

  double r_prof = 0.16; // Initial conditions radial extent.
  double L_prof = 0.04; // Initial conditions length.

  // Collision parameters.
  double nuFrac = 1.0;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(TAr/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mAr)*(TAr*sqrt(TAr)));

  // Simulation parameters.
  int Nx = 18; // Cell count (configuration space: x-direction).
  int Ny = 18; // Cell count (configuration space: y-direction).
  int Nz = 8; // Cell count (configuration space: z-direction).
  int Nvpar = 10; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 5; // Cell count (velocity space: magnetic moment direction).
  double Lx = 0.5; // Lx box size in meters
  double Ly = 0.5; // Ly box size in meters
  double Lz = 0.1; // Lz box size in meters

  double vpar_max_elc = 4.0*vtElc;
  double mu_max_elc = 18*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_Ar = 4.0*vtAr;
  double mu_max_Ar = 18.*mAr*vtAr*vtAr/(2.0*B0);

  double t_end = 1.0e-10; 
  double num_frames = 1;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_mdpx_ctx ctx = {
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massAr = mAr,
    .Te = Te, 
    .TAr = TAr, 
    .vtElc = vtElc,
    .vtAr = vtAr,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuFrac = nuFrac,
    .B0 = B0, 
    .n0 = n0, 
    .r_prof = r_prof,
    .L_prof = L_prof,
    .Nx = Nx,
    .Ny = Ny,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .Lx = Lx,  
    .Ly = Ly,  
    .Lz = Lz, 
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_Ar = vpar_max_Ar, 
    .mu_max_Ar = mu_max_Ar, 
    .t_end = t_end, 
    .num_frames = num_frames, 
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
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = force_write? iot->curr : iot->curr -1;

    gkyl_gyrokinetic_app_write(app, t_curr, frame);

    gkyl_gyrokinetic_app_calc_mom(app);
    gkyl_gyrokinetic_app_write_mom(app, t_curr, frame);
    gkyl_gyrokinetic_app_write_source_mom(app, t_curr, frame);

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
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_mdpx_ctx ctx = create_ctx(); // context for init functions

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], ctx.Nz);
  int NVPAR = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvpar);
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nmu);

  // Create global range.
  int ccells[] = { NX, NY, NZ };
  int cdim = sizeof(ccells) / sizeof(ccells[0]);
  struct gkyl_range cglobal_r;
  gkyl_create_global_range(cdim, ccells, &cglobal_r);

  // Create decomposition.
  int cuts[cdim];
#ifdef GKYL_HAVE_MPI  
  for (int d = 0; d < cdim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < cdim; d++) {
    cuts[d] = 1;
  }
#endif  
    
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(cdim, cuts, &cglobal_r);

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_gpu && app_args.use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
#else
    printf(" Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp,
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank, comm_size;
  gkyl_comm_get_rank(comm, &my_rank);
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < cdim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  for (int d = 0; d < cdim - 1; d++) {
    if (cuts[d] > 1) {
      if (my_rank == 0) {
        fprintf(stderr, "*** Parallelization only allowed in z. Number of ranks, %d, in direction %d cannot be > 1!\n", cuts[d], d);
      }
      goto mpifinalize;
    }
  }

  // electrons
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { NVPAR, NMU },
    .polarization_density = ctx.n0,
    .no_by = true, 

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
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massElc, ctx.massElc, ctx.chargeElc, ctx.chargeElc, ctx.Te, ctx.Te),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massElc, ctx.massAr, ctx.chargeElc, ctx.chargeIon, ctx.Te, ctx.TAr)
      },
      .bmag_mid = ctx.B0,
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 1,
      .collide_with = { "Ar1" },
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

    .bcx = {
      .lower={.type = GKYL_SPECIES_ZERO_FLUX,},
      .upper={.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_ZERO_FLUX,},
      .upper={.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcz = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // Ar1+ ions
  struct gkyl_gyrokinetic_species Ar1 = {
    .name = "Ar1",
    .charge = ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NVPAR, NMU },
    .polarization_density = ctx.n0,
    .no_by = true, 

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massAr, ctx.massAr, ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massAr, ctx.massElc, ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
      },
      .bmag_mid = ctx.B0,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
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

    .bcx = {
      .lower={.type = GKYL_SPECIES_ZERO_FLUX,},
      .upper={.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_ZERO_FLUX,},
      .upper={.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcz = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };


  // neutral Ar
  struct gkyl_gyrokinetic_neut_species Ar0 = {
    .name = "Ar0", .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, -ctx.vpar_max_Ar, -ctx.vpar_max_Ar},
    .upper = { ctx.vpar_max_Ar, ctx.vpar_max_Ar, ctx.vpar_max_Ar },
    .cells = { NVPAR, NVPAR, NVPAR},
    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_arneut,
      .ctx_upar = &ctx,
      .udrift= eval_udrift,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .bcx = { 
      .lower = { .type = GKYL_SPECIES_REFLECT },
      .upper = { .type = GKYL_SPECIES_REFLECT },
    },
    .bcy = { 
      .lower = { .type = GKYL_SPECIES_REFLECT },
      .upper = { .type = GKYL_SPECIES_REFLECT },
    },
    .bcz = { 
      .lower = { .type = GKYL_SPECIES_REFLECT },
      .upper = { .type = GKYL_SPECIES_REFLECT },
    },
    
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2"},
  };


  // field
  struct gkyl_gyrokinetic_field field = {
    .bmag_fac = ctx.B0, 
    .fem_parbc = GKYL_FEM_PARPROJ_NONE, 
    .poisson_bcs = {
      .lo_type = { GKYL_POISSON_DIRICHLET, GKYL_POISSON_DIRICHLET },
      .up_type = { GKYL_POISSON_DIRICHLET, GKYL_POISSON_DIRICHLET },
      .lo_value = { 0.0, 0.0 }, .up_value = { 0.0, 0.0}
    },
  };

  // GK app
  struct gkyl_gk gk = {
    .name = "gk_mdpx_cart_3x2v_p1",

    .cdim = 3, .vdim = 2,
    .lower = { -0.5 * ctx.Lx, -0.5 * ctx.Ly, -0.5 * ctx.Lz },
    .upper = { 0.5 * ctx.Lx, 0.5 * ctx.Ly, 0.5 * ctx.Lz },
    .cells = { NX, NY, NZ },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .mapc2p = mapc2p,
      .c2p_ctx = &ctx,
      .bmag_func = bmag_func,
      .bmag_ctx = &ctx
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {  },

    .num_species = 2,
    .species = { elc, Ar1 },

    .num_neut_species = 1,
    .neut_species = { Ar0 },

    .field = field,

    .use_gpu = app_args.use_gpu,

    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp->ranges[my_rank],
      .comm = comm
    }
  };

  // create app object
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);

  // start, end and initial time-step
  double t_curr = 0.0, t_end = ctx.t_end;
  double dt = t_end-t_curr;
  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc) };
  struct gkyl_tm_trigger trig_write = { .dt = t_end/num_frames };

  // initialize simulation
  gkyl_gyrokinetic_app_apply_ic(app, t_curr);
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
  write_data(&trig_write, app, t_curr, false);
  gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

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

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
    write_data(&trig_write, app, t_curr, false);

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

  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
  write_data(&trig_write, app, t_curr, false);
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
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_gyrokinetic_app_release(app);

  mpifinalize:
  ;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif  

  return 0;
}
