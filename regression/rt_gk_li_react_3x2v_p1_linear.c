#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_gyrokinetic.h>

#include <rt_arg_parse.h>

struct gk_app_ctx {
  int cdim, vdim; // Dimensionality.
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double massLi; // Li mass
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
  double vpar_max_Li; // Velocity space extents in vparallel for Li ions
  double mu_max_Li; // Velocity space extents in mu for Li ions  
  int Nx; // Number of cells along x.
  int Ny; // Number of cells along y.
  int Nz; // Number of cells along z.
  int Nvpar; // Number of cells along vpar.
  int Nmu; // Number of cells along mu.
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double t_end; // end time
  int num_frames; // number of output frames
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
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
void eval_source_density_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
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
  fout[0] = 0.9*n_src*fout[0];
}
void eval_source_density_li(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
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
  fout[0] = 0.05*n_src*fout[0];
}

void eval_source_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
mapc2p_vel_power(const double* vc, double* GKYL_RESTRICT vp, double vmax, double mumax, double vpow, double mupow)
{
  double cvpar = vc[0];
  double cmu = vc[1];

  if ( cvpar < 0.0 )
    vp[0] = -fabs(vmax * pow(fabs(cvpar), vpow));
  else
    vp[0] = fabs(vmax * pow(cvpar, vpow));
  
  vp[1] = mumax * pow(cmu, mupow);
}

void
mapc2p_vel_elc(double t, const double* vc, double* GKYL_RESTRICT vp, void* ctx)
{
  struct gk_app_ctx *app = ctx;
  double vmax = app->vpar_max_elc;
  double mumax = app->mu_max_elc;
  mapc2p_vel_power(vc, vp, vmax, mumax, 1.0, 1.0);
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
void eval_density_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
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
  double nPeak = 0.9*4.*sqrt(5)/3./c_ss*Ls/2.*n_src[0];
  if (fabs(z) <= Ls) {
    fout[0] = nPeak*(1.+sqrt(1.-pow(z/Ls,2)))/2;
  } else {
    fout[0] = nPeak/2.;
  }
}
void eval_density_li(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
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
  double nPeak = 0.05*4.*sqrt(5)/3./c_ss*Ls/2.*n_src[0];
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

void evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
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
  int cdim = 3, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS; // Proton mass.
  double me = GKYL_ELECTRON_MASS; // Electron mass.

  double mi = 2.014*mp; // D ion mass
  double mLi = 6.94*mp; // Li ion mass
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Reference density and temperature.
  double Te = 100.0*eV;
  double Ti = 100.0*eV;
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
  double vtLi = sqrt(Ti/mLi);
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

  double vpar_max_Li = 4.0*vtLi;
  double mu_max_Li = (3./2.)*0.5*mi*pow(4.0*vtLi,2)/(2.0*B0);

  // Number of cells.
  int Nx = 4;
  int Ny = 1;
  int Nz = 8;
  int Nvpar = 6;
  int Nmu = 4;

  double t_end = 1.0e-7; 
  double num_frames = 1;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .massLi = mLi,
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
    .vpar_max_Li = vpar_max_Li, 
    .mu_max_Li = mu_max_Li, 
    .Nx = Nx,
    .Ny = Ny,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Ny, Nz, Nvpar, Nmu},
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
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

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

  // Create decomposition.
  struct gkyl_rect_decomp *decomp = gkyl_gyrokinetic_comms_decomp_new(ctx.cdim, cells_x, app_args.cuts, app_args.use_mpi, stderr);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, decomp, stderr);

  int my_rank = 0;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    gkyl_comm_get_rank(comm, &my_rank);
#endif

  // Electrons.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -1.0, 0.0},
    .upper = {  1.0, 1.0}, 
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

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
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources=1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_source_density,
        .ctx_upar = &ctx,
        .upar= eval_source_upar,
        .ctx_temp = &ctx,
        .temp = eval_source_temp,      
      }, 
    },
    .react = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ELC, 
          .ion_id = GKYL_ION_LI, 
          .elc_nm = "elc", 
          .ion_nm = "Li2", 
          .donor_nm = "Li1", 
          .charge_state = 1, 
          .ion_mass = ctx.massLi, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_LI,
          .elc_nm = "elc",
          .ion_nm = "Li2",
          .recvr_nm = "Li1",
          .charge_state = 1,
          .ion_mass = ctx.massLi,
          .elc_mass = ctx.massElc,
        },
      },
    },
    .bcx = {
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

  // Ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ion,
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
      .num_sources=1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_source_density_ion,
        .ctx_upar = &ctx,
        .upar= eval_source_upar,
        .ctx_temp = &ctx,
        .temp = eval_source_temp,      
      }, 
    },
    .bcx = {
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

  // Li1+ ions
  struct gkyl_gyrokinetic_species Li1 = {
    .name = "Li1",
    .charge = ctx.chargeIon, .mass = ctx.massLi,
    .lower = { -ctx.vpar_max_Li, 0.0},
    .upper = { ctx.vpar_max_Li, ctx.mu_max_Li}, 
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = 0.05*ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_li,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ion,      
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources=1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_source_density_li,
        .ctx_upar = &ctx,
        .upar= eval_source_upar,
        .ctx_temp = &ctx,
        .temp = eval_source_temp,      
      }, 
    },
    .react = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_DONOR, 
          .ion_id = GKYL_ION_LI, 
          .elc_nm = "elc", 
          .ion_nm = "Li2", 
          .donor_nm = "Li1", 
          .charge_state = 1, 
          .ion_mass = ctx.massLi, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_RECVR,
          .ion_id = GKYL_ION_LI,
          .elc_nm = "elc",
          .ion_nm = "Li2",
          .recvr_nm = "Li1",
          .charge_state = 1,
          .ion_mass = ctx.massLi,
          .elc_mass = ctx.massElc,
        },
      },
    },

    .bcx = {
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

  // Li2+ ions
  struct gkyl_gyrokinetic_species Li2 = {
    .name = "Li2",
    .charge = 2.*ctx.chargeIon, .mass = ctx.massLi,
    .lower = { -ctx.vpar_max_Li, 0.0},
    .upper = { ctx.vpar_max_Li, ctx.mu_max_Li}, 
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = 0.05*ctx.n0,

    .projection = { 
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ion,      
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources=1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_source_density_li,
        .ctx_upar = &ctx,
        .upar= eval_source_upar,
        .ctx_temp = &ctx,
        .temp = eval_source_temp,      
      }, 
    },
    .react = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ION, 
          .ion_id = GKYL_ION_LI, 
          .elc_nm = "elc", 
          .ion_nm = "Li2", 
          .donor_nm = "Li1", 
          .charge_state = 1, 
          .ion_mass = ctx.massLi, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_LI,
          .elc_nm = "elc",
          .ion_nm = "Li2",
          .recvr_nm = "Li1",
          .charge_state = 1,
          .ion_mass = ctx.massLi,
          .elc_mass = ctx.massElc,
        },
      },
    },

    .bcx = {
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

  // field
  struct gkyl_gyrokinetic_field field = {
    .fem_parbc = GKYL_FEM_PARPROJ_NONE,
    .poisson_bcs = {.lo_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_PERIODIC},
                    .up_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_PERIODIC},
                    .lo_value = {0.0, 0.0}, .up_value = {0.0, 0.0}},
  };

  // GK app
  struct gkyl_gk app_inp = {
    .name = "gk_li_react_3x2v_p1_linear",

    .cdim = 3, .vdim = 2,
    .lower = { ctx.R-ctx.Lx/2.0, -ctx.Ly/2.0, -ctx.Lz/2.0 },
    .upper = { ctx.R+ctx.Lx/2.0,  ctx.Ly/2.0,  ctx.Lz/2.0 },
    .cells = { cells_x[0], cells_x[1], cells_x[2] },
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

    .num_species = 4,
    .species = { elc, ion, Li1, Li2 },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

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

  // initial time-step
  double dt = t_end-t_curr;

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

  gkyl_gyrokinetic_app_stat_write(app);
  
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

  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  gkyl_gyrokinetic_comms_release(decomp, comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
