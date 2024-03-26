#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_util.h>
#include <gkyl_math.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct ltx_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using non-normalized physical units).
  double epsilon0; // Permittivity of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.

  double Te0; // Reference electron temperature.
  double Ti0; // Reference ion temperature.
  double n0; // Reference number density (1 / m^3).

  double R_axis; // Major radius of the magnetic axis (m).
  double B_axis; // Magnetic field strength at the magnetic axis (Tesla).
  double R_LCFSmid; // Major radius of the Last Closed Flux Surface at the outboard midplane (m).
  double q_sep; // Safety factor at the separatrix.
  double s_sep; // Magnetic shear at the separatrix.
  double kappa; // Elongation (=1 for no elongation).
  double delta; // Triangularity (=0 for no triangularity).

  double nu_frac; // Collision frequency fraction.

  double k_perp_rho_si; // Product of perpendicular wavenumber and ion-sound gyroradius.

  double n_src; // Source number density.
  double Te_src; // Source electron temperature.
  double Ti_src; // Source ion temperature.

  double perturb; // Perturbation amplitude.

  // Derived physical quantities (using non-normalized physical units).
  double R_ref; // Major radius of the domain (m).
  double a_mid; // Minor radius at the outboard midplane (m).
  double r_ref; // Minor radius of the domain (m).
  double B_ref; // Magnetic field strength in the center of the domain (Tesla).

  double q_ref; // Magnetic safety factor in the center of the domain.
  double epsilon_ref; // Inverse aspect ratio of the domain.

  double log_lambda_elc; // Logarithm of electron wavelength.
  double nu_elc; // Electron collision frequency.
  double log_lambda_ion; // Logarithm of ion wavelength.
  double nu_ion; // Ion collision frequency.

  double c_s; // Sound speed.
  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.
  double omega_ci; // Ion cyclotron frequency.
  double rho_si; // Ion-sound gyroradius.

  double k_perp; // Perpendicular wavenumber (for Poisson solver).

  double c_s_src; // Source sound speed.
  double n_peak; // Peak number density.

  // Simulation parameters.
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  double Lz; // Domain size (configuration space: z-direction).
  double Lvpar_elc; // Domain size (electron velocity space: parallel velocity direction).
  double Lmu_elc; // Domain size (electron velocity space: magnetic moment direction).
  double Lvpar_ion; // Domain size (ion velocity space: parallel velocity direction).
  double Lmu_ion; // Domain size (ion velocity space: magnetic moment direction).
  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
};

struct ltx_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using non-normalized physical units).
  double epsilon0 = GKYL_EPSILON0; // Permittivity of free space.
  double mass_elc = GKYL_ELECTRON_MASS; // Electron mass.
  double charge_elc = -GKYL_ELEMENTARY_CHARGE; // Electron charge.
  double mass_ion = GKYL_PROTON_MASS; // Proton mass.
  double charge_ion = GKYL_ELEMENTARY_CHARGE; // Proton charge.

  double Te0 = 178.0 * GKYL_ELEMENTARY_CHARGE; // Reference electron temperature.
  double Ti0 = 70.0 * GKYL_ELEMENTARY_CHARGE; // Reference ion temperature.
  double n0 = 1.78e18; // Reference number density (1 / m^3).

  double R_axis = 0.4; // Major radius of the magnetic axis (m).
  double B_axis = 0.224374548; // Magnetic field strength at the magnetic axis (Tesla).
  double R_LCFSmid = 0.5948; // Major radius of the Last Closed Flux Surface at the outboard midplane (m).
  double q_sep = 3.69546081; // Safety factor at the separatrix.
  double s_sep = 2.27976219; // Magnetic shear at the separatrix.
  double kappa = 1.57; // Elongation (=1 for no elongation).
  double delta = 0.6; // Triangularity (=0 for no triangularity).

  double nu_frac = 1.0; // Collision frequency fraction.

  double k_perp_rho_si = 0.15; // Product of perpendicular wavenumber and ion-sound gyroradius.

  double n_src = 1.95e22; // Source number density.
  double Te_src = 410.0 * GKYL_ELEMENTARY_CHARGE; // Source electron temperature.
  double Ti_src = 40.0 * GKYL_ELEMENTARY_CHARGE; // Source ion temperature.

  double perturb = 0.0; // Perturbation amplitude.

  // Derived physical quantities (using non-normalized physical units).
  double R_ref = R_LCFSmid + 0.025; // Major radius of the domain (m).
  double a_mid = R_LCFSmid - R_axis; // Minor radius at the outboard midplane (m).
  double r_ref = R_ref - R_axis; // Minor radius of the domain (m).
  double B_ref = B_axis * (R_axis / R_ref); // Magnetic field strength in the center of the domain (Tesla).

  double q_ref = q_sep / (1.0 - s_sep * ((0.0 + r_ref - a_mid) / a_mid)); // Magnetic safety factor in the center of the domain.
  double epsilon_ref = r_ref / R_ref; // Inverse aspect ratio of the domain.

  double log_lambda_elc = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Te0 / charge_ion); // Logarithm of electron wavelength.
  double nu_elc = nu_frac * log_lambda_elc * (charge_ion * charge_ion * charge_ion * charge_ion) * n0 /
    (6.0 * sqrt(2.0) * pi * sqrt(pi) * epsilon0 * epsilon0 * sqrt(mass_elc) * (Te0 * sqrt(Te0))); // Electron collision frequency.
  double log_lambda_ion = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Ti0 / charge_ion); // Logarithm of ion wavelength.
  double nu_ion = nu_frac * log_lambda_ion * (charge_ion * charge_ion * charge_ion * charge_ion) * n0 /
    (12.0 * pi * sqrt(pi) * epsilon0 * epsilon0 * sqrt(mass_ion) * (Ti0 * sqrt(Ti0))); // Ion collision frequency.
  
  double c_s = sqrt(Te0 / mass_ion); // Sound speed.
  double vte = sqrt(Te0 / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti0 / mass_ion); // Ion thermal velocity.
  double omega_ci = fabs(charge_ion * B_ref / mass_ion); // Ion cyclotron frequency.
  double rho_si = c_s / omega_ci; // Ion-sound gyroradius.

  double k_perp = k_perp_rho_si / rho_si; // Perpendicular wavenumber (for Poisson solver).

  double c_s_src = sqrt(5.0 / 3.0 * Te_src / mass_ion); // Source sound speed.
  double n_peak = 89.0 * sqrt(5.0) / 3.0 / c_s_src * (0.62 * 2.0 * 0.25 * pi) * 0.5 * n_src; // Peak number density.

  // Simulation parameters.
  int Nz = 64; // Cell count (configuration space: z-direction).
  int Nvpar = 16; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 45; // Cell count (velocity space: magnetic moment direction).
  double Lz = 0.62 * 2.0 * pi; // Domain size (configuration space: z-direction).
  double Lvpar_elc = 8.0 * vte; // Domain size (electron velocity space: parallel velocity direction).
  double Lmu_elc = mass_elc * (1.5 * 4.0 * vte) * (1.5 * 4.0 * vte) / (2.0 * B_ref); // Domain size (electron velocity space: magnetic moment direction).
  double Lvpar_ion = 8.0 * vti; // Domain size (ion velocity space: parallel velocity direction).
  double Lmu_ion = mass_ion * (1.5 * 4.0 * vti) * (1.5 * 4.0 * vti) / (2.0 * B_ref); // Domain size (ion velocity space: magnetic moment direction).
  double t_end = 1.0e-7; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  
  struct ltx_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .Te0 = Te0,
    .Ti0 = Ti0,
    .n0 = n0,
    .R_axis = R_axis,
    .B_axis = B_axis,
    .R_LCFSmid = R_LCFSmid,
    .q_sep = q_sep,
    .s_sep = s_sep,
    .kappa = kappa,
    .delta = delta,
    .nu_frac = nu_frac,
    .k_perp_rho_si = k_perp_rho_si,
    .n_src = n_src,
    .Te_src = Te_src,
    .Ti_src = Ti_src,
    .perturb = perturb,
    .R_ref = R_ref,
    .a_mid = a_mid,
    .r_ref = r_ref,
    .B_ref = B_ref,
    .q_ref = q_ref,
    .epsilon_ref = epsilon_ref,
    .log_lambda_elc = log_lambda_elc,
    .nu_elc = nu_elc,
    .log_lambda_ion = log_lambda_ion,
    .nu_ion = nu_ion,
    .c_s = c_s,
    .vte = vte,
    .vti = vti,
    .omega_ci = omega_ci,
    .rho_si = rho_si,
    .k_perp = k_perp,
    .c_s_src = c_s_src,
    .n_peak = n_peak,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .Lz = Lz,
    .Lvpar_elc = Lvpar_elc,
    .Lmu_elc = Lmu_elc,
    .Lvpar_ion = Lvpar_ion,
    .Lmu_ion = Lmu_ion,
    .t_end = t_end,
    .num_frames = num_frames,
  };

  return ctx;
}

void
evalSourceDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ltx_ctx *app = ctx;
  double z = xn[0];

  double n_src = app -> n_src;
  double Lz = app -> Lz;

  double n = 0.0;

  if (fabs(z) < 0.25 * Lz) {
    n = n_src;
  } 
  else {
    n = 1.0e-40;
  }

  // Set source number density.
  fout[0] = n;
}

void
evalSourceUparElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set electron source parallel velocity.
  fout[0] = 0.0;
}

void
evalSourceTempElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ltx_ctx *app = ctx;

  double Te_src = app -> Te_src;

  // Set electron source temperature.
  fout[0] = Te_src;
}

void
evalSourceUparIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set ion source parallel velocity.
  fout[0] = 0.0;
}

void
evalSourceTempIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ltx_ctx *app = ctx;

  double Ti_src = app -> Ti_src;

  // Set ion source temperature.
  fout[0] = Ti_src;
}

void
evalDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ltx_ctx *app = ctx;
  double z = xn[0];

  double n_peak = app -> n_peak;
  double Lz = app -> Lz;

  double n = 0.0;

  if (fabs(z) <= 0.25 * Lz) {
    n = 0.5 * n_peak * (1.0 + sqrt(1.0 - (z / (0.25 * Lz)) * (z / (0.25 * Lz))));
  }
  else {
    n = 0.5 * n_peak;
  }

  // Set number density.
  fout[0] = n;
}

void
evalUparElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set electron parallel velocity.
  fout[0] = 0.0;
}

void
evalTempElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ltx_ctx *app = ctx;

  double Te0 = app -> Te0;

  // Set electron temperature.
  fout[0] = Te0;
}

void
evalUparIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set ion parallel velocity.
  fout[0] = 0.0;
}

void
evalTempIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ltx_ctx *app = ctx;

  double Ti0 = app -> Ti0;

  // Set ion temperature.
  fout[0] = Ti0;
}

void
evalNuElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ltx_ctx *app = ctx;

  double nu_elc = app -> nu_elc;

  // Set electron collision frequency.
  fout[0] = nu_elc;
}

void
evalNuIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ltx_ctx *app = ctx;

  double nu_ion = app -> nu_ion;

  // Set ion collision frequency.
  fout[0] = nu_ion;
}

// Major radius R, as a function of minor radius r and poloidal angle theta.
double
R_rtheta(double r, double theta, void* ctx)
{
  struct ltx_ctx *app = ctx;

  double R_axis = app -> R_axis;
  double delta = app -> delta;

  // Major radius.
  return R_axis + r * cos(theta + asin(delta) * sin(theta));
}

// Height Z, as a function of minor radius r and poloidal angle theta.
double
Z_rtheta(double r, double theta, void *ctx)
{
  struct ltx_ctx *app = ctx;

  double kappa = app -> kappa;

  // Height.
  return kappa * r * sin(theta);
}

// Partial derivative of major radius R, with respect to minor radius r.
double
dR_dr(double r, double theta, void* ctx)
{
  struct ltx_ctx *app = ctx;

  double delta = app -> delta;

  // dR / dr.
  return cos(theta + asin(delta) * sin(theta));
}

// Partial derivative of major radius R, with respect to poloidal angle theta.
double
dR_dtheta(double r, double theta, void* ctx)
{
  struct ltx_ctx *app = ctx;

  double delta = app -> delta;

  // dR / dtheta.
  return -r * sin(theta + asin(delta) * sin(theta)) * (1.0 + asin(delta) * cos(theta));
}

// Partial derivative of height Z, with respect to minor radius r.
double
dZ_dr(double r, double theta, void* ctx)
{
  struct ltx_ctx *app = ctx;

  double kappa = app -> kappa;

  // dZ / dr.
  return kappa * sin(theta);
}

// Partial derviative of height Z, with respect to poloidal angle theta.
double
dZ_dtheta(double r, double theta, void* ctx)
{
  struct ltx_ctx *app = ctx;

  double kappa = app -> kappa;

  // dZ / dtheta.
  return kappa * r * cos(theta);
}

double
Jr(double r, double theta, void* ctx)
{
  return R_rtheta(r, theta, ctx) * (dR_dr(r, theta, ctx) * dZ_dtheta(r, theta, ctx) - dR_dtheta(r, theta, ctx) * dZ_dr(r, theta, ctx));
}

struct integrand_ctx
{
  struct ltx_ctx *app_ctx;
  double r;
};

double
integrand(double t, void* int_ctx)
{
  struct integrand_ctx *inctx = int_ctx;

  double r = inctx -> r;
  struct ltx_ctx *app = inctx -> app_ctx;

  return Jr(r, t, app) / (R_rtheta(r, t, app) * R_rtheta(r, t, app));
}

double
dPsi_dr(double r, double theta, void* ctx)
{
  struct ltx_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {
    .app_ctx = app,
    .r = r,
  };

  double pi = app -> pi;

  double R_axis = app -> R_axis;
  double q_sep = app -> q_sep;
  double s_sep = app -> s_sep;

  double a_mid = app -> a_mid;
  double B_ref = app -> B_ref;

  struct gkyl_qr_res integral;
  integral = gkyl_dbl_exp(integrand, &tmp_ctx, 0.0, 2.0 * pi, 7, 1.0e-10);

  return ( B_ref * R_axis / (2.0 * pi * (q_sep / (1.0 - s_sep * ((r - a_mid) / a_mid))))) * integral.res;
}

double
alpha(double r, double theta, double phi, void* ctx)
{
  struct ltx_ctx *app = ctx;

  double pi = app -> pi;

  double R_axis = app -> R_axis;
  double B_ref = app -> B_ref;

  double twrap = theta;

  while (twrap < -pi) {
    twrap = twrap + 2.0 * pi;
  }
  while (pi < twrap) {
    twrap = twrap - 2.0 * pi;
  }

  struct integrand_ctx tmp_ctx = {
    .app_ctx = app,
    .r = r
  };
  struct gkyl_qr_res integral;

  if (0.0 < twrap) {
    integral = gkyl_dbl_exp(integrand, &tmp_ctx, 0.0, twrap, 7, 1.0e-10);
  }
  else {
    integral = gkyl_dbl_exp(integrand, &tmp_ctx, twrap, 0.0, 7, 1.0e-10);
    integral.res = -integral.res;
  }

  return phi - B_ref * R_axis * integral.res / dPsi_dr(r, theta, ctx);
}

double
Bphi(double R, void* ctx)
{
  struct ltx_ctx *app = ctx;

  double R_ref = app -> R_ref;
  double B_ref = app -> B_ref;

  return B_ref * R_ref / R;
}

double
grad_r(double r, double theta, void* ctx)
{
  return (R_rtheta(r, theta, ctx) / Jr(r, theta, ctx)) * sqrt((dR_dtheta(r, theta, ctx) * dR_dtheta(r, theta, ctx)) + (dZ_dtheta(r, theta, ctx) * dZ_dtheta(r, theta, ctx)));
}

static inline void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  struct ltx_ctx *app = ctx;
  double x = zc[0], y = zc[1], z = zc[2];

  double kappa = app -> kappa;
  double r_ref = app -> r_ref;
  double q_ref = app -> q_ref;

  double r = x + r_ref;

  // Set cylindrical coordinates (R, Z, phi) from computational coordinates (x, y, z).
  double R = R_rtheta(r, x, ctx);
  double Z = kappa * r * sin(z);
  double phi = - q_ref / r_ref * y - alpha(r, z, 0.0, ctx);

  // Set physical coordinates (X, Y, Z) from cylindrical coordiantes (R, Z, phi).
  xp[0] = R * cos(phi);
  xp[1] = R * sin(phi);
  xp[2] = Z;
}

void
bmag_func(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void* ctx)
{
  struct ltx_ctx *app = ctx;
  double z = zc[2];

  double r_ref = app -> r_ref;

  double Bt = Bphi(R_rtheta(r_ref, z, ctx), ctx);
  double Bp = dPsi_dr(r_ref, z, ctx) / R_rtheta(r_ref, z, ctx) * grad_r(r_ref, z, ctx);

  // Set magnetic field strength.
  fout[0] = sqrt(Bt * Bt + Bp * Bp);
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    gkyl_gyrokinetic_app_write(app, t_curr, iot -> curr - 1);
    gkyl_gyrokinetic_app_calc_mom(app);
    gkyl_gyrokinetic_app_write_mom(app, t_curr, iot -> curr - 1);
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

  struct ltx_ctx ctx = create_ctx(); // Context for initialization functions.

  int NZ = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nz);
  int NVPAR = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvpar);
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nmu);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif  

  // Create global range.
  int ccells[] = { NZ };
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

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
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

  // Electron species.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { -0.5 * ctx.Lvpar_elc, 0.0},
    .upper = { 0.5 * ctx.Lvpar_elc, ctx.Lmu_elc},
    .cells = { NVPAR, NMU },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .density = evalDensityInit,
      .ctx_density = &ctx,
      .upar = evalUparElcInit,
      .ctx_upar = &ctx,
      .temp = evalTempElcInit,
      .ctx_temp = &ctx,
    },
    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalNuElcInit,
      .ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = evalSourceDensityInit,
        .ctx_density = &ctx,
        .upar= evalSourceUparElcInit,
        .ctx_upar = &ctx,
        .temp = evalSourceTempElcInit,
        .ctx_temp = &ctx,
      }, 
    },

    .bcx = {
      .lower = { .type = GKYL_SPECIES_GK_SHEATH, },
      .upper = { .type = GKYL_SPECIES_GK_SHEATH, },
    },

    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ion species.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .lower = { -0.5 * ctx.Lvpar_ion, 0.0},
    .upper = { 0.5 * ctx.Lvpar_ion, ctx.Lmu_ion},
    .cells = { NVPAR, NMU },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .density = evalDensityInit,
      .ctx_density = &ctx,
      .upar = evalUparIonInit,
      .ctx_upar = &ctx,
      .temp = evalTempIonInit,
      .ctx_temp = &ctx,
    },
    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalNuIonInit,
      .ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = evalSourceDensityInit,
        .ctx_density = &ctx,
        .upar = evalSourceUparIonInit,
        .ctx_upar = &ctx,
        .temp = evalSourceTempIonInit,
        .ctx_temp = &ctx,
      }, 
    },

    .bcx = {
      .lower = { .type = GKYL_SPECIES_GK_SHEATH, },
      .upper = { .type = GKYL_SPECIES_GK_SHEATH, },
    },
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .bmag_fac = ctx.B_ref,
    .fem_parbc = GKYL_FEM_PARPROJ_NONE,
    .kperpSq = ctx.k_perp * ctx.k_perp,
  };

  // GK app.
  struct gkyl_gk app_inp = {
    .name = "gk_ltx_1x2v_p1",

    .cdim = 1, .vdim = 2,
    .lower = { -0.5 * ctx.Lz },
    .upper = { 0.5 * ctx.Lz },
    .cells = { NZ },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .world = { 0.0, 0.0 },
      .mapc2p = mapc2p,
      .c2p_ctx = &ctx,
      .bmag_func = bmag_func,
      .bmag_ctx = &ctx
    },

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .use_gpu = app_args.use_gpu,

    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp -> ranges[my_rank],
      .comm = comm
    }
  };

  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // Initialize simulation.
  gkyl_gyrokinetic_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr);

  gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
  gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, t_curr);

    step += 1;
  }

  gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
  gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
  
  write_data(&io_trig, app, t_curr);
  gkyl_gyrokinetic_app_stat_write(app);
  
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);

  mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
