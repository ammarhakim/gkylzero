#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_math.h>

#include <rt_arg_parse.h>

// Define the context of the simulation. This is basically all the globals
struct gk_mirror_ctx
{
  int cdim, vdim; // Dimensionality.

  // Plasma parameters
  double mi;
  double qi;
  double me;
  double qe;
  double Te0;
  double n0;
  double B_p;
  double beta;
  double tau;
  double Ti0;
  double kperpRhos;
  // Electron-electron collision freq.
  double logLambdaElc;
  double nuElc;
  // Ion-ion collision freq.
  double logLambdaIon;
  double nuIon;
  // Thermal speeds.
  double vti;
  double vte;
  double c_s;
  // Gyrofrequencies and gyroradii.
  double omega_ci;
  double rho_s;
  double kperp; // Perpendicular wavenumber in SI units.
  // Axial coordinate Z extents. Endure that Z=0 is not on
  double z_min;
  double z_max;
  double psi_min;
  double psi_max;
  double z_m;
  // Physics parameters at mirror throat
  double n_m;
  double Te_m;
  double Ti_m;
  double Ti_perp0;
  double Ti_par0;
  double Ti_perp_m;
  double Ti_par_m;
  double Te_perp0;
  double Te_par0;
  double Te_perp_m;
  double Te_par_m;
  double cs_m;
  // Source parameters
  double NSrcIon;
  double lineLengthSrcIon;
  double sigSrcIon;
  double NSrcFloorIon;
  double TSrc0Ion;
  double TSrcFloorIon;
  double NSrcElc;
  double lineLengthSrcElc;
  double sigSrcElc;
  double NSrcFloorElc;
  double TSrc0Elc;
  double TSrcFloorElc;
  // Grid parameters
  double vpar_max_ion;
  double vpar_max_elc;
  double mu_max_ion;
  double mu_max_elc;
  int Nx;
  int Ny;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

// -- Source functions.
void
eval_density_elc_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[2];
  double NSrc = app->NSrcElc;
  double zSrc = app->lineLengthSrcElc;
  double sigSrc = app->sigSrcElc;
  double NSrcFloor = app->NSrcFloorElc;
  if (fabs(z) <= app->z_m)
  {
    fout[0] = fmax(NSrcFloor, (NSrc / sqrt(2.0 * M_PI * pow(sigSrc, 2.))) *
                                  exp(-1 * pow((z - zSrc), 2) / (2.0 * pow(sigSrc, 2.))));
  }
  else
  {
    fout[0] = 1e-16;
  }
}

void
eval_upar_elc_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[2];
  double sigSrc = app->sigSrcElc;
  double TSrc0 = app->TSrc0Elc;
  double Tfloor = app->TSrcFloorElc;
  if (fabs(z) <= 2.0 * sigSrc)
  {
    fout[0] = TSrc0;
  }
  else
  {
    fout[0] = Tfloor;
  }
}

void
eval_density_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[2];
  double NSrc = app->NSrcIon;
  double zSrc = app->lineLengthSrcIon;
  double sigSrc = app->sigSrcIon;
  double NSrcFloor = app->NSrcFloorIon;
  if (fabs(z) <= app->z_m)
  {
    fout[0] = fmax(NSrcFloor, (NSrc / sqrt(2.0 * M_PI * pow(sigSrc, 2))) *
                                  exp(-1 * pow((z - zSrc), 2) / (2.0 * pow(sigSrc, 2))));
  }
  else
  {
    fout[0] = 1e-16;
  }
}

void
eval_upar_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[2];
  double sigSrc = app->sigSrcIon;
  double TSrc0 = app->TSrc0Ion;
  double Tfloor = app->TSrcFloorIon;
  if (fabs(z) <= 2.0 * sigSrc)
  {
    fout[0] = TSrc0;
  }
  else
  {
    fout[0] = Tfloor;
  }
}

// Potential initial condition
void
eval_potential(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[2];
  double z_m = app->z_m;
  double sigma = 0.9*z_m;
  double center_potential = 5.0 * app->Te0 / app->qi;
  if (fabs(z) <= sigma)
  {
    fout[0] = 0.5*center_potential*(1. + tanh(10. * sigma * fabs(sigma - fabs(z))));
  }
  else
  {
    fout[0] = 0.0;
  }
}

// Electrons initial conditions
void
eval_density_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[2];
  double z_m = app->z_m;
  double sigma = 0.9*z_m;
  if (fabs(z) <= sigma)
  {
    fout[0] = 0.5*app->n0*(1. + tanh(10. * sigma * fabs(sigma - fabs(z))));
  }
  else
  {
    fout[0] = 0.5*app->n0*exp(-5 * (fabs(sigma - fabs(z))));
  }
}

void
eval_upar_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[2];
  double cs_m = app->cs_m;
  double z_m = app->z_m;
  double z_max = app->z_max;
  if (fabs(z) <= z_m)
  {
    fout[0] = 0.0;
  }
  else
  {
    fout[0] = fabs(z) / z * cs_m * tanh(3 * (z_max - z_m) * fabs(fabs(z) - z_m)); // Maybe put a 5 here
  }
}

void
eval_temp_par_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[2];
  double z_m = app->z_m;
  double Te_par0 = app->Te_par0;
  double Te_par_m = app->Te_par_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Te_par_m+(Te_par0-Te_par_m)*tanh(4 * fabs(z_m - fabs(z)));
  }
  else
  {
    fout[0] = Te_par_m;
  }
}

void
eval_temp_perp_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[2];
  double z_m = app->z_m;
  double Te_perp0 = app->Te_perp0;
  double Te_perp_m = app->Te_perp_m;
  if (fabs(z) <= z_m)
  {
   fout[0] = Te_perp_m - Te_perp0*tanh(3.*fabs(z_m-fabs(z)));
  }
  else
  {
    fout[0] = Te_perp_m * GKYL_MAX2(1.e-3, exp(-5. * (fabs(z_m - fabs(z)))));
  }
}

// Ion initial conditions
void
eval_density_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[2];
  double z_m = app->z_m;
  double sigma = 0.9*z_m;
  if (fabs(z) <= sigma)
  {
    fout[0] = 0.5*app->n0*(1. + tanh(10. * sigma * fabs(sigma - fabs(z))));
  }
  else
  {
    fout[0] = 0.5*app->n0* exp(-5 * (fabs(sigma - fabs(z))));
  }
}

void
eval_upar_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[2];
  double cs_m = app->cs_m;
  double z_m = app->z_m;
  double z_max = app->z_max;
  if (fabs(z) <= z_m)
  {
    fout[0] = 0.0;
  }
  else
  {
    fout[0] = fabs(z) / z * cs_m * tanh(3 * (z_max - z_m) * fabs(fabs(z) - z_m)); // Maybe put a 5 here
  }
}
void
eval_temp_par_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[2];
  double z_m = app->z_m;
  double Ti_par0 = app->Ti_par0;
  double Ti_par_m = app->Ti_par_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Ti_par_m + (Ti_par0 - Ti_par_m) * tanh(4 * fabs(z_m - fabs(z)));
  }
  else
  {
    fout[0] = Ti_par_m * GKYL_MAX2(1.e-2, 4 * log(fabs(fabs(z) - z_m) + 1));
  }
}

void
eval_temp_perp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[2];
  double z_m = app->z_m;
  double Ti_perp0 = app->Ti_perp0;
  double Ti_perp_m = app->Ti_perp_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Ti_perp_m + (Ti_perp0 - Ti_perp_m) * tanh(3. * fabs(z_m - fabs(z)));
  }
  else
  {
    fout[0] = Ti_perp_m * GKYL_MAX2(1.e-3, exp(-5. * (fabs(z_m - fabs(z)))));
  }
}

// Evaluate collision frequencies
void
evalNuElc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->nuIon;
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  double b = 1.45;
  double linear_velocity_threshold = 1./3.;
  double frac_linear = 1/b*atan(linear_velocity_threshold*tan(b));
  if (fabs(cvpar) < frac_linear) {
    double func_frac = tan(frac_linear*b) / tan(b);
    vp[0] = vpar_max_ion*func_frac*cvpar/frac_linear;
  }
  else {
    vp[0] = vpar_max_ion*tan(cvpar*b)/tan(b);
  }
  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}

void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;

  double cvpar = vc[0], cmu = vc[1];
  double b = 1.45;
  double linear_velocity_threshold = 1./6.;
  double frac_linear = 1/b*atan(linear_velocity_threshold*tan(b));
  if (fabs(cvpar) < frac_linear) {
    double func_frac = tan(frac_linear*b) / tan(b);
    vp[0] = vpar_max_elc*func_frac*cvpar/frac_linear;
  }
  else {
    vp[0] = vpar_max_elc*tan(cvpar*b)/tan(b);
  }
  // Quadratic map in mu.
  vp[1] = mu_max_elc*pow(cmu,2);
}

struct gk_mirror_ctx
create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0;
  double mu0 = GKYL_MU0; // Not sure if this is right
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV;  // ion charge
  double qe = -eV; // electron charge

  // Plasma parameters.
  double mi = 2.014 * mp;
  double Te0 = 940 * eV;
  double n0 = 3e19;
  double B_p = 0.53;
  double beta = 0.4;
  double tau = pow(B_p, 2.) * beta / (2.0 * mu0 * n0 * Te0) - 1.;
  double Ti0 = tau * Te0;
  double kperpRhos = 0.1;

  double nuFrac = 1.0;
  // Electron-electron collision freq.
  double logLambdaElc = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Te0 / eV);
  double nuElc = nuFrac * logLambdaElc * pow(eV, 4.) * n0 /
                 (6. * sqrt(2.) * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(me) * pow(Te0, 3. / 2.));
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Ti0 / eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4.) * n0 /
                 (12 * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(mi) * pow(Ti0, 3. / 2.));

  // Thermal speeds.
  double vti = sqrt(Ti0 / mi);
  double vte = sqrt(Te0 / me);
  double c_s = sqrt(Te0 / mi);

  // Gyrofrequencies and gyroradii.
  double omega_ci = eV * B_p / mi;
  double rho_s = c_s / omega_ci;

  // Perpendicular wavenumber in SI units:
  double kperp = kperpRhos / rho_s;

  // Geometry parameters.
  double z_min = -2.0;
  double z_max =  2.0;
  double z_m = 0.982544;
  double psi_min = 1e-5;
  double psi_max = 1e-3;

  // Source parameters
  double NSrcIon = 3.1715e23 / 8.0;
  double lineLengthSrcIon = 0.0;
  double sigSrcIon = z_m / 4.0;
  double NSrcFloorIon = 0.05 * NSrcIon;
  double TSrc0Ion = Ti0 * 1.25;
  double TSrcFloorIon = TSrc0Ion / 8.0;
  double NSrcElc = NSrcIon;
  double lineLengthSrcElc = lineLengthSrcIon;
  double sigSrcElc = sigSrcIon;
  double NSrcFloorElc = NSrcFloorIon;
  double TSrc0Elc = TSrc0Ion / tau;
  double TSrcFloorElc = TSrcFloorIon / tau;

  // Initial conditions parameters
  double Ti_perp0 = 10000 * eV;
  double Ti_perp_m = 15000 * eV;
  double Ti_par0 = 7500 * eV;
  double Ti_par_m = 1000 * eV;

  double Te_par0 = 1800 * eV;  
  double Te_par_m = 300 * eV;
  double Te_perp0 = 2000 * eV;
  double Te_perp_m = 3000 * eV;

  // Physics parameters at mirror throat
  double n_m = 1.105617e19;
  double Te_m = 346.426583 * eV;
  double Ti_m = 3081.437703 * eV;
  double cs_m = 4.037740e5;

  // Grid parameters
  double vpar_max_elc = 20 * vte;
  double mu_max_elc = me * pow(3. * vte, 2.) / (2. * B_p);
  double vpar_max_ion = 20 * vti;
  double mu_max_ion = mi * pow(3. * vti, 2.) / (2. * B_p);
  int Nx = 2;
  int Ny = 2;
  int Nz = 32;
  int Nvpar = 10; // Number of cells in the paralell velocity direction 96
  int Nmu = 10;  // Number of cells in the mu direction 192
  int poly_order = 1;

  double t_end = 1e-9;
  int num_frames = 1;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_mirror_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .mi = mi,
    .qi = qi,
    .me = me,
    .qe = qe,
    .Te0 = Te0,
    .n0 = n0,
    .B_p = B_p,
    .beta = beta,
    .tau = tau,
    .Ti0 = Ti0,
    .kperpRhos = kperpRhos,
    .logLambdaElc = logLambdaElc,
    .nuElc = nuElc,
    .logLambdaIon = logLambdaIon,
    .nuIon = nuIon,
    .vti = vti,
    .vte = vte,
    .c_s = c_s,
    .omega_ci = omega_ci,
    .rho_s = rho_s,
    .kperp = kperp, 
    .z_min = z_min,
    .z_max = z_max,
    .z_m = z_m,
    .psi_min = psi_min,
    .psi_max = psi_max,
    .Ti_perp0 = Ti_perp0,
    .Ti_par0 = Ti_par0,
    .Ti_perp_m = Ti_perp_m,
    .Ti_par_m = Ti_par_m,
    .Te_par0 = Te_par0,
    .Te_par_m = Te_par_m,
    .Te_perp0 = Te_perp0,
    .Te_perp_m = Te_perp_m,
    .cs_m = cs_m,
    .NSrcIon = NSrcIon,
    .lineLengthSrcIon = lineLengthSrcIon,
    .sigSrcIon = sigSrcIon,
    .NSrcFloorIon = NSrcFloorIon,
    .TSrc0Ion = TSrc0Ion,
    .TSrcFloorIon = TSrcFloorIon,
    .NSrcElc = NSrcElc,
    .lineLengthSrcElc = lineLengthSrcElc,
    .sigSrcElc = sigSrcElc,
    .NSrcFloorElc = NSrcFloorElc,
    .TSrc0Elc = TSrc0Elc,
    .TSrcFloorElc = TSrcFloorElc,
    .vpar_max_ion = vpar_max_ion,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_ion = mu_max_ion,
    .mu_max_elc = mu_max_elc,
    .Nx = Nx,
    .Ny = Ny,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Ny, Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .t_end = t_end,
    .num_frames = num_frames,
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app,
  double t_curr, bool is_restart_IC, bool force_calc, double dt)
{
  if (!is_restart_IC && (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc)) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);

    if ( !(dt < 0.0) )
      gkyl_gyrokinetic_app_save_dt(app, t_curr, dt);
  }
}

void
write_data(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_app* app, double t_curr, bool is_restart_IC, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;
    gkyl_gyrokinetic_app_write_conf(app, t_curr, frame);

    if (!is_restart_IC) {
      gkyl_gyrokinetic_app_write_field_energy(app);
      gkyl_gyrokinetic_app_write_integrated_mom(app);
      gkyl_gyrokinetic_app_write_dt(app);
    }
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

    gkyl_gyrokinetic_app_write_phase(app, t_curr, frame);
  }
}

int main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_mirror_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  struct gkyl_gyrokinetic_projection elc_ic = {
    .proj_id = GKYL_PROJ_BIMAXWELLIAN, 
    .ctx_density = &ctx,
    .density = eval_density_elc,
    .ctx_upar = &ctx,
    .upar= eval_upar_elc,
    .ctx_temppar = &ctx,
    .temppar = eval_temp_par_elc,      
    .ctx_tempperp = &ctx,
    .tempperp = eval_temp_perp_elc,   
  };

  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe,
    .mass = ctx.me,
    .lower = {-1.0, 0.0},
    .upper = { 1.0, 1.0},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,
    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },
    .projection = elc_ic,
    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_elc_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_elc_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_elc_source,      
      }, 
    },
    .bcx = {
      .lower = {
        .type = GKYL_SPECIES_FIXED_FUNC,
        .projection = elc_ic,
      },
      .upper = {
        .type = GKYL_SPECIES_FIXED_FUNC,
        .projection = elc_ic,
      },
    },
    .bcz = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    .num_diag_moments = 8,
    .diag_moments = {GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
  };

  struct gkyl_gyrokinetic_projection ion_ic = {
    .proj_id = GKYL_PROJ_BIMAXWELLIAN, 
    .ctx_density = &ctx,
    .density = eval_density_ion,
    .ctx_upar = &ctx,
    .upar= eval_upar_ion,
    .ctx_temppar = &ctx,
    .temppar = eval_temp_par_ion,      
    .ctx_tempperp = &ctx,
    .tempperp = eval_temp_perp_ion,   
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .lower = {-1.0, 0.0},
    .upper = { 1.0, 1.0},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,
    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },
    .projection = ion_ic,
    .scale_with_polarization = true,
    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_ion_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_ion_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_ion_source,      
      }, 
    },
    .bcx = {
      .lower = {
        .type = GKYL_SPECIES_FIXED_FUNC,
        .projection = ion_ic,
      },
      .upper = {
        .type = GKYL_SPECIES_FIXED_FUNC,
        .projection = ion_ic,
      },
    },
    .bcz = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    .num_diag_moments = 8,
    .diag_moments = {GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
  };
  struct gkyl_gyrokinetic_field field =
  {
    .polarization_bmag = ctx.B_p,
    .poisson_bcs = {
      .lo_type = {GKYL_POISSON_NEUMANN, GKYL_POISSON_PERIODIC},
      .up_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_PERIODIC},
      .lo_value = {0.0, 0.0},
      .up_value = {0.0, 0.0},
    },
    .polarization_potential = eval_potential,
    .polarization_potential_ctx = &ctx,
  };

  struct gkyl_mirror_geo_grid_inp grid_inp = {
    .filename_psi = "core/data/unit/wham_hires.geqdsk_psi.gkyl", // psi file to use
    .rclose = 0.2, // closest R to region of interest
    .zmin = -2.0,  // Z of lower boundary
    .zmax =  2.0,  // Z of upper boundary 
    .include_axis = false, // Include R=0 axis in grid
    .fl_coord = GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z, // coordinate system for psi grid
  };

  struct gkyl_gk app_inp = {
    .name = "gk_wham_nonuniformx_3x2v_p1_numeric",
    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .lower = {ctx.psi_min, - M_PI, ctx.z_min},
    .upper = {ctx.psi_max,   M_PI, ctx.z_max},
    .cells = { cells_x[0], cells_x[1], cells_x[2] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .geometry = {
      .geometry_id = GKYL_MIRROR,
      .mirror_grid_info = grid_inp,
      .position_map_info = {
        .id = GKYL_PMAP_CONSTANT_DB_NUMERIC,
        .map_strength = 0.666,
        .maximum_slope_at_max_B = 0.6,
        .maximum_slope_at_min_B = 4,
      },
    },
    .num_periodic_dir = 1,
    .periodic_dirs = {1},
    .num_species = 2,
    .species = {elc, ion},
    .field = field,
    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1], app_args.cuts[2] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

  double t_curr = 0.0, t_end = ctx.t_end; // Initial and final simulation times.
  int frame_curr = 0; // Initialize simulation.

  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
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
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, app_args.is_restart, false, -1.0);
  write_data(&trig_write_conf, &trig_write_phase, app, t_curr, app_args.is_restart, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false, t_curr > t_end, status.dt_actual);
    write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false, t_curr > t_end);

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
        calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false, true, status.dt_actual);
        write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  gkyl_gyrokinetic_app_stat_write(app);

  // Fetch simulation statistics.
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
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_gyrokinetic_app_print_timings(app, stdout);

freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
