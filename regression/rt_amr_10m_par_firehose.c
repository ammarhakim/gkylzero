#include <gkyl_amr_core.h>

struct amr_10m_par_firehose_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double vAe; // Electron Alfven velocity.
  double n0; // Reference number density.

  // Derived physical quantities (using normalized code units).
  double light_speed; // Speed of light.
  double B0; // Reference magnetic field strength.
  double beta; // Trace proton plasma beta.
  double dbeta; // Parallel proton plasma beta - perpendicular proton plasma beta.
  double beta_par; // Parallel proton plasma beta.
  double beta_perp; // Perpendicular proton plasma beta.

  double vte; // Electron thermal velocity.
  double Te; // Electron temperature.
  
  double Ti_par; // Parallel ion temperature.
  double Ti_perp; // Perpendicular ion temperature.

  double omega_ci; // Ion cyclotron frequency.
  double omega_pe; // Electron plasma frequency.
  double de; // Electron skin depth.
  double omega_pi; // Ion plasma frequency.
  double di; // Ion skin depth.
  double lambdaD; // Electron Debye length.

  double noise_amp; // Noise level for perturbation.
  int mode_init; // Initial wave mode to perturb with noise.
  int mode_final; // Final wave mode to perturb with noise.

  double k0_elc; // Electron closure parameter.
  double k0_ion; // Ion closure parameter.

  // Simulation parameters.
  int Nx; // Coarse cell count (x-direction).
  int ref_factor; // Refinement factor.
  double Lx; // Coarse domain size (x-direction).
  double fine_Lx; // Fine domain size (x-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct amr_10m_par_firehose_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 1836.0; // Proton mass.
  double charge_ion = 1.0; // Proton charge.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.

  double vAe = 0.0125; // Electron Alfven velocity.
  double n0 = 1.0; // Reference number density.

  // Derived physical quantities (using normalized code units).
  double light_speed = 1.0 / sqrt(mu0 * epsilon0); // Speed of light.
  double B0 = vAe * sqrt(mu0 * n0 * mass_elc); // Reference magnetic field strength.
  double beta = 300.0 / pi; // Trace proton plasma beta.
  double dbeta = 100.0; // Parallel proton plasma beta - perpendicular proton plasma beta.
  double beta_par = beta + 2.0 * dbeta / 3.0; // Parallel proton plasma beta.
  double beta_perp = beta - dbeta / 3.0; // Perpendicular proton plasma beta.

  double vte = vAe * sqrt(beta); // Electron thermal velocity.
  double Te = vte * vte * mass_elc / 2.0; // Electron temperature.
  
  double Ti_par = vAe * vAe * (beta_par * mass_elc / 2.0); // Parallel ion temperature.
  double Ti_perp = vAe * vAe * (beta_perp * mass_elc / 2.0); // Perpendicular ion temperature.

  double omega_ci = charge_ion * B0 / mass_ion; // Ion cyclotron frequency.
  double omega_pe = sqrt(n0 * charge_elc * charge_elc / (epsilon0 * mass_elc)); // Electron plasma frequency.
  double de = light_speed / omega_pe; // Electron skin depth.
  double omega_pi = sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_ion)); // Ion plasma frequency.
  double di = light_speed / omega_pi; // Ion skin depth.
  double lambdaD = vte / omega_pe; // Electron Debye length.

  double noise_amp = 1.0e-6 * B0; // Noise level for perturbation.
  int mode_init = 1; // Initial wave mode to perturb with noise.
  int mode_final = 48; // Final wave mode to perturb with noise.

  double k0_elc = 0.1 / de; // Electron closure parameter.
  double k0_ion = 0.1 / di; // Ion closure parameter.

  // Simulation parameters.
  int Nx = 64; // Coarse cell count (x-direction).
  int ref_factor = 2; // Refinement factor.
  double Lx = 300.0 * di; // Coarse domain size (x-direction).
  double fine_Lx = 0.5 * (300.0 * di); // Fine domain size (x-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 10.0 / omega_ci; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct amr_10m_par_firehose_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .vAe = vAe,
    .n0 = n0,
    .light_speed = light_speed,
    .B0 = B0,
    .beta = beta,
    .dbeta = dbeta,
    .beta_par = beta_par,
    .beta_perp = beta_perp,
    .vte = vte,
    .Te = Te,
    .Ti_par = Ti_par,
    .Ti_perp = Ti_perp,
    .omega_ci = omega_ci,
    .omega_pe = omega_pe,
    .de = de,
    .omega_pi = omega_pi,
    .di = di,
    .lambdaD = lambdaD,
    .noise_amp = noise_amp,
    .mode_init = mode_init,
    .mode_final = mode_final,
    .k0_elc = k0_elc,
    .k0_ion = k0_ion,
    .Nx = Nx,
    .ref_factor = ref_factor,
    .Lx = Lx,
    .fine_Lx = fine_Lx,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct amr_10m_par_firehose_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_10m_par_firehose_ctx *app = &new_ctx;

  double mass_elc = app->mass_elc;

  double n0 = app->n0;

  double Te = app->Te;

  double rhoe = mass_elc * n0; // Electron mass density.
  double momxe = 0.0; // Electron momentum density (x-direction).
  double momye = 0.0; // Electron momentum density (y-direction).
  double momze = 0.0; // Electron momentum density (z-direction).
  double pxxe = n0 * Te + momxe * momxe / rhoe; // Electron pressure tensor (x-x component).
  double pxye = momxe * momye / rhoe; // Electron pressure tensor (x-y/y-x component).
  double pxze = momxe * momze / rhoe; // Electron pressure tensor (x-z/z-x component).
  double pyye = n0 * Te + momye * momye / rhoe; // Electron pressure tensor (y-y component).
  double pyze = momye * momye / rhoe; // Electron pressure tensor (y-z/z-y component).
  double pzze = n0 * Te + momze * momze / rhoe; // Electron pressure tensor (z-z component).

  // Set electron mass density.
  fout[0] = rhoe;
  // Set electron momentum density.
  fout[1] = momxe; fout[2] = momye; fout[3] = momze;
  // Set electron pressure tensor.
  fout[4] = pxxe; fout[5] = pxye; fout[6] = pxze;  
  fout[7] = pyye; fout[8] = pyze; fout[9] = pzze;
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct amr_10m_par_firehose_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_10m_par_firehose_ctx *app = &new_ctx;

  double mass_ion = app->mass_ion;

  double n0 = app->n0;

  double Ti_par = app->Ti_par;
  double Ti_perp = app->Ti_perp;

  double rhoi = mass_ion * n0; // Ion mass density.
  double momxi = 0.0; // Ion momentum density (x-direction).
  double momyi = 0.0; // Ion momentum density (y-direction).
  double momzi = 0.0; // Ion momentum density (z-direction).
  double pxxi = n0 * Ti_par + momxi * momxi / rhoi; // Ion pressure tensor (x-x component).
  double pxyi = momxi * momyi / rhoi; // Ion pressure tensor (x-y/y-x component).
  double pxzi = momxi * momzi / rhoi; // Ion pressure tensor (x-z/z-x component).
  double pyyi = n0 * Ti_perp + momyi * momyi / rhoi; // Ion pressure tensor (y-y component).
  double pyzi = momyi * momyi / rhoi; // Ion pressure tensor (y-z/z-y component).
  double pzzi = n0 * Ti_perp + momzi * momzi / rhoi; // Ion pressure tensor (z-z component).

  // Set ion mass density.
  fout[0] = rhoi;
  // Set ion momentum density.
  fout[1] = momxi; fout[2] = momyi; fout[3] = momzi;
  // Set ion pressure tensor.
  fout[4] = pxxi; fout[5] = pxyi; fout[6] = pxzi;  
  fout[7] = pyyi; fout[8] = pyzi; fout[9] = pzzi;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct amr_10m_par_firehose_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_10m_par_firehose_ctx *app = &new_ctx;

  double pi = app->pi;

  double B0 = app->B0;

  double noise_amp = app->noise_amp;
  double mode_init = app->mode_init;
  double mode_final = app->mode_final;

  double Lx = app->Lx;

  double Bx = B0; // Total magnetic field (x-direction).
  double By = 0.0;
  double Bz = 0.0;

  double alpha = noise_amp * Bx; // Applied amplitude.
  double kx = 2.0 * pi / Lx; // Wave number (x-direction).

  pcg64_random_t rng = gkyl_pcg64_init(0); // Random number generator.

  for (int i = mode_init; i < mode_final; i++)
  {
    By -= alpha * gkyl_pcg64_rand_double(&rng) * sin(i * kx * x + 2.0 * pi * gkyl_pcg64_rand_double(&rng)); // Total magnetic field (y-direction).
    Bz -= alpha * gkyl_pcg64_rand_double(&rng) * sin(i * kx * x + 2.0 * pi * gkyl_pcg64_rand_double(&rng)); // Total magnetic field (z-direction).
  }

  // Set electric field.
  fout[0] = 0.0; fout[1] = 0.0; fout[2] = 0.0;
  // Set magnetic field.
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

int main(int argc, char **argv)
{
  struct amr_10m_par_firehose_ctx ctx = create_ctx(); // Context for initialization functions.

  struct ten_moment_1d_single_init init = {
    .base_Nx = ctx.Nx,
    .ref_factor = ctx.ref_factor,

    .coarse_x1 = 0.0,
    .coarse_x2 = ctx.Lx,

    .refined_x1 = (0.5 * ctx.Lx) - (0.5 * ctx.fine_Lx),
    .refined_x2 = (0.5 * ctx.Lx) + (0.5 * ctx.fine_Lx),

    .eval_elc = evalElcInit,
    .eval_ion = evalIonInit,
    .eval_field = evalFieldInit,

    .k0_elc = ctx.k0_elc,
    .k0_ion = ctx.k0_ion,

    .light_speed = 1.0,
    .e_fact = 0.0,
    .b_fact = 1.0,

    .epsilon0 = ctx.epsilon0,
    .mass_elc = ctx.mass_elc,
    .charge_elc = ctx.charge_elc,
    .mass_ion = ctx.mass_ion,
    .charge_ion = ctx.charge_ion,

    .ten_moment_output = "amr_10m_par_firehose",

    .low_order_flux = false,
    .cfl_frac = ctx.cfl_frac,

    .t_end = ctx.t_end,
    .num_frames = ctx.num_frames,
    .dt_failure_tol = ctx.dt_failure_tol,
    .num_failures_max = ctx.num_failures_max,
  };

  ten_moment_1d_run_single(argc, argv, &init);
}