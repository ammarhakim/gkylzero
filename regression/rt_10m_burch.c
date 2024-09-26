// Burch et al. magnetic reconnection test for the 10-moment equations.
// Input parameters match the initial conditions in Section 3, from the article:
// J. M. TenBarge, J. Ng, J. Juno, L. Wang, A. M. Hakim and A. Bhattacharjee (2019), "An Extended MHD Study of the 16 October 2015 MMS Diffusion Region Crossing",
// Journal of Geophysical Research: Space Physics, Volume 124 (11): 8474-8487.
// https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JA026731

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_moment.h>
#include <gkyl_util.h>
#include <gkyl_wv_ten_moment.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct burch_ctx
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
  double beta2; // Magnetosheath plasma beta.
  double Ti1_over_Ti2; // Magnetospheric ion temperature / magnetosheath ion temperature.
  double Te1_over_Te2; // Magnetospheric electron temperature / magnetosheath ion temperature.

  // Derived physical quantities (using normalized code units).
  double light_speed; // Speed of light.
  double B0; // Reference magnetic field strength.
  double omega_pi; // Ion plasma frequency.
  double di; // Ion skin depth.
  double omega0; // Reference frequency.
  double psi0; // Reference magnetic scalar potential.
  double guide1; // Magnetospheric guide field strength.
  double guide2; // Magnetosheath guide field strength.
  double b1; // Magnetospheric magnetic field strength.
  double b2; // Magnetosheath magnetic field strength.
  double n1; // Magnetospheric number density.
  double n2; // Magnetosheath number density.
  double Ti2; // Magnetosheath ion temperature.

  double Te1; // Magnetospheric electron temperature.
  double Ti1; // Magnetospheric ion temperature.

  double Te2; // Magnetosheath electron temperature (so that the system is in force balance).

  // Simulation parameters.
  int Nx; // Cell count (x-direction);
  int Ny; // Cell count (y-direction);
  double Lx; // Domain size (x-direction).
  double Ly; // Domain size (y-direction).
  double k0_elc; // Closure parameter for electrons.
  double k0_ion; // Closure parameter for ions.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct burch_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 1.0; // Proton mass.
  double charge_ion = 1.0; // Proton charge.
  double mass_elc = mass_ion / 100.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double vAe = 0.2; // Electron Alfven velocity.
  double n0 = 1.0; // Reference number density.
  double beta2 = 2.748; // Magnetosheath plasma beta.
  double Ti1_over_Ti2 = 7.73 / 1.374; // Magnetospheric ion temperature / magnetosheath ion temperature.
  double Te1_over_Te2 = 1.288 / 1.374; // Magnetospheric electron temperature / magnetosheath ion temperature.

  // Derived physical quantities (using normalized code units).
  double light_speed = 1.0 / sqrt(mu0 * epsilon0); // Speed of light.
  double B0 = vAe * sqrt(n0 * mass_elc); // Reference magnetic field strength.
  double omega_pi = sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_ion)); // Ion plasma frequency.
  double di = light_speed / omega_pi; // Ion skin depth.
  double omega0 = 1.0 * di; // Reference frequency.
  double psi0 = 0.1 * B0 * di; // Reference magnetic scalar potential.
  double guide1 = 0.099 * B0; // Magnetospheric guide field strength.
  double guide2 = guide1; // Magnetosheath guide field strength.
  double b1 = 1.696 * B0; // Magnetospheric magnetic field strength.
  double b2 = 1.00 * B0; // Magnetosheath magnetic field strength.
  double n1 = 0.062 * n0; // Magnetospheric number density.
  double n2 = 1.0 * n0; // Magnetosheath number density.
  double Ti2 = beta2 * (b2 * b2) / (2.0 * n2 * mu0); // Magnetosheath ion temperature.

  double Te1 = Ti2 * Te1_over_Te2; // Magnetospheric electron temperature.
  double Ti1 = Ti2 * Ti1_over_Ti2; // Magnetospheric ion temperature.

  double Te2 = (0.5 * (b1 * b1 - b2 * b2) + 0.5 * (guide1 * guide1 - guide2 * guide2)
    + n1 * (Ti1 + Te1) - n2 * Ti2) / n2; // Magnetosheath electron temperature (so that the system is in force balance).

  // Simulation parameters.
  int Nx = 256; // Cell count (x-direction).
  int Ny = 128; // Cell count (y-direction).
  double Lx = 40.96 * di; // Domain size (x-direction).
  double Ly = 20.48 * di; // Domain size (y-direction).
  double k0_elc = 1.0; // Closure parameter for electrons.
  double k0_ion = 0.1; // Closure parameter for ions.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 250.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct burch_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .vAe = vAe,
    .n0 = n0,
    .beta2 = beta2,
    .Ti1_over_Ti2 = Ti1_over_Ti2,
    .Te1_over_Te2 = Te1_over_Te2,
    .light_speed = light_speed,
    .B0 = B0,
    .omega_pi = omega_pi,
    .di = di,
    .omega0 = omega0,
    .psi0 = psi0,
    .guide1 = guide1,
    .guide2 = guide2,
    .b1 = b1,
    .b2 = b2,
    .n1 = n1,
    .n2 = n2,
    .Ti2 = Ti2,
    .Te1 = Te1,
    .Ti1 = Ti1,
    .Te2 = Te2,
    .Nx = Nx,
    .Ny = Ny,
    .Lx = Lx,
    .Ly = Ly,
    .k0_elc = k0_elc,
    .k0_ion = k0_ion,
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
  double x = xn[0], y = xn[1];
  struct burch_ctx *app = ctx;

  double pi = app->pi;

  double mass_elc = app->mass_elc;
  double charge_elc = app->charge_elc;

  double omega0 = app->omega0;
  double psi0 = app->psi0;
  double guide1 = app->guide1;
  double guide2 = app->guide2;
  double b1 = app->b1;
  double b2 = app->b2;
  double n1 = app->n1;
  double Ti2 = app->Ti2;

  double Te1 = app->Te1;
  double Ti1 = app->Ti1;

  double Te2 = app->Te2;

  double Lx = app->Lx;
  double Ly = app->Ly;

  double b1x = 0.5 * (b2 + b1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (b2 - b1); // Magnetospheric magnetic field (x-direction).
  double b1y = 0.0; // Magnetospheric magnetic field (y-direction).
  double b1z = 0.5 * (guide2 - guide1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (guide2 + guide1); // Magnetospheric magnetic field (z-direction).

  double Ti_tot = 0.5 * (Ti2 - Ti1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Ti2 + Ti1); // Total ion temperature.
  double Te_tot = 0.5 * (Te2 - Te1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Te2 + Te1); // Total electron temperature.
  double n = (0.5 * (b1 * b1 - b1x * b1x) + 0.5 * (guide1 * guide1 - b1z * b1z) + n1 * (Ti1 + Te1)) / (Ti_tot + Te_tot); // Total number density.

  double Bx = b1x - psi0 * 4.0 * pi / Ly * sin(2.0 * pi * x / Lx) * sin(4.0 * pi * y / Ly); // Total magnetic field (x-direction).
  double By = b1y + psi0 * 2.0 * pi / Lx * cos(2.0 * pi * x / Lx) * (1.0 - cos(4.0 * pi * y / Ly)); // Total magnetic field (y-direction).
  double Bz = b1z; // Total magnetic field (z-direction).

  double Te_frac = Te_tot / (Te_tot + Ti_tot); // Fraction of total temperature from electrons.
  double Ti_frac = Ti_tot / (Te_tot + Ti_tot); // Fraction of total temperature from ions;

  double Jx = 0.5 * (guide2 - guide1) / omega0 * ((1.0 / cosh((y - Ly * 0.25) / omega0)) * (1.0 / cosh((y - Ly * 0.25) / omega0)) 
    - (1.0 / cosh((y - Ly * 0.75) / omega0)) * (1.0 / cosh((y - Ly * 0.75) / omega0)) 
    + (1.0 / cosh((y - Ly * 1.25) / omega0)) * (1.0 / cosh((y - Ly * 1.25) / omega0)) 
    - (1.0 / cosh((y + Ly * 0.25) / omega0)) * (1.0 / cosh((y + Ly * 0.25) / omega0))); // Total current density (x-direction).
  double Jy = 0.0; // Total current density (y-direction).
  double Jz  = -0.5 * (b2 + b1) / omega0 * ((1.0 / cosh((y - Ly * 0.25) / omega0)) * (1.0 / cosh((y - Ly * 0.25) / omega0)) 
    - (1.0 / cosh((y - Ly * 0.75) / omega0)) * (1.0 / cosh((y - Ly * 0.75) / omega0)) 
    + (1.0 / cosh((y - Ly * 1.25) / omega0)) * (1.0 / cosh((y - Ly * 1.25) / omega0)) 
    - (1.0 / cosh((y + Ly * 0.25) / omega0)) * (1.0 / cosh((y + Ly * 0.25) / omega0))) 
    - psi0 * sin(2.0 * pi * x / Lx) * ((2.0 * pi / Lx) * (2.0 * pi / Lx) * (1.0 - cos(4.0 * pi * y / Ly)) +
      (4.0 * pi / Ly) * (4.0 * pi / Ly) * cos(4.0 * pi * y / Ly)); // Total current density (z-direction).
   
  double Jxe = Jx * Te_frac; // Electron current density (x-direction).
  double Jye = Jy * Te_frac; // Electron current density (y-direction).
  double Jze = Jz * Te_frac; // Electron current density (z-direction).

  double rhoe = mass_elc * n; // Electron mass density.
  double momxe = (mass_elc / charge_elc) * Jxe; // Electron momentum density (x-direction).
  double momye = (mass_elc / charge_elc) * Jye; // Electron momentum density (y-direction).
  double momze = (mass_elc / charge_elc) * Jze; // Electron momentum density (z-direction).
  double pxxe = n * Te_tot + momxe * momxe / rhoe; // Electron pressure tensor (x-x component).
  double pxye = momxe * momye / rhoe; // Electron pressure tensor (x-y/y-x component).
  double pxze = momxe * momze / rhoe; // Electron pressure tensor (x-z/z-x component).
  double pyye = n*Te_tot + momye * momye / rhoe; // Electron pressure tensor (y-y component).
  double pyze = momye * momye / rhoe; // Electron pressure tensor (y-z/z-y component).
  double pzze = n*Te_tot + momze * momze / rhoe; // Electron pressure tensor (z-z component).

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
  double x = xn[0], y = xn[1];
  struct burch_ctx *app = ctx;

  double pi = app->pi;

  double mass_ion = app->mass_ion;
  double charge_ion = app->charge_ion;

  double omega0 = app->omega0;
  double psi0 = app->psi0;
  double guide1 = app->guide1;
  double guide2 = app->guide2;
  double b1 = app->b1;
  double b2 = app->b2;
  double n1 = app->n1;
  double Ti2 = app->Ti2;

  double Te1 = app->Te1;
  double Ti1 = app->Ti1;

  double Te2 = app->Te2;

  double Lx = app->Lx;
  double Ly = app->Ly;

  double b1x = 0.5 * (b2 + b1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (b2 - b1); // Magnetospheric magnetic field (x-direction).
  double b1y = 0.0; // Magnetospheric magnetic field (y-direction).
  double b1z = 0.5 * (guide2 - guide1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (guide2 + guide1); // Magnetospheric magnetic field (z-direction).

  double Ti_tot = 0.5 * (Ti2 - Ti1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Ti2 + Ti1); // Total ion temperature.
  double Te_tot = 0.5 * (Te2 - Te1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Te2 + Te1); // Total electron temperature.
  double n = (0.5 * (b1 * b1 - b1x * b1x) + 0.5 * (guide1 * guide1 - b1z * b1z) + n1 * (Ti1 + Te1)) / (Ti_tot + Te_tot); // Total number density.

  double Bx = b1x - psi0 * 4.0 * pi / Ly * sin(2.0 * pi * x / Lx) * sin(4.0 * pi * y / Ly); // Total magnetic field (x-direction).
  double By = b1y + psi0 * 2.0 * pi / Lx * cos(2.0 * pi * x / Lx) * (1.0 - cos(4.0 * pi * y / Ly)); // Total magnetic field (y-direction).
  double Bz = b1z; // Total magnetic field (z-direction).

  double Te_frac = Te_tot / (Te_tot + Ti_tot); // Fraction of total temperature from electrons.
  double Ti_frac = Ti_tot / (Te_tot + Ti_tot); // Fraction of total temperature from ions;

  double Jx = 0.5 * (guide2 - guide1) / omega0 * ((1.0 / cosh((y - Ly * 0.25) / omega0)) * (1.0 / cosh((y - Ly * 0.25) / omega0)) 
    - (1.0 / cosh((y - Ly * 0.75) / omega0)) * (1.0 / cosh((y - Ly * 0.75) / omega0)) 
    + (1.0 / cosh((y - Ly * 1.25) / omega0)) * (1.0 / cosh((y - Ly * 1.25) / omega0)) 
    - (1.0 / cosh((y + Ly * 0.25) / omega0)) * (1.0 / cosh((y + Ly * 0.25) / omega0))); // Total current density (x-direction).
  double Jy = 0.0; // Total current density (y-direction).
  double Jz  = -0.5 * (b2 + b1) / omega0 * ((1.0 / cosh((y - Ly * 0.25) / omega0)) * (1.0 / cosh((y - Ly * 0.25) / omega0)) 
    - (1.0 / cosh((y - Ly * 0.75) / omega0)) * (1.0 / cosh((y - Ly * 0.75) / omega0)) 
    + (1.0 / cosh((y - Ly * 1.25) / omega0)) * (1.0 / cosh((y - Ly * 1.25) / omega0)) 
    - (1.0 / cosh((y + Ly * 0.25) / omega0)) * (1.0 / cosh((y + Ly * 0.25) / omega0))) 
    - psi0 * sin(2.0 * pi * x / Lx) * ((2.0 * pi / Lx) * (2.0 * pi / Lx) * (1.0 - cos(4.0 * pi * y / Ly)) +
      (4.0 * pi / Ly) * (4.0 * pi / Ly) * cos(4.0 * pi * y / Ly)); // Total current density (z-direction).
   
  double Jxi = Jx * Ti_frac; // Ion current density (x-direction).
  double Jyi = Jy * Ti_frac; // Ion current density (y-direction).
  double Jzi = Jz * Ti_frac; // Ion current density (z-direction).

  double rhoi = mass_ion * n; // Ion mass density.
  double momxi = (mass_ion / charge_ion) * Jxi; // Ion momentum density (x-direction).
  double momyi = (mass_ion / charge_ion) * Jyi; // Ion momentum density (y-direction).
  double momzi = (mass_ion / charge_ion) * Jzi; // Ion momentum density (z-direction).
  double pxxi = n * Ti_tot + momxi * momxi / rhoi; // Ion pressure tensor (x-x component).
  double pxyi = momxi * momyi / rhoi; // Ion pressure tensor (x-y/y-x component).
  double pxzi = momxi * momzi / rhoi; // Ion pressure tensor (x-z/z-x component).
  double pyyi = n * Ti_tot + momyi * momyi / rhoi; // Ion pressure tensor (y-y component).
  double pyzi = momyi * momyi / rhoi; // Ion pressure tensor (y-z/z-y component).
  double pzzi = n * Ti_tot + momzi * momzi / rhoi; // Ion pressure tensor (z-z component).

  // Set ion mass density,
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
  double x = xn[0], y = xn[1];
  struct burch_ctx *app = ctx;

  double pi = app->pi;

  double omega0 = app->omega0;
  double psi0 = app->psi0;
  double guide1 = app->guide1;
  double guide2 = app->guide2;
  double b1 = app->b1;
  double b2 = app->b2;
  double n1 = app->n1;
  double Ti2 = app->Ti2;

  double Te1 = app->Te1;
  double Ti1 = app->Ti1;

  double Te2 = app->Te2;

  double Lx = app->Lx;
  double Ly = app->Ly;

  double b1x = 0.5 * (b2 + b1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (b2 - b1); // Magnetospheric magnetic field (x-direction).
  double b1y = 0.0; // Magnetospheric magnetic field (y-direction).
  double b1z = 0.5 * (guide2 - guide1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (guide2 + guide1); // Magnetospheric magnetic field (z-direction).

  double Ti_tot = 0.5 * (Ti2 - Ti1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Ti2 + Ti1); // Total ion temperature.
  double Te_tot = 0.5 * (Te2 - Te1) * (tanh((y - Ly * 0.25) / omega0) - tanh((y - Ly * 0.75) / omega0)
    + tanh((y - Ly * 1.25) / omega0) - tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Te2 + Te1); // Total electron temperature.
  double n = (0.5 * (b1 * b1 - b1x * b1x) + 0.5 * (guide1 * guide1 - b1z * b1z) + n1 * (Ti1 + Te1)) / (Ti_tot + Te_tot); // Total number density.

  double Bx = b1x - psi0 * 4.0 * pi / Ly * sin(2.0 * pi * x / Lx) * sin(4.0 * pi * y / Ly); // Total magnetic field (x-direction).
  double By = b1y + psi0 * 2.0 * pi / Lx * cos(2.0 * pi * x / Lx) * (1.0 - cos(4.0 * pi * y / Ly)); // Total magnetic field (y-direction).
  double Bz = b1z; // Total magnetic field (z-direction).

  // Set electric field.
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // Set magnetic field.
  fout[3] = Bx, fout[4] = By; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_moment_app_write(app, t_curr, frame);
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

  struct burch_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);

  // Electron/ion equations.
  struct gkyl_wv_eqn *elc_ten_moment = gkyl_wv_ten_moment_new(ctx.k0_elc, false, app_args.use_gpu);
  struct gkyl_wv_eqn *ion_ten_moment = gkyl_wv_ten_moment_new(ctx.k0_ion, false, app_args.use_gpu);

  struct gkyl_moment_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .equation = elc_ten_moment,
    .evolve = true,
    .init = evalElcInit,
    .ctx = &ctx,
  };

  struct gkyl_moment_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .equation = ion_ten_moment,
    .evolve = true,
    .init = evalIonInit,
    .ctx = &ctx,
  };

  // Field.
  struct gkyl_moment_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .mag_error_speed_fact = 1.0,
    
    .evolve = true,
    .init = evalFieldInit,
    .ctx = &ctx,
  };

  int nrank = 1; // Number of processes in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  // Create global range.
  int cells[] = { NX, NY };
  int dim = sizeof(cells) / sizeof(cells[0]);

  int cuts[dim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < dim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < dim; d++) {
    cuts[d] = 1;
  }
#endif

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < dim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Moment app.
  struct gkyl_moment app_inp = {
    .name = "10m_burch",

    .ndim = 2,
    .lower = { 0.0, 0.0 },
    .upper = { ctx.Lx / ctx.di, ctx.Ly / ctx.di }, 
    .cells = { NX, NY },

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },
    .cfl_frac = ctx.cfl_frac,

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
  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // Initialize simulation.
  gkyl_moment_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr, false);

  // Compute estimate of maximum stable time-step.
  double dt = gkyl_moment_app_max_dt(app);

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_moment_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_moment_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_moment_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_moment_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_moment_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  write_data(&io_trig, app, t_curr, false);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  // Free resources after simulation completion.
  gkyl_wv_eqn_release(elc_ten_moment);
  gkyl_wv_eqn_release(ion_ten_moment);
  gkyl_comm_release(comm);
  gkyl_moment_app_release(app);  
  
mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
