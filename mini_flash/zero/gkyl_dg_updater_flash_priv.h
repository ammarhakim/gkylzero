#pragma once

#include <gkyl_dg_eqn.h>
#include <gkyl_flash_compute_cfl_kernels.h>
#include <gkyl_hyper_dg.h>

// Function types for CFL frequency calculation kernels
typedef double (*flash_compute_cfl_vlasov_pkpm_boundary_accel_t)(const double *w, const double *dxv, 
  const double *div_b, const double *pkpm_accel_vars, const int edge);

typedef double (*flash_compute_cfl_vlasov_pkpm_surf_accel_t)(const double *w, const double *dxv, 
  const double *pkpm_accel_vars);

typedef double (*flash_compute_cfl_vlasov_pkpm_surf_t)(const double *w, const double *dxv, 
  const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, const double *pkpm_accel_vars);

typedef double (*flash_compute_cfl_vlasov_pkpm_vol_t)(const double *w, const double *dxv, const double* pkpm_accel_vars);

typedef double (*flash_compute_cfl_lbo_pkpm_vol_t)(const double *w, const double *dxv, const double *nuSum, const double* nuPrimMomsSum);

// Object containing kernels for timestep calculation
struct gkyl_flash_compute_cfl_kernels {
  flash_compute_cfl_vlasov_pkpm_boundary_accel_t boundary_accel;
  flash_compute_cfl_vlasov_pkpm_surf_accel_t surf_accel;
  flash_compute_cfl_vlasov_pkpm_surf_t surf_stream;
  flash_compute_cfl_vlasov_pkpm_vol_t vol;
  flash_compute_cfl_lbo_pkpm_vol_t lbo_vol;
};

// Updater object
struct gkyl_dg_updater_flash {
  struct gkyl_dg_eqn *eqn_vlasov; // Equation object for flow-frame kinetic equation.
  struct gkyl_hyper_dg *up_vlasov; // Solver for flow-frame kinetic equation.

  struct gkyl_flash_compute_cfl_kernels compute_cfl_kernels;

  double vlasov_tm; // Total time spent in solving kinetic equation.
};

// Timestep calculation kernel lists
static const flash_compute_cfl_vlasov_pkpm_boundary_accel_t flash_compute_cfl_vlasov_pkpm_boundary_accel_kernels[3] = {
  flash_compute_cfl_vlasov_pkpm_boundary_accel_1x1v_ser_p1,
  NULL,
  NULL
};

static const flash_compute_cfl_vlasov_pkpm_surf_accel_t flash_compute_cfl_vlasov_pkpm_surf_accel_kernels[3] = {
  flash_compute_cfl_vlasov_pkpm_surf_accel_1x1v_ser_p1,
  NULL,
  NULL
};

static const flash_compute_cfl_vlasov_pkpm_surf_t flash_compute_cfl_vlasov_pkpm_surf_kernels[3] = {
  flash_compute_cfl_vlasov_pkpm_surfx_1x1v_ser_p1,
  NULL, 
  NULL
};

static const flash_compute_cfl_vlasov_pkpm_vol_t flash_compute_cfl_vlasov_pkpm_vol_kernels[3] = {
  flash_compute_cfl_vlasov_pkpm_vol_1x1v_ser_p1,
  NULL,
  NULL
};

static const flash_compute_cfl_lbo_pkpm_vol_t flash_compute_cfl_lbo_pkpm_vol_kernels[3] = {
  flash_compute_cfl_lbo_pkpm_vol_1x1v_ser_p1,
  NULL,
  NULL
};

// Kernel selector functions
static flash_compute_cfl_vlasov_pkpm_boundary_accel_t
choose_flash_compute_cfl_vlasov_pkpm_boundary_accel_kernel(int cdim)
{
  return flash_compute_cfl_vlasov_pkpm_boundary_accel_kernels[cdim-1];
}

static flash_compute_cfl_vlasov_pkpm_surf_accel_t
choose_flash_compute_cfl_vlasov_pkpm_surf_accel_kernel(int cdim)
{
  return flash_compute_cfl_vlasov_pkpm_surf_accel_kernels[cdim-1];
}

static flash_compute_cfl_vlasov_pkpm_surf_t
choose_flash_compute_cfl_vlasov_pkpm_surf_kernel(int cdim)
{
  return flash_compute_cfl_vlasov_pkpm_surf_kernels[cdim-1];
}

static flash_compute_cfl_vlasov_pkpm_vol_t
choose_flash_compute_cfl_vlasov_pkpm_vol_kernel(int cdim)
{
  return flash_compute_cfl_vlasov_pkpm_vol_kernels[cdim-1];
}

static flash_compute_cfl_lbo_pkpm_vol_t
choose_flash_compute_cfl_lbo_pkpm_vol_kernel(int cdim)
{
  return flash_compute_cfl_lbo_pkpm_vol_kernels[cdim-1];
}