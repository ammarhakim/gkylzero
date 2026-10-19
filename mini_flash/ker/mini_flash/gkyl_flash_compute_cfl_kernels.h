#pragma once

#include <gkyl_util.h>

double flash_compute_cfl_vlasov_pkpm_boundary_accel_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *div_b, const double *pkpm_accel_vars, const int edge);

double flash_compute_cfl_vlasov_pkpm_surf_accel_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *pkpm_accel_vars);

double flash_compute_cfl_vlasov_pkpm_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, const double *pkpm_accel_vars);

double flash_compute_cfl_vlasov_pkpm_vol_1x1v_ser_p1(const double *w, const double *dxv, const double* pkpm_accel_vars);

double flash_compute_cfl_lbo_pkpm_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double* nuPrimMomsSum);