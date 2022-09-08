// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_euler_kernels.h>
#include <gkyl_mat.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*euler_pressure_t)(const double gas_gamma, 
  const double *uvar, const double *statevec, 
  double* GKYL_RESTRICT out);

typedef void (*euler_pkpm_pressure_t)(const double *u_i, const double *bvar, 
  const double *vlasov_pkpm_moms, const double *statevec, 
  double* GKYL_RESTRICT p_ij);

// for use in kernel tables
typedef struct { euler_pressure_t kernels[3]; } gkyl_dg_euler_pressure_kern_list;
typedef struct { euler_pkpm_pressure_t kernels[3]; } gkyl_dg_euler_pkpm_pressure_kern_list;

// Pressure from state variables kernel list
GKYL_CU_D
static const gkyl_dg_euler_pressure_kern_list ser_pressure_kernels[] = {
  { NULL, euler_pressure_1x_ser_p1, euler_pressure_1x_ser_p2 }, // 0
  { NULL, euler_pressure_2x_ser_p1, euler_pressure_2x_ser_p2 }, // 1
  { NULL, euler_pressure_3x_ser_p1, euler_pressure_3x_ser_p2 }, // 2
};

// PKPM pressure from state variables kernel list
GKYL_CU_D
static const gkyl_dg_euler_pkpm_pressure_kern_list ser_pkpm_pressure_kernels[] = {
  { NULL, euler_pkpm_pressure_1x_ser_p1, euler_pkpm_pressure_1x_ser_p2 }, // 0
  { NULL, euler_pkpm_pressure_2x_ser_p1, euler_pkpm_pressure_2x_ser_p2 }, // 1
  { NULL, euler_pkpm_pressure_3x_ser_p1, euler_pkpm_pressure_3x_ser_p2 }, // 2
};

GKYL_CU_D
static euler_pressure_t
choose_ser_euler_pressure_kern(int cdim, int poly_order)
{
  return ser_pressure_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static euler_pkpm_pressure_t
choose_ser_euler_pkpm_pressure_kern(int cdim, int poly_order)
{
  return ser_pkpm_pressure_kernels[cdim-1].kernels[poly_order];
}
