// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_euler_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*euler_pressure_t)(const double gas_gamma, 
  const double *uvar, const double *statevec, 
  double* out);

// for use in kernel tables
typedef struct { euler_pressure_t kernels[3]; } gkyl_dg_euler_pressure_kern_list;

// Pressure from state variables kernel list
GKYL_CU_D
static const gkyl_dg_euler_pressure_kern_list ser_pressure_kernels[] = {
  { NULL, euler_pressure_1x_ser_p1, euler_pressure_1x_ser_p2 }, // 0
  { NULL, euler_pressure_2x_ser_p1, NULL }, // 1
  { NULL, euler_pressure_3x_ser_p1, NULL }, // 2
};

GKYL_CU_D
static euler_pressure_t
choose_ser_euler_pressure_kern(int cdim, int poly_order)
{
  return ser_pressure_kernels[cdim-1].kernels[poly_order];
}
