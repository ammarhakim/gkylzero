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

typedef void (*euler_pkpm_pressure_t)(const double *bvar, 
  const double *vlasov_pkpm_moms, const double *statevec, 
  double* GKYL_RESTRICT p_ij);

typedef void (*euler_pkpm_pressure_source_t)(const double* qmem, 
  const double *nu, const double *nu_vth_sq, 
  const double *vlasov_pkpm_moms, const double *u_i, const double *euler_pkpm, 
  double* GKYL_RESTRICT out);

typedef void (*euler_pkpm_pressure_div_t)(const double *dxv, 
  const double *p_ijl, const double *p_ijc, const double *p_ijr, 
  double* GKYL_RESTRICT div_p);

typedef void (*euler_pkpm_bb_grad_u_t)(const double *dxv, 
  const double *bvar, const double *u_il, const double *u_ic, const double *u_ir, 
  double* GKYL_RESTRICT bb_grad_u);

typedef void (*euler_pkpm_p_force_t)(const double *bvar, const double *div_p, 
  const double *vlasov_pkpm_moms, 
  double* GKYL_RESTRICT p_force);

// for use in kernel tables
typedef struct { euler_pressure_t kernels[4]; } gkyl_dg_euler_pressure_kern_list;
typedef struct { euler_pkpm_pressure_t kernels[4]; } gkyl_dg_euler_pkpm_pressure_kern_list;
typedef struct { euler_pkpm_pressure_source_t kernels[4]; } gkyl_dg_euler_pkpm_pressure_source_kern_list;
typedef struct { euler_pkpm_pressure_div_t kernels[4]; } gkyl_dg_euler_pkpm_pressure_div_kern_list;
typedef struct { euler_pkpm_bb_grad_u_t kernels[4]; } gkyl_dg_euler_pkpm_bb_grad_u_kern_list;
typedef struct { euler_pkpm_p_force_t kernels[4]; } gkyl_dg_euler_pkpm_p_force_kern_list;


// Pressure from state variables kernel list
GKYL_CU_D
static const gkyl_dg_euler_pressure_kern_list ser_pressure_kernels[] = {
  { NULL, euler_pressure_1x_ser_p1, euler_pressure_1x_ser_p2, NULL }, // 0
  { NULL, euler_pressure_2x_ser_p1, euler_pressure_2x_ser_p2, NULL }, // 1
  { NULL, euler_pressure_3x_ser_p1, euler_pressure_3x_ser_p2, NULL }, // 2
};

// PKPM pressure from state variables kernel list
GKYL_CU_D
static const gkyl_dg_euler_pkpm_pressure_kern_list ser_pkpm_pressure_kernels[] = {
  { NULL, euler_pkpm_pressure_1x_ser_p1, euler_pkpm_pressure_1x_ser_p2, euler_pkpm_pressure_1x_ser_p3 }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
};

// PKPM pressure sourceation from state variables kernel list
GKYL_CU_D
static const gkyl_dg_euler_pkpm_pressure_source_kern_list ser_pkpm_pressure_source_kernels[] = {
  { NULL, euler_pkpm_pressure_source_1x_ser_p1, euler_pkpm_pressure_source_1x_ser_p2, euler_pkpm_pressure_source_1x_ser_p3 }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
};

// PKPM divergence of pressure
GKYL_CU_D
static const gkyl_dg_euler_pkpm_pressure_div_kern_list ser_pkpm_pressure_div_kernels[] = {
  { NULL, euler_pkpm_pressure_div_x_1x_ser_p1, euler_pkpm_pressure_div_x_1x_ser_p2, euler_pkpm_pressure_div_x_1x_ser_p3 }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
};

// PKPM bb : grad(u)
GKYL_CU_D
static const gkyl_dg_euler_pkpm_bb_grad_u_kern_list ser_pkpm_bb_grad_u_kernels[] = {
  { NULL, euler_pkpm_bb_grad_u_x_1x_ser_p1, euler_pkpm_bb_grad_u_x_1x_ser_p2, euler_pkpm_bb_grad_u_x_1x_ser_p3 }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
};

// PKPM total pressure force, p_force = 1/rho (b . div(P) + p_perp div(b))
GKYL_CU_D
static const gkyl_dg_euler_pkpm_p_force_kern_list ser_pkpm_p_force_kernels[] = {
  { NULL, euler_pkpm_p_force_1x_ser_p1, euler_pkpm_p_force_1x_ser_p2, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
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

GKYL_CU_D
static euler_pkpm_pressure_source_t
choose_ser_euler_pkpm_pressure_source_kern(int cdim, int poly_order)
{
  return ser_pkpm_pressure_source_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static euler_pkpm_pressure_div_t
choose_ser_euler_pkpm_pressure_div_kern(int cdim, int poly_order)
{
  return ser_pkpm_pressure_div_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static euler_pkpm_bb_grad_u_t
choose_ser_euler_pkpm_bb_grad_u_kern(int cdim, int poly_order)
{
  return ser_pkpm_bb_grad_u_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static euler_pkpm_p_force_t
choose_ser_euler_pkpm_p_force_kern(int cdim, int poly_order)
{
  return ser_pkpm_p_force_kernels[cdim-1].kernels[poly_order];
}