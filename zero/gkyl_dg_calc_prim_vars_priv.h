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

typedef void (*euler_pkpm_prim_vars_t)(const double *bvar, 
  const double *vlasov_pkpm_moms, const double *statevec, 
  double* u_i, double* p_ij, double* T_ij, double* T_perp_over_m); 

typedef void (*euler_pkpm_source_t)(const double* qmem, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  double* out);

typedef void (*euler_pkpm_recovery_t)(const double *dxv, double nuHyp, 
  const double *bvarl, const double *bvarc, const double *bvarr, 
  const double *u_il, const double *u_ic, const double *u_ir, 
  const double *p_ijl, const double *p_ijc, const double *p_ijr, 
  const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr, 
  const double *statevecl, const double *statevecc, const double *statevecr, 
  const double *T_perp_over_m, const double *nu, const double *nu_vthsq, 
  double* div_b, double* bb_grad_u, double* div_p, double* p_force, double* p_perp_source, double* p_perp_div_b); 

typedef void (*pkpm_dist_mirror_force_t)(const double* T_perp_over_m, const double* f, const double* F_k_p_1, 
  double* g_dist_source, double* F_k_m_1); 

// for use in kernel tables
typedef struct { euler_pressure_t kernels[3]; } gkyl_dg_euler_pressure_kern_list;
typedef struct { euler_pkpm_prim_vars_t kernels[3]; } gkyl_dg_euler_pkpm_prim_vars_kern_list;
typedef struct { euler_pkpm_source_t kernels[3]; } gkyl_dg_euler_pkpm_source_kern_list;
typedef struct { euler_pkpm_recovery_t kernels[3]; } gkyl_dg_euler_pkpm_recovery_kern_list;
typedef struct { pkpm_dist_mirror_force_t kernels[3]; } gkyl_dg_pkpm_dist_mirror_force_kern_list;


// Pressure from state variables kernel list
GKYL_CU_D
static const gkyl_dg_euler_pressure_kern_list ser_pressure_kernels[] = {
  { NULL, euler_pressure_1x_ser_p1, euler_pressure_1x_ser_p2 }, // 0
  { NULL, euler_pressure_2x_ser_p1, NULL }, // 1
  { NULL, euler_pressure_3x_ser_p1, euler_pressure_3x_ser_p2 }, // 2
};

// PKPM primitive variables (u, p_ij, T_ij, T_perp/m)
GKYL_CU_D
static const gkyl_dg_euler_pkpm_prim_vars_kern_list ser_pkpm_prim_vars_kernels[] = {
  { NULL, euler_pkpm_prim_vars_1x_ser_p1, euler_pkpm_prim_vars_1x_ser_p2 }, // 0
  { NULL, euler_pkpm_prim_vars_2x_ser_p1, NULL }, // 1
  { NULL, euler_pkpm_prim_vars_3x_ser_p1, NULL }, // 2
};

// PKPM explicit source solve 
GKYL_CU_D
static const gkyl_dg_euler_pkpm_source_kern_list ser_pkpm_source_kernels[] = {
  { NULL, euler_pkpm_source_1x_ser_p1, euler_pkpm_source_1x_ser_p2 }, // 0
  { NULL, euler_pkpm_source_2x_ser_p1, NULL }, // 1
  { NULL, euler_pkpm_source_3x_ser_p1, NULL }, // 2
};

// PKPM recovery (in x) kernels
GKYL_CU_D
static const gkyl_dg_euler_pkpm_recovery_kern_list ser_pkpm_recovery_x_kernels[] = {
  { NULL, euler_pkpm_recovery_x_1x_ser_p1, euler_pkpm_recovery_x_1x_ser_p2 }, // 0
  { NULL, euler_pkpm_recovery_x_2x_ser_p1, NULL }, // 1
  { NULL, euler_pkpm_recovery_x_3x_ser_p1, NULL }, // 2
};

// PKPM recovery (in y) kernels
GKYL_CU_D
static const gkyl_dg_euler_pkpm_recovery_kern_list ser_pkpm_recovery_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, euler_pkpm_recovery_y_2x_ser_p1, NULL }, // 1
  { NULL, euler_pkpm_recovery_y_3x_ser_p1, NULL }, // 2
};

// PKPM recovery (in z) kernels
GKYL_CU_D
static const gkyl_dg_euler_pkpm_recovery_kern_list ser_pkpm_recovery_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, euler_pkpm_recovery_z_3x_ser_p1, NULL }, // 2
};

GKYL_CU_D
static const gkyl_dg_pkpm_dist_mirror_force_kern_list ser_pkpm_dist_mirror_force_kernels[] = {
  { NULL, pkpm_dist_mirror_force_1x1v_ser_p1, pkpm_dist_mirror_force_1x1v_ser_p2 }, // 0
  { NULL, pkpm_dist_mirror_force_2x1v_ser_p1, NULL }, // 1
  { NULL, pkpm_dist_mirror_force_3x1v_ser_p1, NULL }, // 2
};

GKYL_CU_D
static euler_pressure_t
choose_ser_euler_pressure_kern(int cdim, int poly_order)
{
  return ser_pressure_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static euler_pkpm_prim_vars_t
choose_ser_euler_pkpm_prim_vars_kern(int cdim, int poly_order)
{
  return ser_pkpm_prim_vars_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static euler_pkpm_source_t
choose_ser_euler_pkpm_source_kern(int cdim, int poly_order)
{
  return ser_pkpm_source_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static euler_pkpm_recovery_t
choose_ser_euler_pkpm_recovery_kern(int dir, int cdim, int poly_order)
{
  if (dir == 0)
    return ser_pkpm_recovery_x_kernels[cdim-1].kernels[poly_order];
  else if (dir == 1)
    return ser_pkpm_recovery_y_kernels[cdim-1].kernels[poly_order];
  else if (dir == 2)
    return ser_pkpm_recovery_z_kernels[cdim-1].kernels[poly_order];
  else 
    return NULL;
}

GKYL_CU_D
static pkpm_dist_mirror_force_t
choose_ser_pkpm_dist_mirror_force_kern(int cdim, int poly_order)
{
  return ser_pkpm_dist_mirror_force_kernels[cdim-1].kernels[poly_order];
}