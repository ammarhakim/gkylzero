// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_euler_pkpm_kernels.h>
#include <gkyl_vlasov_pkpm_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*pkpm_prim_t)(const double *bvar, 
  const double *vlasov_pkpm_moms, const double *statevec, 
  double* u_i, double* p_ij, double* T_ij, 
  double* rho_inv, double* T_perp_over_m, double* T_perp_over_m_inv); 

typedef void (*pkpm_source_t)(const double* qmem, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  double* out);

typedef void (*pkpm_recovery_t)(const double *dxv, 
  const double *bvarl, const double *bvarc, const double *bvarr, 
  const double *u_il, const double *u_ic, const double *u_ir, 
  const double *p_ijl, const double *p_ijc, const double *p_ijr, 
  const double *pkpm_div_ppar, const double *rho_inv, const double *T_perp_over_m, const double *nu, 
  double* div_p, double* pkpm_accel_vars); 

typedef void (*pkpm_dist_mirror_force_t)(const double *w, const double *dxv, 
  const double* T_perp_over_m, const double* T_perp_over_m_inv, 
  const double *nu_vthsq, const double* pkpm_accel_vars, 
  const double* f, const double* F_k_p_1, 
  double* g_dist_source, double* F_k_m_1); 

typedef void (*pkpm_pressure_t)(const double *w, const double *dxv, 
     const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

// for use in kernel tables
typedef struct { pkpm_prim_t kernels[3]; } gkyl_dg_pkpm_prim_kern_list;
typedef struct { pkpm_source_t kernels[3]; } gkyl_dg_pkpm_source_kern_list;
typedef struct { pkpm_recovery_t kernels[3]; } gkyl_dg_pkpm_recovery_kern_list;
typedef struct { pkpm_dist_mirror_force_t kernels[3]; } gkyl_dg_pkpm_dist_mirror_force_kern_list;
typedef struct { pkpm_pressure_t kernels[3]; } gkyl_dg_pkpm_pressure_kern_list;

// PKPM primitive variables (u, p_ij, T_ij, T_perp/m, (T_perp/m)^-1)
GKYL_CU_D
static const gkyl_dg_pkpm_prim_kern_list ser_pkpm_prim_kernels[] = {
  { NULL, euler_pkpm_prim_vars_1x_ser_p1, euler_pkpm_prim_vars_1x_ser_p2 }, // 0
  { NULL, euler_pkpm_prim_vars_2x_ser_p1, NULL }, // 1
  { NULL, euler_pkpm_prim_vars_3x_ser_p1, NULL }, // 2
};

// PKPM explicit source solve 
GKYL_CU_D
static const gkyl_dg_pkpm_source_kern_list ser_pkpm_source_kernels[] = {
  { NULL, euler_pkpm_source_1x_ser_p1, euler_pkpm_source_1x_ser_p2 }, // 0
  { NULL, euler_pkpm_source_2x_ser_p1, NULL }, // 1
  { NULL, euler_pkpm_source_3x_ser_p1, NULL }, // 2
};

// PKPM recovery (in x) kernels
GKYL_CU_D
static const gkyl_dg_pkpm_recovery_kern_list ser_pkpm_recovery_x_kernels[] = {
  { NULL, euler_pkpm_recovery_x_1x_ser_p1, euler_pkpm_recovery_x_1x_ser_p2 }, // 0
  { NULL, euler_pkpm_recovery_x_2x_ser_p1, NULL }, // 1
  { NULL, euler_pkpm_recovery_x_3x_ser_p1, NULL }, // 2
};

// PKPM recovery (in y) kernels
GKYL_CU_D
static const gkyl_dg_pkpm_recovery_kern_list ser_pkpm_recovery_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, euler_pkpm_recovery_y_2x_ser_p1, NULL }, // 1
  { NULL, euler_pkpm_recovery_y_3x_ser_p1, NULL }, // 2
};

// PKPM recovery (in z) kernels
GKYL_CU_D
static const gkyl_dg_pkpm_recovery_kern_list ser_pkpm_recovery_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, euler_pkpm_recovery_z_3x_ser_p1, NULL }, // 2
};

// PKPM distribution function source in mirror force (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_pkpm_dist_mirror_force_kern_list ser_pkpm_dist_mirror_force_kernels[] = {
  { NULL, pkpm_dist_mirror_force_1x1v_ser_p1, pkpm_dist_mirror_force_1x1v_ser_p2 }, // 0
  { NULL, pkpm_dist_mirror_force_2x1v_ser_p1, NULL }, // 1
  { NULL, pkpm_dist_mirror_force_3x1v_ser_p1, NULL }, // 2
};

// PKPM distribution function source in mirror force (Tensor basis)
GKYL_CU_D
static const gkyl_dg_pkpm_dist_mirror_force_kern_list ten_pkpm_dist_mirror_force_kernels[] = {
  { NULL, pkpm_dist_mirror_force_1x1v_ser_p1, pkpm_dist_mirror_force_1x1v_tensor_p2 }, // 0
  { NULL, pkpm_dist_mirror_force_2x1v_ser_p1, NULL }, // 1
  { NULL, pkpm_dist_mirror_force_3x1v_ser_p1, NULL }, // 2
};

// PKPM pressure (in x) kernels (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_pkpm_pressure_kern_list ser_pkpm_pressure_x_kernels[] = {
  { NULL, vlasov_pkpm_pressure_x_1x1v_ser_p1, vlasov_pkpm_pressure_x_1x1v_ser_p2 }, // 0
  { NULL, vlasov_pkpm_pressure_x_2x1v_ser_p1, NULL }, // 1
  { NULL, vlasov_pkpm_pressure_x_3x1v_ser_p1, NULL }, // 2
};

// PKPM pressure (in y) kernels (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_pkpm_pressure_kern_list ser_pkpm_pressure_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_pkpm_pressure_y_2x1v_ser_p1, NULL }, // 1
  { NULL, vlasov_pkpm_pressure_y_3x1v_ser_p1, NULL }, // 2
};

// PKPM pressure (in z) kernels (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_pkpm_pressure_kern_list ser_pkpm_pressure_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_pkpm_pressure_z_3x1v_ser_p1, NULL }, // 2
};

// PKPM pressure (in x) kernels (Tensor basis)
GKYL_CU_D
static const gkyl_dg_pkpm_pressure_kern_list ten_pkpm_pressure_x_kernels[] = {
  { NULL, vlasov_pkpm_pressure_x_1x1v_ser_p1, vlasov_pkpm_pressure_x_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_pkpm_pressure_x_2x1v_ser_p1, NULL }, // 1
  { NULL, vlasov_pkpm_pressure_x_3x1v_ser_p1, NULL }, // 2
};

// PKPM pressure (in y) kernels (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_pkpm_pressure_kern_list ten_pkpm_pressure_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_pkpm_pressure_y_2x1v_ser_p1, NULL }, // 1
  { NULL, vlasov_pkpm_pressure_y_3x1v_ser_p1, NULL }, // 2
};

// PKPM pressure (in z) kernels (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_pkpm_pressure_kern_list ten_pkpm_pressure_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_pkpm_pressure_z_3x1v_ser_p1, NULL }, // 2
};

GKYL_CU_D
static pkpm_prim_t
choose_ser_pkpm_prim_kern(int cdim, int poly_order)
{
  return ser_pkpm_prim_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static pkpm_source_t
choose_ser_pkpm_source_kern(int cdim, int poly_order)
{
  return ser_pkpm_source_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static pkpm_recovery_t
choose_ser_pkpm_recovery_kern(int dir, int cdim, int poly_order)
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

GKYL_CU_D
static pkpm_dist_mirror_force_t
choose_ten_pkpm_dist_mirror_force_kern(int cdim, int poly_order)
{
  return ten_pkpm_dist_mirror_force_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static pkpm_pressure_t
choose_ser_pkpm_pressure_kern(int dir, int cdim, int poly_order)
{
  if (dir == 0)
    return ser_pkpm_pressure_x_kernels[cdim-1].kernels[poly_order];
  else if (dir == 1)
    return ser_pkpm_pressure_y_kernels[cdim-1].kernels[poly_order];
  else if (dir == 2)
    return ser_pkpm_pressure_z_kernels[cdim-1].kernels[poly_order];
  else 
    return NULL;
}

GKYL_CU_D
static pkpm_pressure_t
choose_ten_pkpm_pressure_kern(int dir, int cdim, int poly_order)
{
  if (dir == 0)
    return ten_pkpm_pressure_x_kernels[cdim-1].kernels[poly_order];
  else if (dir == 1)
    return ten_pkpm_pressure_y_kernels[cdim-1].kernels[poly_order];
  else if (dir == 2)
    return ten_pkpm_pressure_z_kernels[cdim-1].kernels[poly_order];
  else 
    return NULL;
}
