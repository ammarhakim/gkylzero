// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_maxwell_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>
#include <assert.h>

typedef void (*em_calc_temp_t)(const double *em, double* GKYL_RESTRICT out);

typedef void (*em_set_t)(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *temp);

typedef void (*em_copy_t)(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_bb, 
  double* GKYL_RESTRICT out, double* GKYL_RESTRICT out_surf);

typedef void (*em_set_diag_t)(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *em, const double *temp);

typedef void (*em_copy_diag_t)(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_bb, 
  double* GKYL_RESTRICT out);

typedef void (*em_div_b_t)(const double *dxv, 
  const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
  const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b); 

typedef void (*em_limiter_t)(double limiter_fac, 
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
  double *ql, double *qc, double *qr);

// for use in kernel tables
typedef struct { em_calc_temp_t kernels[3]; } gkyl_dg_em_calc_BB_kern_list;
typedef struct { em_set_t kernels[3]; } gkyl_dg_em_set_bvar_kern_list;
typedef struct { em_copy_t kernels[3]; } gkyl_dg_em_copy_bvar_kern_list;
typedef struct { em_set_diag_t kernels[3]; } gkyl_dg_em_set_diag_kern_list;
typedef struct { em_copy_diag_t kernels[3]; } gkyl_dg_em_copy_diag_kern_list;
typedef struct { em_div_b_t kernels[3]; } gkyl_dg_em_div_b_kern_list;
typedef struct { em_limiter_t kernels[3]; } gkyl_dg_em_limiter_kern_list;

struct gkyl_dg_calc_em_vars {
  struct gkyl_rect_grid conf_grid; // Configuration space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  struct gkyl_range mem_range; // Configuration space range for linear solve

  const struct gkyl_wv_eqn *wv_eqn; // Wave equation for characteristic limiting of solution
  const struct gkyl_wave_geom *geom; // Wave geometry for rotating solution
  double limiter_fac; // Factor for relationship between cell slopes and cell average differences (by default: 1/sqrt(3))

  struct gkyl_array *temp_var; // intermediate variable for storing weak multiplications of B_i B_j on order 2*p basis

  struct gkyl_nmat *As, *xs; // matrices for LHS and RHS
  gkyl_nmat_mem *mem; // memory for use in batched linear solve
  int Ncomp; // number of components in the linear solve (6 for bb)

  struct gkyl_nmat *As_diag, *xs_diag; // matrices for LHS and RHS for diagnostic variable solve
  gkyl_nmat_mem *mem_diag; // memory for use in batched linear solve for diagnostic variable solve
  int Ncomp_diag; // number of components in the linear solve (9, 6 for bb + 3 for E x B)

  em_calc_temp_t em_calc_temp; // kernel for intermediate variable computation
  em_set_t em_set;  // kernel for setting matrices for linear solve for bvar
  em_copy_t em_copy; // kernel for copying solution to output for bvar, also computes needed surface expansions

  em_set_diag_t em_diag_set;  // kernel for setting matrices for linear solve for diagnostic variables 
  em_copy_diag_t em_diag_copy; // kernel for copying solution to output for diagnostic variables 

  em_div_b_t em_div_b[3]; // kernel for computing div(b) and max(|b_i|) penalization
  em_limiter_t em_limiter[3]; // kernel for limiting slopes of em variables

  uint32_t flags;
  struct gkyl_dg_calc_em_vars *on_dev; // pointer to itself or device data
};

// Compute BB tensor for computing bb and b (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_calc_BB_kern_list ten_em_calc_BB_kernels[] = {
  { NULL, NULL, em_calc_BB_1x_tensor_p2 }, // 0
  { NULL, NULL, em_calc_BB_2x_tensor_p2 }, // 1
  { NULL, NULL, em_calc_BB_3x_tensor_p2 }, // 2
};

// Set matrices for computing bb (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_set_bvar_kern_list ten_em_set_bvar_kernels[] = {
  { NULL, NULL, em_set_bvar_1x_tensor_p2 }, // 0
  { NULL, NULL, em_set_bvar_2x_tensor_p2 }, // 1
  { NULL, NULL, em_set_bvar_3x_tensor_p2 }, // 2
};

// Magnetic field unit vector and unit tensor kernel list copy solution (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_copy_bvar_kern_list ten_em_copy_bvar_kernels[] = {
  { NULL, NULL, em_copy_bvar_1x_tensor_p2 }, // 0
  { NULL, NULL, em_copy_bvar_2x_tensor_p2 }, // 1
  { NULL, NULL, em_copy_bvar_3x_tensor_p2 }, // 2
};

// Set matrices for computing diagnostic EM variables (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_set_diag_kern_list ten_em_set_diag_kernels[] = {
  { NULL, NULL, em_set_diag_1x_tensor_p2 }, // 0
  { NULL, NULL, em_set_diag_2x_tensor_p2 }, // 1
  { NULL, NULL, em_set_diag_3x_tensor_p2 }, // 2
};

// EM variables kernel list copy solution (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_copy_diag_kern_list ten_em_copy_diag_kernels[] = {
  { NULL, NULL, em_copy_diag_1x_tensor_p2 }, // 0
  { NULL, NULL, em_copy_diag_2x_tensor_p2 }, // 1
  { NULL, NULL, em_copy_diag_3x_tensor_p2 }, // 2
};

// div(b) and max(|b_i|) penalization (in x) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_em_div_b_kern_list ten_em_div_b_x_kernels[] = {
  { NULL, NULL, em_div_b_x_1x_tensor_p2 }, // 0
  { NULL, NULL, em_div_b_x_2x_tensor_p2 }, // 1
  { NULL, NULL, em_div_b_x_3x_tensor_p2 }, // 2
};

// div(b) and max(|b_i|) penalization (in y) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_em_div_b_kern_list ten_em_div_b_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, em_div_b_y_2x_tensor_p2 }, // 1
  { NULL, NULL, em_div_b_y_3x_tensor_p2 }, // 2
};

// div(b) and max(|b_i|) penalization (in z) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_em_div_b_kern_list ten_em_div_b_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, em_div_b_z_3x_tensor_p2 }, // 2
};

// Characteristic limiter in x (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_em_limiter_kern_list ten_em_limiter_x_kernels[] = {
  { NULL, em_vars_limiter_x_1x_tensor_p1, NULL }, // 0
  { NULL, em_vars_limiter_x_2x_tensor_p1, NULL }, // 1
  { NULL, em_vars_limiter_x_3x_tensor_p1, NULL }, // 2
};

// Characteristic limiter in y (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_em_limiter_kern_list ten_em_limiter_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, em_vars_limiter_y_2x_tensor_p1, NULL }, // 1
  { NULL, em_vars_limiter_y_3x_tensor_p1, NULL }, // 2
};

// Characteristic limiter in z (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_em_limiter_kern_list ten_em_limiter_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, em_vars_limiter_z_3x_tensor_p1, NULL }, // 2
};

GKYL_CU_D
static em_calc_temp_t
choose_em_calc_BB_kern(int cdim, int poly_order)
{
  return ten_em_calc_BB_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_set_t
choose_em_set_bvar_kern(int cdim, int poly_order)
{
  return ten_em_set_bvar_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_copy_t
choose_em_copy_bvar_kern(int cdim, int poly_order)
{
  return ten_em_copy_bvar_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_set_diag_t
choose_em_set_diag_kern(int cdim, int poly_order)
{
  return ten_em_set_diag_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_copy_diag_t
choose_em_copy_diag_kern(int cdim, int poly_order)
{
  return ten_em_copy_diag_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_div_b_t
choose_em_div_b_kern(int dir, int cdim, int poly_order)
{
  if (dir == 0)
    return ten_em_div_b_x_kernels[cdim-1].kernels[poly_order];
  else if (dir == 1)
    return ten_em_div_b_y_kernels[cdim-1].kernels[poly_order];
  else if (dir == 2)
    return ten_em_div_b_z_kernels[cdim-1].kernels[poly_order];
  else
    return NULL;
}

GKYL_CU_D
static em_limiter_t
choose_em_limiter_kern(int dir, int cdim, int poly_order)
{
  if (dir == 0)
    return ten_em_limiter_x_kernels[cdim-1].kernels[poly_order];
  else if (dir == 1)
    return ten_em_limiter_y_kernels[cdim-1].kernels[poly_order];
  else if (dir == 2)
    return ten_em_limiter_z_kernels[cdim-1].kernels[poly_order];
  else
    return NULL;  
}
