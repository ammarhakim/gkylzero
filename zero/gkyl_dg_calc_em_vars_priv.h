// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_maxwell_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

typedef void (*em_basis_inv_t)(const double *em, double* out);
typedef int (*em_set_t)(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *em);
typedef void (*em_copy_t)(struct gkyl_mat *x, const double *em, int* cell_avg_magB2, double* out);

// for use in kernel tables
typedef struct { em_basis_inv_t kernels[3]; } gkyl_dg_em_bvar_basis_inv_kern_list;
typedef struct { em_basis_inv_t kernels[3]; } gkyl_dg_em_ExB_basis_inv_kern_list;
typedef struct { em_set_t kernels[3]; } gkyl_dg_em_set_magB2_kern_list;
typedef struct { em_copy_t kernels[3]; } gkyl_dg_em_copy_bvar_kern_list;
typedef struct { em_copy_t kernels[3]; } gkyl_dg_em_copy_ExB_kern_list;

// Magnetic field unit vector and unit tensor kernel list using basis_inv
GKYL_CU_D
static const gkyl_dg_em_bvar_basis_inv_kern_list ser_em_bvar_basis_inv_kernels[] = {
  { NULL, em_bvar_basis_inv_1x_ser_p1, em_bvar_basis_inv_1x_ser_p2 }, // 0
  { NULL, em_bvar_basis_inv_2x_ser_p1, NULL }, // 1
  { NULL, em_bvar_basis_inv_3x_ser_p1, NULL }, // 2
};

// E x B velocity kernel list
GKYL_CU_D
static const gkyl_dg_em_ExB_basis_inv_kern_list ser_em_ExB_basis_inv_kernels[] = {
  { NULL, em_ExB_basis_inv_1x_ser_p1, em_ExB_basis_inv_1x_ser_p2 }, // 0
  { NULL, em_ExB_basis_inv_2x_ser_p1, NULL }, // 1
  { NULL, em_ExB_basis_inv_3x_ser_p1, NULL }, // 2
};

// Set matrices for computing 1/|B|^2 (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_em_set_magB2_kern_list ser_em_set_magB2_kernels[] = {
  { NULL, em_set_magB2_1x_ser_p1, em_set_magB2_1x_ser_p2 }, // 0
  { NULL, em_set_magB2_2x_ser_p1, em_set_magB2_2x_ser_p2 }, // 1
  { NULL, em_set_magB2_3x_ser_p1, em_set_magB2_3x_ser_p2 }, // 2
};

// Set matrices for computing 1/|B|^2 (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_set_magB2_kern_list ten_em_set_magB2_kernels[] = {
  { NULL, em_set_magB2_1x_ser_p1, em_set_magB2_1x_ser_p2 }, // 0
  { NULL, em_set_magB2_2x_ser_p1, em_set_magB2_2x_tensor_p2 }, // 1
  { NULL, em_set_magB2_3x_ser_p1, em_set_magB2_3x_tensor_p2 }, // 2
};

// Magnetic field unit vector and unit tensor kernel list using weak division (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_em_copy_bvar_kern_list ser_em_copy_bvar_kernels[] = {
  { NULL, em_copy_bvar_1x_ser_p1, em_copy_bvar_1x_ser_p2 }, // 0
  { NULL, em_copy_bvar_2x_ser_p1, em_copy_bvar_2x_ser_p2 }, // 1
  { NULL, em_copy_bvar_3x_ser_p1, em_copy_bvar_3x_ser_p2 }, // 2
};

// Magnetic field unit vector and unit tensor kernel list using weak division (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_copy_bvar_kern_list ten_em_copy_bvar_kernels[] = {
  { NULL, em_copy_bvar_1x_ser_p1, em_copy_bvar_1x_ser_p2 }, // 0
  { NULL, em_copy_bvar_2x_ser_p1, em_copy_bvar_2x_tensor_p2 }, // 1
  { NULL, em_copy_bvar_3x_ser_p1, em_copy_bvar_3x_tensor_p2 }, // 2
};

// E x B velocity kernel list using weak division (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_em_copy_ExB_kern_list ser_em_copy_ExB_kernels[] = {
  { NULL, em_copy_ExB_1x_ser_p1, em_copy_ExB_1x_ser_p2 }, // 0
  { NULL, em_copy_ExB_2x_ser_p1, em_copy_ExB_2x_ser_p2 }, // 1
  { NULL, em_copy_ExB_3x_ser_p1, em_copy_ExB_3x_ser_p2  }, // 2
};
// E x B velocity kernel list using weak division (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_copy_ExB_kern_list ten_em_copy_ExB_kernels[] = {
  { NULL, em_copy_ExB_1x_ser_p1, em_copy_ExB_1x_ser_p2 }, // 0
  { NULL, em_copy_ExB_2x_ser_p1, em_copy_ExB_2x_tensor_p2 }, // 1
  { NULL, em_copy_ExB_3x_ser_p1, em_copy_ExB_3x_tensor_p2 }, // 2
};

GKYL_CU_D
static em_basis_inv_t
choose_ser_em_bvar_basis_inv_kern(int cdim, int poly_order)
{
  return ser_em_bvar_basis_inv_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_basis_inv_t
choose_ser_em_ExB_basis_inv_kern(int cdim, int poly_order)
{
  return ser_em_ExB_basis_inv_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_set_t
choose_ser_em_set_magB2_kern(int cdim, int poly_order)
{
  return ser_em_set_magB2_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_set_t
choose_ten_em_set_magB2_kern(int cdim, int poly_order)
{
  return ten_em_set_magB2_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_copy_t
choose_ser_em_copy_bvar_kern(int cdim, int poly_order)
{
  return ser_em_copy_bvar_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_copy_t
choose_ten_em_copy_bvar_kern(int cdim, int poly_order)
{
  return ten_em_copy_bvar_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_copy_t
choose_ser_em_copy_ExB_kern(int cdim, int poly_order)
{
  return ser_em_copy_ExB_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_copy_t
choose_ten_em_copy_ExB_kern(int cdim, int poly_order)
{
  return ten_em_copy_ExB_kernels[cdim-1].kernels[poly_order];
}
