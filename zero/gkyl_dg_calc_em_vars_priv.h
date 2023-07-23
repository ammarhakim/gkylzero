// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_maxwell_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

typedef void (*em_calc_temp_t)(const double *em, double* out);
typedef int (*em_set_t)(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *temp);
typedef void (*em_copy_t)(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, double* out);

// for use in kernel tables
typedef struct { em_calc_temp_t kernels[3]; } gkyl_dg_em_calc_BB_kern_list;
typedef struct { em_calc_temp_t kernels[3]; } gkyl_dg_em_calc_num_ExB_kern_list;
typedef struct { em_set_t kernels[3]; } gkyl_dg_em_set_bvar_kern_list;
typedef struct { em_set_t kernels[3]; } gkyl_dg_em_set_ExB_kern_list;
typedef struct { em_copy_t kernels[3]; } gkyl_dg_em_copy_bvar_kern_list;
typedef struct { em_copy_t kernels[3]; } gkyl_dg_em_copy_ExB_kern_list;

struct gkyl_dg_calc_em_vars {
  struct gkyl_range conf_range; // Configuration space range

  struct gkyl_nmat *As, *xs; // matrices for LHS and RHS
  gkyl_nmat_mem *mem; // memory for use in batched linear solve
  struct gkyl_array *temp_var; // intermediate variable for storing weak multiplications
  int Ncomp; // number of components in the linear solve (3 for E x B, 6 for bb)
  int poly_order; // polynomial order (needed because only p>1 solves linear system)

  em_calc_temp_t em_calc_temp; // kernel for intermediate variable computation
  em_set_t em_set;  // kernel for setting matrices for linear solve
  em_copy_t em_copy; // kernel for copying solution to output 

  uint32_t flags;
  struct gkyl_dg_calc_em_vars *on_dev; // pointer to itself or device data
};

// Compute BB tensor for computing bb and b (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_em_calc_BB_kern_list ser_em_calc_BB_kernels[] = {
  { NULL, em_calc_BB_1x_ser_p1, em_calc_BB_1x_ser_p2 }, // 0
  { NULL, em_calc_BB_2x_ser_p1, NULL }, // 1
  { NULL, em_calc_BB_3x_ser_p1, NULL }, // 2
};

// Compute BB tensor for computing bb and b (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_calc_BB_kern_list ten_em_calc_BB_kernels[] = {
  { NULL, em_calc_BB_1x_ser_p1, em_calc_BB_1x_ser_p2 }, // 0
  { NULL, em_calc_BB_2x_ser_p1, em_calc_BB_2x_tensor_p2 }, // 1
  { NULL, em_calc_BB_3x_ser_p1, em_calc_BB_3x_tensor_p2 }, // 2
};

// Compute (E x B)_i and B_i^2 (numerator and denominator of E x B velocity) (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_em_calc_num_ExB_kern_list ser_em_calc_num_ExB_kernels[] = {
  { NULL, em_calc_num_ExB_1x_ser_p1, em_calc_num_ExB_1x_ser_p2 }, // 0
  { NULL, em_calc_num_ExB_2x_ser_p1, NULL }, // 1
  { NULL, em_calc_num_ExB_3x_ser_p1, NULL }, // 2
};

// Compute (E x B)_i and B_i^2 (numerator and denominator of E x B velocity) (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_calc_num_ExB_kern_list ten_em_calc_num_ExB_kernels[] = {
  { NULL, em_calc_num_ExB_1x_ser_p1, em_calc_num_ExB_1x_ser_p2 }, // 0
  { NULL, em_calc_num_ExB_2x_ser_p1, em_calc_num_ExB_2x_tensor_p2 }, // 1
  { NULL, em_calc_num_ExB_3x_ser_p1, em_calc_num_ExB_3x_tensor_p2 }, // 2
};

// Set matrices for computing bb, p=1 analytically solved (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_em_set_bvar_kern_list ser_em_set_bvar_kernels[] = {
  { NULL, em_set_bvar_1x_ser_p1, em_set_bvar_1x_ser_p2 }, // 0
  { NULL, em_set_bvar_2x_ser_p1, NULL }, // 1
  { NULL, em_set_bvar_3x_ser_p1, NULL }, // 2
};

// Set matrices for computing bb, p=1 analytically solved (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_set_bvar_kern_list ten_em_set_bvar_kernels[] = {
  { NULL, em_set_bvar_1x_ser_p1, em_set_bvar_1x_ser_p2 }, // 0
  { NULL, em_set_bvar_2x_ser_p1, em_set_bvar_2x_tensor_p2 }, // 1
  { NULL, em_set_bvar_3x_ser_p1, em_set_bvar_3x_tensor_p2 }, // 2
};

// Set matrices for computing ExB, p=1 analytically solved (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_em_set_ExB_kern_list ser_em_set_ExB_kernels[] = {
  { NULL, em_set_ExB_1x_ser_p1, em_set_ExB_1x_ser_p2 }, // 0
  { NULL, em_set_ExB_2x_ser_p1, NULL }, // 1
  { NULL, em_set_ExB_3x_ser_p1, NULL }, // 2
};

// Set matrices for computing ExB, p=1 analytically solved (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_set_ExB_kern_list ten_em_set_ExB_kernels[] = {
  { NULL, em_set_ExB_1x_ser_p1, em_set_ExB_1x_ser_p2 }, // 0
  { NULL, em_set_ExB_2x_ser_p1, em_set_ExB_2x_tensor_p2 }, // 1
  { NULL, em_set_ExB_3x_ser_p1, em_set_ExB_3x_tensor_p2 }, // 2
};

// Magnetic field unit vector and unit tensor kernel list copy solution (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_em_copy_bvar_kern_list ser_em_copy_bvar_kernels[] = {
  { NULL, em_copy_bvar_1x_ser_p1, em_copy_bvar_1x_ser_p2 }, // 0
  { NULL, em_copy_bvar_2x_ser_p1, NULL }, // 1
  { NULL, em_copy_bvar_3x_ser_p1, NULL }, // 2
};

// Magnetic field unit vector and unit tensor kernel list copy solution (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_copy_bvar_kern_list ten_em_copy_bvar_kernels[] = {
  { NULL, em_copy_bvar_1x_ser_p1, em_copy_bvar_1x_ser_p2 }, // 0
  { NULL, em_copy_bvar_2x_ser_p1, em_copy_bvar_2x_tensor_p2 }, // 1
  { NULL, em_copy_bvar_3x_ser_p1, em_copy_bvar_3x_tensor_p2 }, // 2
};

// E x B velocity kernel list copy solution (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_em_copy_ExB_kern_list ser_em_copy_ExB_kernels[] = {
  { NULL, em_copy_ExB_1x_ser_p1, em_copy_ExB_1x_ser_p2 }, // 0
  { NULL, em_copy_ExB_2x_ser_p1, NULL }, // 1
  { NULL, em_copy_ExB_3x_ser_p1, NULL }, // 2
};
// E x B velocity kernel list copy solution (Tensor basis)
GKYL_CU_D
static const gkyl_dg_em_copy_ExB_kern_list ten_em_copy_ExB_kernels[] = {
  { NULL, em_copy_ExB_1x_ser_p1, em_copy_ExB_1x_ser_p2 }, // 0
  { NULL, em_copy_ExB_2x_ser_p1, em_copy_ExB_2x_tensor_p2 }, // 1
  { NULL, em_copy_ExB_3x_ser_p1, em_copy_ExB_3x_tensor_p2 }, // 2
};

GKYL_CU_D
static em_calc_temp_t
choose_em_calc_BB_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_em_calc_BB_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_em_calc_BB_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static em_calc_temp_t
choose_em_calc_num_ExB_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_em_calc_num_ExB_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_em_calc_num_ExB_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static em_set_t
choose_em_set_bvar_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_em_set_bvar_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_em_set_bvar_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static em_set_t
choose_em_set_ExB_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_em_set_ExB_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_em_set_ExB_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static em_copy_t
choose_em_copy_bvar_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_em_copy_bvar_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_em_copy_bvar_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static em_copy_t
choose_em_copy_ExB_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_em_copy_ExB_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_em_copy_ExB_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}
