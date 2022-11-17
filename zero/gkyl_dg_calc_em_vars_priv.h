// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_maxwell_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*em_t)(const double *em, double* out);
typedef void (*em_pkpm_kappa_inv_b_t)(const double *bvar, const double *ExB, double* out);

// for use in kernel tables
typedef struct { em_t kernels[3]; } gkyl_dg_em_bvar_kern_list;
typedef struct { em_t kernels[3]; } gkyl_dg_em_ExB_kern_list;
typedef struct { em_pkpm_kappa_inv_b_t kernels[3]; } gkyl_dg_em_pkpm_kappa_inv_b_kern_list;

// Magnetic field unit vector and unit tensor kernel list
GKYL_CU_D
static const gkyl_dg_em_bvar_kern_list ser_em_bvar_kernels[] = {
  { NULL, em_bvar_1x_ser_p1, em_bvar_1x_ser_p2 }, // 0
  { NULL, em_bvar_2x_ser_p1, em_bvar_2x_ser_p2 }, // 1
  { NULL, em_bvar_3x_ser_p1, NULL }, // 2
};

// E x B velocity kernel list
GKYL_CU_D
static const gkyl_dg_em_ExB_kern_list ser_em_ExB_kernels[] = {
  { NULL, em_ExB_1x_ser_p1, em_ExB_1x_ser_p2 }, // 0
  { NULL, em_ExB_2x_ser_p1, em_ExB_2x_ser_p2 }, // 1
  { NULL, em_ExB_3x_ser_p1, NULL }, // 2
};

// b_hat/kappa (1/kappa = 1 - (E x B)^2/(c^2 |B|^4)) for pkpm model kernel list
GKYL_CU_D
static const gkyl_dg_em_pkpm_kappa_inv_b_kern_list ser_em_pkpm_kappa_inv_b_kernels[] = {
  { NULL, em_pkpm_kappa_inv_b_1x_ser_p1, em_pkpm_kappa_inv_b_1x_ser_p2 }, // 0
  { NULL, em_pkpm_kappa_inv_b_2x_ser_p1, em_pkpm_kappa_inv_b_2x_ser_p2 }, // 1
  { NULL, em_pkpm_kappa_inv_b_3x_ser_p1, NULL }, // 2
};

GKYL_CU_D
static em_t
choose_ser_em_bvar_kern(int cdim, int poly_order)
{
  return ser_em_bvar_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_t
choose_ser_em_ExB_kern(int cdim, int poly_order)
{
  return ser_em_ExB_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static em_pkpm_kappa_inv_b_t
choose_ser_em_pkpm_kappa_inv_b_kern(int cdim, int poly_order)
{
  return ser_em_pkpm_kappa_inv_b_kernels[cdim-1].kernels[poly_order];
}