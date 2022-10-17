// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_sr_Gamma_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*sr_t)(const double *V_i, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1,  3,  4,  5}, // 2x kernel indices
  {-1,  6,  7,  8}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { sr_t kernels[3]; } gkyl_dg_sr_Gamma2_kern_list;
typedef struct { sr_t kernels[3]; } gkyl_dg_sr_Gamma_kern_list;
typedef struct { sr_t kernels[3]; } gkyl_dg_sr_Gamma_inv_kern_list;

// Lorentz boost factor *squared* kernel list
GKYL_CU_D
static const gkyl_dg_sr_Gamma2_kern_list ser_sr_Gamma2_kernels[] = {
  // 1x kernels
  { NULL, sr_Gamma2_1x1v_ser_p1, sr_Gamma2_1x1v_ser_p2 }, // 0
  { NULL, sr_Gamma2_1x2v_ser_p1, sr_Gamma2_1x2v_ser_p2 }, // 1
  { NULL, sr_Gamma2_1x3v_ser_p1, sr_Gamma2_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, sr_Gamma2_2x1v_ser_p1, sr_Gamma2_2x1v_ser_p2 }, // 3
  { NULL, sr_Gamma2_2x2v_ser_p1, sr_Gamma2_2x2v_ser_p2 }, // 4
  { NULL, sr_Gamma2_2x3v_ser_p1, sr_Gamma2_2x3v_ser_p2 }, // 5
  // 3x kernels
  { NULL, sr_Gamma2_3x1v_ser_p1, NULL                  }, // 6
  { NULL, sr_Gamma2_3x2v_ser_p1, NULL                  }, // 7
  { NULL, sr_Gamma2_3x3v_ser_p1, NULL                  }, // 8
};

// Lorentz boost factor kernel list
GKYL_CU_D
static const gkyl_dg_sr_Gamma_kern_list ser_sr_Gamma_kernels[] = {
  // 1x kernels
  { NULL, sr_Gamma_1x1v_ser_p1, sr_Gamma_1x1v_ser_p2 }, // 0
  { NULL, sr_Gamma_1x2v_ser_p1, sr_Gamma_1x2v_ser_p2 }, // 1
  { NULL, sr_Gamma_1x3v_ser_p1, sr_Gamma_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, sr_Gamma_2x1v_ser_p1, sr_Gamma_2x1v_ser_p2 }, // 3
  { NULL, sr_Gamma_2x2v_ser_p1, sr_Gamma_2x2v_ser_p2 }, // 4
  { NULL, sr_Gamma_2x3v_ser_p1, sr_Gamma_2x3v_ser_p2 }, // 5
  // 3x kernels
  { NULL, sr_Gamma_3x1v_ser_p1, NULL                 }, // 6
  { NULL, sr_Gamma_3x2v_ser_p1, NULL                 }, // 7
  { NULL, sr_Gamma_3x3v_ser_p1, NULL                 }, // 8
};

// *inverse* Lorentz boost factor kernel list
GKYL_CU_D
static const gkyl_dg_sr_Gamma_inv_kern_list ser_sr_Gamma_inv_kernels[] = {
  // 1x kernels
  { NULL, sr_Gamma_inv_1x1v_ser_p1, sr_Gamma_inv_1x1v_ser_p2 }, // 0
  { NULL, sr_Gamma_inv_1x2v_ser_p1, sr_Gamma_inv_1x2v_ser_p2 }, // 1
  { NULL, sr_Gamma_inv_1x3v_ser_p1, sr_Gamma_inv_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, sr_Gamma_inv_2x1v_ser_p1, sr_Gamma_inv_2x1v_ser_p2 }, // 3
  { NULL, sr_Gamma_inv_2x2v_ser_p1, sr_Gamma_inv_2x2v_ser_p2 }, // 4
  { NULL, sr_Gamma_inv_2x3v_ser_p1, sr_Gamma_inv_2x3v_ser_p2 }, // 5
  // 3x kernels
  { NULL, sr_Gamma_inv_3x1v_ser_p1, NULL                     }, // 6
  { NULL, sr_Gamma_inv_3x2v_ser_p1, NULL                     }, // 7
  { NULL, sr_Gamma_inv_3x3v_ser_p1, NULL                     }, // 8
};

GKYL_CU_D
static sr_t
choose_ser_sr_Gamma2_kern(int cdim, int vdim, int poly_order)
{
  return ser_sr_Gamma2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

GKYL_CU_D
static sr_t
choose_ser_sr_Gamma_kern(int cdim, int vdim, int poly_order)
{
  return ser_sr_Gamma_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

GKYL_CU_D
static sr_t
choose_ser_sr_Gamma_inv_kern(int cdim, int vdim, int poly_order)
{
  return ser_sr_Gamma_inv_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}