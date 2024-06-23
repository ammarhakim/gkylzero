// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_sr_Gamma_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*p_vars_t)(const double *w, const double *dv, 
  double* GKYL_RESTRICT gamma, double* GKYL_RESTRICT gamma_inv);

// Function pointer type for different GammaV functions (GammaV^2, 1/GammaV)
typedef void (*sr_t)(const double *V_i, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
GKYL_CU_D
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { p_vars_t kernels[3]; } gkyl_dg_sr_p_vars_kern_list;
typedef struct { sr_t kernels[3]; } gkyl_dg_sr_GammaV2_kern_list;
typedef struct { sr_t kernels[3]; } gkyl_dg_sr_GammaV_inv_kern_list;

// Particle Lorentz boost factor gamma = sqrt(1 + p^2) (also 1/gamma) kernel list
GKYL_CU_D
static const gkyl_dg_sr_p_vars_kern_list ser_sr_p_vars_kernels[] = {
  // 1x kernels
  { NULL, NULL, sr_vars_lorentz_1v_ser_p2 }, // 0
  { NULL, NULL, sr_vars_lorentz_2v_ser_p2 }, // 1
  { NULL, NULL, sr_vars_lorentz_3v_ser_p2 }, // 2
};

// Lorentz boost factor *squared* kernel list
GKYL_CU_D
static const gkyl_dg_sr_GammaV2_kern_list ser_sr_GammaV2_kernels[] = {
  // 1x kernels
  { NULL, sr_Gamma2_1x1v_ser_p1, sr_Gamma2_1x1v_ser_p2 }, // 0
  { NULL, sr_Gamma2_1x2v_ser_p1, sr_Gamma2_1x2v_ser_p2 }, // 1
  { NULL, sr_Gamma2_1x3v_ser_p1, sr_Gamma2_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, sr_Gamma2_2x2v_ser_p1, NULL }, // 3
  { NULL, sr_Gamma2_2x3v_ser_p1, NULL }, // 4
  // 3x kernels
  { NULL, sr_Gamma2_3x3v_ser_p1, NULL }, // 5
};

// *inverse* Lorentz boost factor kernel list
GKYL_CU_D
static const gkyl_dg_sr_GammaV_inv_kern_list ser_sr_GammaV_inv_kernels[] = {
  // 1x kernels
  { NULL, sr_Gamma_inv_1x1v_ser_p1, sr_Gamma_inv_1x1v_ser_p2 }, // 0
  { NULL, sr_Gamma_inv_1x2v_ser_p1, sr_Gamma_inv_1x2v_ser_p2 }, // 1
  { NULL, sr_Gamma_inv_1x3v_ser_p1, sr_Gamma_inv_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, sr_Gamma_inv_2x2v_ser_p1, NULL }, // 3
  { NULL, sr_Gamma_inv_2x3v_ser_p1, NULL }, // 4
  // 3x kernels
  { NULL, sr_Gamma_inv_3x3v_ser_p1, NULL }, // 5
};

GKYL_CU_D
static p_vars_t
choose_ser_sr_p_vars_kern(int vdim, int poly_order)
{
  return ser_sr_p_vars_kernels[vdim-1].kernels[poly_order];
}

GKYL_CU_D
static sr_t
choose_ser_sr_GammaV2_kern(int cdim, int vdim, int poly_order)
{
  return ser_sr_GammaV2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

GKYL_CU_D
static sr_t
choose_ser_sr_GammaV_inv_kern(int cdim, int vdim, int poly_order)
{
  return ser_sr_GammaV_inv_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}
