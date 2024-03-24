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

// Projection functions for p/(gamma) = v in special relativistic systems
// Simplifies to p/sqrt(1 + p^2) where c = 1
static void 
ev_p_over_gamma_1p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = xn[0]/sqrt(1.0 + xn[0]*xn[0]);
}
static void 
ev_p_over_gamma_2p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = xn[0]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1]);
  out[1] = xn[1]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1]);
}
static void 
ev_p_over_gamma_3p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = xn[0]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
  out[1] = xn[1]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
  out[2] = xn[2]/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
}
static const evalf_t p_over_gamma_func[3] = {ev_p_over_gamma_1p, ev_p_over_gamma_2p, ev_p_over_gamma_3p};

// Projection functions for gamma = sqrt(1 + p^2) in special relativistic systems
static void 
ev_gamma_1p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = sqrt(1.0 + xn[0]*xn[0]);
}
static void 
ev_gamma_2p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1]);
}
static void 
ev_gamma_3p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
}
static const evalf_t gamma_func[3] = {ev_gamma_1p, ev_gamma_2p, ev_gamma_3p};

// Projection functions for gamma_inv = 1/sqrt(1 + p^2) in special relativistic systems
static void 
ev_gamma_inv_1p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = 1.0/sqrt(1.0 + xn[0]*xn[0]);
}
static void 
ev_gamma_inv_2p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = 1.0/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1]);
}
static void 
ev_gamma_inv_3p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = 1.0/sqrt(1.0 + xn[0]*xn[0] + xn[1]*xn[1] + xn[2]*xn[2]);
}
static const evalf_t gamma_inv_func[3] = {ev_gamma_inv_1p, ev_gamma_inv_2p, ev_gamma_inv_3p};

// Function pointer type for different GammaV functions (GammaV^2, GammaV, 1/GammaV)
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
typedef struct { sr_t kernels[3]; } gkyl_dg_sr_GammaV2_kern_list;
typedef struct { sr_t kernels[3]; } gkyl_dg_sr_GammaV_inv_kern_list;

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
