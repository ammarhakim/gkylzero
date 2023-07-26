#pragma once 

#include <math.h> 
#include <gkyl_util.h> 
 
EXTERN_C_BEG 
 
GKYL_CU_DH double dg_euleriso_diffusion_vol_1x_ser_p1(const double *w, const double *dxv, const double *D, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_euleriso_diffusion_surfx_1x_ser_p1(const double *w, const double *dxv, const double *D, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_euleriso_diffusion_vol_1x_ser_p2(const double *w, const double *dxv, const double *D, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_euleriso_diffusion_surfx_1x_ser_p2(const double *w, const double *dxv, const double *D, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_euleriso_diffusion_vol_2x_ser_p1(const double *w, const double *dxv, const double *D, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_euleriso_diffusion_surfx_2x_ser_p1(const double *w, const double *dxv, const double *D, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_euleriso_diffusion_surfy_2x_ser_p1(const double *w, const double *dxv, const double *D, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_euleriso_diffusion_vol_2x_ser_p2(const double *w, const double *dxv, const double *D, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_euleriso_diffusion_surfx_2x_ser_p2(const double *w, const double *dxv, const double *D, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_euleriso_diffusion_surfy_2x_ser_p2(const double *w, const double *dxv, const double *D, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_euleriso_diffusion_vol_3x_ser_p1(const double *w, const double *dxv, const double *D, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_euleriso_diffusion_surfx_3x_ser_p1(const double *w, const double *dxv, const double *D, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_euleriso_diffusion_surfy_3x_ser_p1(const double *w, const double *dxv, const double *D, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_euleriso_diffusion_surfz_3x_ser_p1(const double *w, const double *dxv, const double *D, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_euleriso_diffusion_vol_3x_ser_p2(const double *w, const double *dxv, const double *D, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_euleriso_diffusion_surfx_3x_ser_p2(const double *w, const double *dxv, const double *D, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_euleriso_diffusion_surfy_3x_ser_p2(const double *w, const double *dxv, const double *D, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_euleriso_diffusion_surfz_3x_ser_p2(const double *w, const double *dxv, const double *D, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

EXTERN_C_END 
