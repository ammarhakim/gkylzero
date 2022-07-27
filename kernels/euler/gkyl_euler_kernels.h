#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double euler_vol_1x_ser_p1(const double *w, const double *dxv, double gas_gamma, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_1x_ser_p1(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_vol_1x_ser_p2(const double *w, const double *dxv, double gas_gamma, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_1x_ser_p2(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_vol_2x_ser_p1(const double *w, const double *dxv, double gas_gamma, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_2x_ser_p1(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfy_2x_ser_p1(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_vol_2x_ser_p2(const double *w, const double *dxv, double gas_gamma, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_2x_ser_p2(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfy_2x_ser_p2(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_vol_3x_ser_p1(const double *w, const double *dxv, double gas_gamma, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_3x_ser_p1(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfy_3x_ser_p1(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfz_3x_ser_p1(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_vol_3x_ser_p2(const double *w, const double *dxv, double gas_gamma, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_3x_ser_p2(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfy_3x_ser_p2(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfz_3x_ser_p2(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_vol_1x_tensor_p2(const double *w, const double *dxv, double gas_gamma, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_1x_tensor_p2(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_vol_2x_tensor_p2(const double *w, const double *dxv, double gas_gamma, const double *u, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_2x_tensor_p2(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfy_2x_tensor_p2(const double *w, const double *dxv, double gas_gamma, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

EXTERN_C_END 
