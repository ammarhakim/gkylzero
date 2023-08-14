#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void euler_pressure_1x_ser_p1(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_1x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfx_1x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_1x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_iso_surfx_1x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pressure_1x_ser_p2(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_1x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfx_1x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_1x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_iso_surfx_1x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pressure_2x_ser_p1(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_2x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfx_2x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfy_2x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_2x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_iso_surfx_2x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_iso_surfy_2x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pressure_2x_ser_p2(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_2x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfx_2x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfy_2x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_2x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_iso_surfx_2x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_iso_surfy_2x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pressure_3x_ser_p1(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_3x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfx_3x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfy_3x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfz_3x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_3x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_iso_surfx_3x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_iso_surfy_3x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_iso_surfz_3x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pressure_3x_ser_p2(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_3x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfx_3x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfy_3x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_surfz_3x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_3x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_iso_surfx_3x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_iso_surfy_3x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double euler_iso_surfz_3x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

EXTERN_C_END 
