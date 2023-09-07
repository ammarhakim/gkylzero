#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 
typedef struct { double c, chi, gamma; } gkyl_maxwell_inp; 

GKYL_CU_DH double maxwell_vol_1x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_1x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void em_bvar_1x_ser_p1(const double *em, double* bvar); 
GKYL_CU_DH void em_set_bvar_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *em); 
GKYL_CU_DH void em_copy_bvar_1x_ser_p1(int count, struct gkyl_nmat *x, const double *em, double* bvar); 
GKYL_CU_DH void em_ExB_1x_ser_p1(const double *em, double* ExB); 
GKYL_CU_DH void em_pkpm_kappa_inv_b_1x_ser_p1(const double *bvar, const double *ExB, double* kappa_inv_b); 

GKYL_CU_DH double maxwell_vol_1x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_1x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void em_bvar_1x_ser_p2(const double *em, double* bvar); 
GKYL_CU_DH void em_set_bvar_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *em); 
GKYL_CU_DH void em_copy_bvar_1x_ser_p2(int count, struct gkyl_nmat *x, const double *em, double* bvar); 
GKYL_CU_DH void em_ExB_1x_ser_p2(const double *em, double* ExB); 
GKYL_CU_DH void em_pkpm_kappa_inv_b_1x_ser_p2(const double *bvar, const double *ExB, double* kappa_inv_b); 

GKYL_CU_DH double maxwell_vol_2x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_2x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfy_2x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void em_bvar_2x_ser_p1(const double *em, double* bvar); 
GKYL_CU_DH void em_set_bvar_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *em); 
GKYL_CU_DH void em_copy_bvar_2x_ser_p1(int count, struct gkyl_nmat *x, const double *em, double* bvar); 
GKYL_CU_DH void em_ExB_2x_ser_p1(const double *em, double* ExB); 
GKYL_CU_DH void em_pkpm_kappa_inv_b_2x_ser_p1(const double *bvar, const double *ExB, double* kappa_inv_b); 

GKYL_CU_DH double maxwell_vol_2x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_2x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfy_2x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void em_bvar_2x_ser_p2(const double *em, double* bvar); 
GKYL_CU_DH void em_set_bvar_2x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *em); 
GKYL_CU_DH void em_copy_bvar_2x_ser_p2(int count, struct gkyl_nmat *x, const double *em, double* bvar); 
GKYL_CU_DH void em_ExB_2x_ser_p2(const double *em, double* ExB); 
GKYL_CU_DH void em_pkpm_kappa_inv_b_2x_ser_p2(const double *bvar, const double *ExB, double* kappa_inv_b); 

GKYL_CU_DH double maxwell_vol_3x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_3x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfy_3x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfz_3x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void em_bvar_3x_ser_p1(const double *em, double* bvar); 
GKYL_CU_DH void em_set_bvar_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *em); 
GKYL_CU_DH void em_copy_bvar_3x_ser_p1(int count, struct gkyl_nmat *x, const double *em, double* bvar); 
GKYL_CU_DH void em_ExB_3x_ser_p1(const double *em, double* ExB); 
GKYL_CU_DH void em_pkpm_kappa_inv_b_3x_ser_p1(const double *bvar, const double *ExB, double* kappa_inv_b); 

GKYL_CU_DH double maxwell_vol_1x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_1x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void em_bvar_1x_tensor_p2(const double *em, double* bvar); 
GKYL_CU_DH void em_set_bvar_1x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *em); 
GKYL_CU_DH void em_copy_bvar_1x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, double* bvar); 
GKYL_CU_DH void em_ExB_1x_tensor_p2(const double *em, double* ExB); 
GKYL_CU_DH void em_pkpm_kappa_inv_b_1x_tensor_p2(const double *bvar, const double *ExB, double* kappa_inv_b); 

GKYL_CU_DH double maxwell_vol_2x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_2x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfy_2x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void em_bvar_2x_tensor_p2(const double *em, double* bvar); 
GKYL_CU_DH void em_set_bvar_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *em); 
GKYL_CU_DH void em_copy_bvar_2x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, double* bvar); 
GKYL_CU_DH void em_ExB_2x_tensor_p2(const double *em, double* ExB); 
GKYL_CU_DH void em_pkpm_kappa_inv_b_2x_tensor_p2(const double *bvar, const double *ExB, double* kappa_inv_b); 

EXTERN_C_END 
