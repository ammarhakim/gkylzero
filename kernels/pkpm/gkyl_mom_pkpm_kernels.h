#pragma once 
#include <math.h> 
#include <gkyl_eqn_type.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void mom_pkpm_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_pkpm_diag_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_pkpm_1x1v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, double mass, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_self_prim_moments_1x1v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *vlasov_pkpm_moms, const double *boundary_corrections); 

GKYL_CU_DH void mom_pkpm_2x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_pkpm_diag_2x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_pkpm_2x1v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, double mass, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_self_prim_moments_2x1v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *vlasov_pkpm_moms, const double *boundary_corrections); 

GKYL_CU_DH void mom_pkpm_3x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_pkpm_diag_3x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_bcorr_lbo_pkpm_3x1v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, double mass, const double *fIn, double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_self_prim_moments_3x1v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *vlasov_pkpm_moms, const double *boundary_corrections); 

EXTERN_C_END 
