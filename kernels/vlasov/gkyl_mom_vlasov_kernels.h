#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void vlasov_M0_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_vlasov_pkpm_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_vlasov_pkpm_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_vlasov_pkpm_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_1x1v_ser_p3(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_2x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_vlasov_pkpm_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfy_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_vlasov_pkpm_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfy_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_vlasov_pkpm_2x1v_ser_p3(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_2x1v_ser_p3(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfy_2x1v_ser_p3(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_vlasov_pkpm_3x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_3x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfy_3x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfz_3x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_vlasov_pkpm_3x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_3x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfy_3x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfz_3x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_vlasov_pkpm_3x1v_ser_p3(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_3x1v_ser_p3(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfy_3x1v_ser_p3(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfz_3x1v_ser_p3(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_vlasov_pkpm_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_vlasov_pkpm_1x1v_tensor_p3(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_1x1v_tensor_p3(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_M0_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M1i_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M2ij_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3i_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_M3ijk_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_five_moments_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_int_mom_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void vlasov_sr_M0_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_M1i_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_int_mom_2x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_vlasov_pkpm_2x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_2x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfy_2x1v_tensor_p2(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void mom_vlasov_pkpm_2x1v_tensor_p3(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_2x1v_tensor_p3(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void mom_vlasov_pkpm_surfy_2x1v_tensor_p3(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out); 

EXTERN_C_END 
