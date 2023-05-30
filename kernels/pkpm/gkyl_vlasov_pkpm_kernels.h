#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double vlasov_pkpm_vol_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *bvar, const double *u_i, 
  const double *pkpm_accel_vars, const double *g_dist_source, 
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *T_ijl, const double *T_ijc, const double *T_ijr, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *bvar, const double *u_i, 
  const double *pkpm_accel_vars, const double *g_dist_source, 
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *T_ijl, const double *T_ijc, const double *T_ijr, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_2x1v_ser_p1(const double *w, const double *dxv, 
  const double *bvar, const double *u_i, 
  const double *pkpm_accel_vars, const double *g_dist_source, 
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfx_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *T_ijl, const double *T_ijc, const double *T_ijr, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfy_2x1v_ser_p1(const double *w, const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *T_ijl, const double *T_ijc, const double *T_ijr, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_2x1v_ser_p2(const double *w, const double *dxv, 
  const double *bvar, const double *u_i, 
  const double *pkpm_accel_vars, const double *g_dist_source, 
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfx_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *T_ijl, const double *T_ijc, const double *T_ijr, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfy_2x1v_ser_p2(const double *w, const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *T_ijl, const double *T_ijc, const double *T_ijr, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfvpar_2x1v_ser_p2(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_2x1v_ser_p2(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_3x1v_ser_p1(const double *w, const double *dxv, 
  const double *bvar, const double *u_i, 
  const double *pkpm_accel_vars, const double *g_dist_source, 
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfx_3x1v_ser_p1(const double *w, const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *T_ijl, const double *T_ijc, const double *T_ijr, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfy_3x1v_ser_p1(const double *w, const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *T_ijl, const double *T_ijc, const double *T_ijr, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfz_3x1v_ser_p1(const double *w, const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *T_ijl, const double *T_ijc, const double *T_ijr, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfvpar_3x1v_ser_p1(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_3x1v_ser_p1(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

EXTERN_C_END 
