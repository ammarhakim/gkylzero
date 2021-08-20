#pragma once 
#include <math.h> 
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH void vlasov_lbo_constNu_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_lbo_constNu_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_lbo_constNu_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_lbo_constNu_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_lbo_constNu_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_lbo_constNu_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_lbo_constNu_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_lbo_constNu_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_lbo_constNu_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_lbo_constNu_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum,
  const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

EXTERN_C_END
