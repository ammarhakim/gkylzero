#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double vlasov_poisson_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_1x2v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_1x3v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x3v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_1x3v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x3v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfy_2x2v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfy_2x2v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_2x3v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfy_2x3v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfy_2x3v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_2x3v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_2x3v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfy_2x3v_ser_p2(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfy_2x3v_ser_p2(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfx_3x3v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_3x3v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfy_3x3v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfy_3x3v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfz_3x3v_ser_p1(const double *w, const double *dxv, 
      const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfz_3x3v_ser_p1(const double *w, const double *dxv, 
      const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_boundary_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_boundary_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

EXTERN_C_END 