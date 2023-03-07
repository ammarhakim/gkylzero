#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double vlasov_poisson_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_stream_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_gen_geo_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *cot_vec, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_stream_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_gen_geo_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *cot_vec, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_stream_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_gen_geo_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *cot_vec, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_stream_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_gen_geo_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *cot_vec, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfx_1x2v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfx_1x2v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfx_1x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_stream_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_gen_geo_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *cot_vec, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfx_1x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfx_1x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfx_1x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfx_1x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_stream_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_gen_geo_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *cot_vec, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfx_1x3v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfx_1x3v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfx_1x3v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfx_1x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvy_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvz_1x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *bvar, const double *u_i, 
  const double *pkpm_accel_vars, const double *g_dist_source, 
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, const double *vth_sql, const double *vth_sqc, const double *vth_sqr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
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
GKYL_CU_DH void vlasov_pkpm_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, const double *vth_sql, const double *vth_sqc, const double *vth_sqr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_stream_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_gen_geo_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *cot_vec, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_stream_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_gen_geo_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *cot_vec, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_stream_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_gen_geo_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *cot_vec, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfx_2x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfx_2x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfx_2x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfx_2x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfy_2x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfy_2x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfy_2x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfy_2x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_stream_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_gen_geo_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *cot_vec, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfx_2x3v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfx_2x3v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfx_2x3v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfx_2x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfy_2x3v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfy_2x3v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfy_2x3v_ser_p2(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfy_2x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_2x1v_ser_p1(const double *w, const double *dxv, 
  const double *bvar, const double *u_i, 
  const double *pkpm_accel_vars, const double *g_dist_source, 
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfx_2x1v_ser_p1(const double *w, const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, const double *vth_sql, const double *vth_sqc, const double *vth_sqr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfy_2x1v_ser_p1(const double *w, const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, const double *vth_sql, const double *vth_sqc, const double *vth_sqr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
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
GKYL_CU_DH void vlasov_pkpm_surfx_2x1v_ser_p2(const double *w, const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, const double *vth_sql, const double *vth_sqc, const double *vth_sqr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfy_2x1v_ser_p2(const double *w, const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, const double *vth_sql, const double *vth_sqc, const double *vth_sqr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfvpar_2x1v_ser_p2(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_2x1v_ser_p2(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_poisson_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_poisson_extem_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_stream_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_stream_gen_geo_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *cot_vec, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_sr_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double vlasov_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfx_3x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfx_3x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfx_3x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfx_3x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfy_3x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfy_3x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfy_3x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfy_3x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfz_3x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfz_3x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_gen_geo_surfz_3x3v_ser_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfz_3x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_poisson_boundary_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_sr_boundary_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_boundary_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *ext_field, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_3x1v_ser_p1(const double *w, const double *dxv, 
  const double *bvar, const double *u_i, 
  const double *pkpm_accel_vars, const double *g_dist_source, 
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfx_3x1v_ser_p1(const double *w, const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, const double *vth_sql, const double *vth_sqc, const double *vth_sqr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfy_3x1v_ser_p1(const double *w, const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, const double *vth_sql, const double *vth_sqc, const double *vth_sqr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfz_3x1v_ser_p1(const double *w, const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, const double *vth_sql, const double *vth_sqc, const double *vth_sqr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfvpar_3x1v_ser_p1(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_3x1v_ser_p1(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double vlasov_pkpm_vol_3x1v_ser_p2(const double *w, const double *dxv, 
  const double *bvar, const double *u_i, 
  const double *pkpm_accel_vars, const double *g_dist_source, 
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfx_3x1v_ser_p2(const double *w, const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, const double *vth_sql, const double *vth_sqc, const double *vth_sqr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfy_3x1v_ser_p2(const double *w, const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, const double *vth_sql, const double *vth_sqc, const double *vth_sqr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfz_3x1v_ser_p2(const double *w, const double *dxv, const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, const double *vth_sql, const double *vth_sqc, const double *vth_sqr, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_surfvpar_3x1v_ser_p2(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_3x1v_ser_p2(const double *w, const double *dxv, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

EXTERN_C_END 
