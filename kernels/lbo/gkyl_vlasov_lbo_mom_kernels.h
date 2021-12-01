#pragma once
#include <gkyl_mom_bcorr.h>
#include <math.h> 
#include <gkyl_util.h>

EXTERN_C_BEG 

GKYL_CU_DH void vlasov_boundary_integral_1x1v_ser_f_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_1x1v_ser_vf_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_1x1v_ser_f_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_1x1v_ser_vf_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_boundary_integral_1x2v_ser_f_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_1x2v_ser_vf_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_1x2v_ser_f_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_1x2v_ser_vf_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_boundary_integral_1x3v_ser_f_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_1x3v_ser_vf_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_1x3v_ser_f_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_1x3v_ser_vf_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_boundary_integral_2x2v_ser_f_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_2x2v_ser_vf_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_2x2v_ser_f_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_2x2v_ser_vf_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_boundary_integral_2x3v_ser_f_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_2x3v_ser_vf_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_boundary_integral_3x3v_ser_f_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_boundary_integral_3x3v_ser_vf_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);

EXTERN_C_END 
