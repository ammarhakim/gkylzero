#pragma once 
#include <math.h> 
#include <gkyl_util.h>
#include <gkyl_mom_calc_bcorr.h>

EXTERN_C_BEG 

GKYL_CU_DH void vlasov_bflux_1x1v_ser_p1(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_bflux_1x2v_ser_p1(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_bflux_1x3v_ser_p1(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_bflux_2x2v_ser_p1(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_bflux_2x3v_ser_p1(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_bflux_3x3v_ser_p1(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);

GKYL_CU_DH void vlasov_bflux_1x1v_ser_p2(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_bflux_1x2v_ser_p2(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_bflux_1x3v_ser_p2(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_bflux_2x2v_ser_p2(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);
GKYL_CU_DH void vlasov_bflux_2x3v_ser_p2(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out);

EXTERN_C_END 
