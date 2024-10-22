#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void canonical_pb_MEnergy_1x1v_ser_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_int_mom_1x1v_ser_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void canonical_pb_MEnergy_1x1v_ser_p2(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_int_mom_1x1v_ser_p2(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void canonical_pb_MEnergy_2x2v_ser_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_int_mom_2x2v_ser_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void canonical_pb_MEnergy_2x2v_ser_p2(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_int_mom_2x2v_ser_p2(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void canonical_pb_MEnergy_3x3v_ser_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_int_mom_3x3v_ser_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void canonical_pb_MEnergy_1x1v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_int_mom_1x1v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void canonical_pb_MEnergy_1x1v_tensor_p2(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_int_mom_1x1v_tensor_p2(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void canonical_pb_MEnergy_2x2v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_int_mom_2x2v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void canonical_pb_MEnergy_2x2v_tensor_p2(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_int_mom_2x2v_tensor_p2(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void canonical_pb_MEnergy_3x3v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void canonical_pb_int_mom_3x3v_tensor_p1(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out); 

EXTERN_C_END 
