#pragma once 
#include <math.h> 
#include <stdbool.h> 
#include <gkyl_util.h> 
#include <gkyl_basis_ser_kernels.h>
EXTERN_C_BEG 

GKYL_CU_DH void tok_derived_geo_3x_Ser_p1(const double *gij, const double *bmag, const double *J, double *Jinv, double *grij, double *bi, double *cmag, double *Jtot, double *Jtotinv, double *bmaginv, double *bmaginvsq, double *gxxJ, double *gxyJ, double *gyyJ, double *gxzJ, double *eps2); 
GKYL_CU_DH void adjust_gzz_3x_Ser_p1(const double *cmag, const double *cmag_ref, double *gzz, const double *J, const double *bmag, double *gij); 
GKYL_CU_DH void tok_derived_geo_3x_Ser_p2(const double *gij, const double *bmag, const double *J, double *Jinv, double *grij, double *bi, double *cmag, double *Jtot, double *Jtotinv, double *bmaginv, double *bmaginvsq, double *gxxJ, double *gxyJ, double *gyyJ, double *gxzJ, double *eps2); 
GKYL_CU_DH void adjust_gzz_3x_Ser_p2(const double *cmag, const double *cmag_ref, double *gzz, const double *J, const double *bmag, double *gij); 
EXTERN_C_END 
