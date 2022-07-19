#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void fpo_vlasov_diff_surfvxvx_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvxvx_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);


EXTERN_C_END 
