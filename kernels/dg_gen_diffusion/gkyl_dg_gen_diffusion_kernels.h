#pragma once 
#include <math.h> 
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH void dg_gen_diffusion_surfx_2x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfy_2x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);

EXTERN_C_END
