#pragma once 
#include <math.h> 
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double diffusion_vol_1x_ser_p1(const double* w, const double* dx, const double* D, const double* q, double* GKYL_RESTRICT out);

GKYL_CU_DH void diffusion_surfx_1x_ser_p1(const double* w, const double* dx, const double* Dl, const double* Dc, const double* Dr, const double* ql, const double* qc, const double* qr, double* GKYL_RESTRICT out);

EXTERN_C_END
