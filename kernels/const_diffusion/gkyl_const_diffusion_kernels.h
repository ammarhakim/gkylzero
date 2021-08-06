#pragma once 
#include <math.h> 
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double const_diffusion_vol_1x_ser_p1(const double* w, const double* dx, const double* D, const double* f, double* GKYL_RESTRICT out);
GKYL_CU_DH double const_diffusion_vol_1x_ser_p2(const double* w, const double* dx, const double* D, const double* f, double* GKYL_RESTRICT out);

GKYL_CU_DH void const_diffusion_surfx_1x_ser_p1(const double* w, const double* dx, const double* D, const double* fl, const double* fc, const double* fr, double* GKYL_RESTRICT out);
GKYL_CU_DH void const_diffusion_surfx_1x_ser_p2(const double* w, const double* dx, const double* D, const double* fl, const double* fc, const double* fr, double* GKYL_RESTRICT out);

EXTERN_C_END
