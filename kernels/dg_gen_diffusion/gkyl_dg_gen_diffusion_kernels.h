#pragma once 
#include <math.h> 
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double dg_gen_diffusion_vol_2x_ser_p1(const double* w, const double* dx, const double* D, const double* q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_gen_diffusion_vol_2x_ser_p2(const double* w, const double* dx, const double* D, const double* q, double* GKYL_RESTRICT out);

GKYL_CU_DH void dg_gen_diffusion_surfxx_2x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfxx_2x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfxy_2x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfxy_2x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfyx_2x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfyx_2x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfyy_2x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfyy_2x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_gen_diffusion_vol_3x_ser_p1(const double* w, const double* dx, const double* D, const double* q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_gen_diffusion_vol_3x_ser_p2(const double* w, const double* dx, const double* D, const double* q, double* GKYL_RESTRICT out);

GKYL_CU_DH void dg_gen_diffusion_surfxx_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfxx_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfxy_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfxy_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfxz_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfxz_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);

GKYL_CU_DH void dg_gen_diffusion_surfyx_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfyx_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfyy_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfyy_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfyz_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfyz_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);

GKYL_CU_DH void dg_gen_diffusion_surfzx_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfzx_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfzy_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfzy_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfzz_3x_ser_p2(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);
GKYL_CU_DH void dg_gen_diffusion_surfzz_3x_ser_p1(const double* w, const double* dx, const double* Dij, const double* q[], double* GKYL_RESTRICT out);

EXTERN_C_END