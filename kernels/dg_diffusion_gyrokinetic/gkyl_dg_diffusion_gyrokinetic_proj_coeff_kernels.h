#pragma once
#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH void dg_diffusion_gyrokinetic_proj_coeff_1x1v_ser_p1_diffdirsx(const double *xc, const double *dx, const double *nu, const double *xi, double mass, double vtsq_min, const double *gijJ, const double *bmag, const double *vtsq, const double *vmap, const double *vmapSq, double *out);

GKYL_CU_DH void dg_diffusion_gyrokinetic_proj_coeff_1x2v_ser_p1_diffdirsx(const double *xc, const double *dx, const double *nu, const double *xi, double mass, double vtsq_min, const double *gijJ, const double *bmag, const double *vtsq, const double *vmap, const double *vmapSq, double *out);

GKYL_CU_DH void dg_diffusion_gyrokinetic_proj_coeff_2x2v_ser_p1_diffdirsx(const double *xc, const double *dx, const double *nu, const double *xi, double mass, double vtsq_min, const double *gijJ, const double *bmag, const double *vtsq, const double *vmap, const double *vmapSq, double *out);
GKYL_CU_DH void dg_diffusion_gyrokinetic_proj_coeff_2x2v_ser_p1_diffdirsxz(const double *xc, const double *dx, const double *nu, const double *xi, double mass, double vtsq_min, const double *gijJ, const double *bmag, const double *vtsq, const double *vmap, const double *vmapSq, double *out);
GKYL_CU_DH void dg_diffusion_gyrokinetic_proj_coeff_2x2v_ser_p1_diffdirsz(const double *xc, const double *dx, const double *nu, const double *xi, double mass, double vtsq_min, const double *gijJ, const double *bmag, const double *vtsq, const double *vmap, const double *vmapSq, double *out);


EXTERN_C_END
