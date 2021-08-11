#pragma once 
#include <math.h> 
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH void lbo_constNu_surf_1x1v_ser_vx_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out);
GKYL_CU_DH void lbo_constNu_surf_1x1v_ser_vx_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out);
GKYL_CU_DH void lbo_constNu_surf_1x2v_ser_vx_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out);
GKYL_CU_DH void lbo_constNu_surf_1x2v_ser_vy_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out);
GKYL_CU_DH void lbo_constNu_surf_1x2v_ser_vx_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out);
GKYL_CU_DH void lbo_constNu_surf_1x2v_ser_vy_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out);
GKYL_CU_DH void lbo_constNu_surf_1x3v_ser_vx_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out);
GKYL_CU_DH void lbo_constNu_surf_1x3v_ser_vy_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out);
GKYL_CU_DH void lbo_constNu_surf_1x3v_ser_vz_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out);
GKYL_CU_DH void lbo_constNu_surf_1x3v_ser_vx_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out);
GKYL_CU_DH void lbo_constNu_surf_1x3v_ser_vy_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out);
GKYL_CU_DH void lbo_constNu_surf_1x3v_ser_vz_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out);

EXTERN_C_END
