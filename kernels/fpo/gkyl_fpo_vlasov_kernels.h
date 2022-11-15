#pragma once
#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double fpo_vlasov_drag_vol_1x3v_ser_p1(const double* w, const double* dx, const double* a, const double* f, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_vol_1x3v_ser_p1(const double* w, const double* dx, const double* D, const double* f, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_drag_surfvx_1x3v_ser_p1(const double* w, const double* dx, const double* al, const double* ac, const double* au, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvxvx_1x3v_ser_p1(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvxvy_1x3v_ser_p1(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvxvz_1x3v_ser_p1(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_Dvxvx_surf_1x3v_ser_p1(const double* w, const double* dx, const double* g[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_Dvxvy_surf_1x3v_ser_p1(const double* w, const double* dx, const double* g[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_Dvxvz_surf_1x3v_ser_p1(const double* w, const double* dx, const double* g[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_drag_surfvy_1x3v_ser_p1(const double* w, const double* dx, const double* al, const double* ac, const double* au, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvyvx_1x3v_ser_p1(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvyvy_1x3v_ser_p1(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvyvz_1x3v_ser_p1(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_Dvyvy_surf_1x3v_ser_p1(const double* w, const double* dx, const double* g[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_Dvyvz_surf_1x3v_ser_p1(const double* w, const double* dx, const double* g[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_drag_surfvz_1x3v_ser_p1(const double* w, const double* dx, const double* al, const double* ac, const double* au, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvzvx_1x3v_ser_p1(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvzvy_1x3v_ser_p1(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvzvz_1x3v_ser_p1(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_Dvzvz_surf_1x3v_ser_p1(const double* w, const double* dx, const double* g[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_vol_1x3v_ser_p2(const double* w, const double* dx, const double* a, const double* f, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_vol_1x3v_ser_p2(const double* w, const double* dx, const double* D, const double* f, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_drag_surfvx_1x3v_ser_p2(const double* w, const double* dx, const double* al, const double* ac, const double* au, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvxvx_1x3v_ser_p2(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvxvy_1x3v_ser_p2(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvxvz_1x3v_ser_p2(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_Dvxvx_surf_1x3v_ser_p2(const double* w, const double* dx, const double* g[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_Dvxvy_surf_1x3v_ser_p2(const double* w, const double* dx, const double* g[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_Dvxvz_surf_1x3v_ser_p2(const double* w, const double* dx, const double* g[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_drag_surfvy_1x3v_ser_p2(const double* w, const double* dx, const double* al, const double* ac, const double* au, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvyvx_1x3v_ser_p2(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvyvy_1x3v_ser_p2(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvyvz_1x3v_ser_p2(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_Dvyvy_surf_1x3v_ser_p2(const double* w, const double* dx, const double* g[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_Dvyvz_surf_1x3v_ser_p2(const double* w, const double* dx, const double* g[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_drag_surfvz_1x3v_ser_p2(const double* w, const double* dx, const double* al, const double* ac, const double* au, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvzvx_1x3v_ser_p2(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvzvy_1x3v_ser_p2(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvzvz_1x3v_ser_p2(const double* w, const double* dx, const double* D[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_Dvzvz_surf_1x3v_ser_p2(const double* w, const double* dx, const double* g[], double* GKYL_RESTRICT out);


EXTERN_C_END