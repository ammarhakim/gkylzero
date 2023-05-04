#pragma once
#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double fpo_vlasov_drag_vol_1x3v_ser_p1(const double* w, const double* dx, const double* h, const double* f, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_vol_1x3v_ser_p1(const double* w, const double* dx, const double* g, const double* f, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_drag_surfvx_1x3v_ser_p1(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvxvx_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvxvy_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvxvz_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_drag_surfvy_1x3v_ser_p1(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvyvx_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvyvy_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvyvz_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_drag_surfvz_1x3v_ser_p1(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvzvx_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvzvy_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvzvz_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_vol_1x3v_ser_p2(const double* w, const double* dx, const double* h, const double* f, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_vol_1x3v_ser_p2(const double* w, const double* dx, const double* g, const double* f, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_drag_surfvx_1x3v_ser_p2(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvxvx_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvxvy_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvxvz_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_drag_surfvy_1x3v_ser_p2(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvyvx_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvyvy_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvyvz_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_drag_surfvz_1x3v_ser_p2(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvzvx_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvzvy_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH void fpo_vlasov_diff_surfvzvz_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);

GKYL_CU_DH void fpo_drag_coeff_recov_surf_vy_3x3v_ser_p1(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vx_3x3v_ser_p1(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vz_2x3v_ser_p2(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vz_2x3v_ser_p1(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vz_3x3v_ser_p1(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vx_2x3v_ser_p1(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vy_2x3v_ser_p1(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vy_2x3v_ser_p2(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vz_1x3v_ser_p1(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vx_2x3v_ser_p2(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vy_1x3v_ser_p2(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vz_1x3v_ser_p2(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vy_1x3v_ser_p1(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vx_1x3v_ser_p1(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vx_1x3v_ser_p2(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff);

GKYL_CU_DH void fpo_drag_coeff_recov_vy_3x3v_ser_p1(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vz_3x3v_ser_p1(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vx_3x3v_ser_p1(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vz_1x3v_ser_p1(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vy_1x3v_ser_p1(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vz_2x3v_ser_p1(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vy_2x3v_ser_p1(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vx_2x3v_ser_p1(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vz_1x3v_ser_p2(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vy_1x3v_ser_p2(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vx_1x3v_ser_p1(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vz_2x3v_ser_p2(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vx_1x3v_ser_p2(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vx_2x3v_ser_p2(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);
GKYL_CU_DH void fpo_drag_coeff_recov_vy_2x3v_ser_p2(const double *dxv, const double *H_l, const double *H_c, const double *H_r, double *drag_coeff);

GKYL_CU_DH void fpo_diff_coeff_recov_surf_vz_3x3v_ser_p1(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vy_3x3v_ser_p1(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vx_3x3v_ser_p1(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vz_2x3v_ser_p2(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vz_2x3v_ser_p1(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vy_2x3v_ser_p1(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vy_2x3v_ser_p2(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vx_2x3v_ser_p1(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vz_1x3v_ser_p2(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vz_1x3v_ser_p1(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vx_2x3v_ser_p2(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vy_1x3v_ser_p2(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vy_1x3v_ser_p1(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vx_1x3v_ser_p1(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vx_1x3v_ser_p2(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff);

GKYL_CU_DH void fpo_diff_coeff_recov_vy_3x3v_ser_p1(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vz_3x3v_ser_p1(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vx_3x3v_ser_p1(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vz_2x3v_ser_p2(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vz_2x3v_ser_p1(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vy_2x3v_ser_p1(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vy_2x3v_ser_p2(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vx_2x3v_ser_p1(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vz_1x3v_ser_p2(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vz_1x3v_ser_p1(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vx_2x3v_ser_p2(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vy_1x3v_ser_p2(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vy_1x3v_ser_p1(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vx_1x3v_ser_p2(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 
GKYL_CU_DH void fpo_diff_coeff_recov_vx_1x3v_ser_p1(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff); 

EXTERN_C_END
