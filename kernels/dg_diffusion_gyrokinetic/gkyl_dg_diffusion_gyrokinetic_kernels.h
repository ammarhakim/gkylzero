#pragma once
#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_1x1v_ser_p1_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_1x1v_ser_p1_varcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_vol_1x1v_ser_p1_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_vol_1x1v_ser_p1_varcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfx_1x1v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x1v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfx_1x1v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x1v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfx_1x1v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_1x1v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfx_1x1v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_1x1v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_1x2v_ser_p1_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_1x2v_ser_p1_varcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_vol_1x2v_ser_p1_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_vol_1x2v_ser_p1_varcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfx_1x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfx_1x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfx_1x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_1x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfx_1x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_1x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_2x2v_ser_p1_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_2x2v_ser_p1_constcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_2x2v_ser_p1_constcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_2x2v_ser_p1_varcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_2x2v_ser_p1_varcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_vol_2x2v_ser_p1_varcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_vol_2x2v_ser_p1_constcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_vol_2x2v_ser_p1_constcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_vol_2x2v_ser_p1_constcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_vol_2x2v_ser_p1_varcoeff_diffdirsx(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_vol_2x2v_ser_p1_varcoeff_diffdirsxz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_vol_2x2v_ser_p1_varcoeff_diffdirsz(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *q, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfx_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfx_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfx_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfx_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfy_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfy_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_surfy_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfy_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfy_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfy_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfy_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);
GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfy_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);


EXTERN_C_END
