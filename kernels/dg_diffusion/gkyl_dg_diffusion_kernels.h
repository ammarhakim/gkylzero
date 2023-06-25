#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double dg_diffusion_vol_1x_ser_p1(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion4_vol_1x_ser_p1(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion6_vol_1x_ser_p1(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_diffusion_vol_1x_ser_p2(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion4_vol_1x_ser_p2(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion6_vol_1x_ser_p2(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_diffusion_vol_2x_ser_p1(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion4_vol_2x_ser_p1(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion6_vol_2x_ser_p1(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_diffusion_vol_2x_ser_p2(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion4_vol_2x_ser_p2(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion6_vol_2x_ser_p2(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_diffusion_vol_3x_ser_p1(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion4_vol_3x_ser_p1(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion6_vol_3x_ser_p1(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 

GKYL_CU_DH double dg_diffusion_vol_3x_ser_p2(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion4_vol_3x_ser_p2(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double dg_diffusion6_vol_3x_ser_p2(const double *w, const double *dxv, double D, const double *q, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfx_1x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfx_1x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfx_1x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfx_1x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfx_1x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfx_1x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfx_1x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfx_1x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfx_1x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfx_1x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfx_1x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfx_1x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfx_1x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfx_1x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfx_1x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfx_1x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfx_1x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfx_1x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfx_1x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfx_1x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfx_1x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfx_1x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfx_1x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfx_1x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfx_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfx_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfx_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfx_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfx_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfx_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfx_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfx_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfx_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfx_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfx_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfx_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfy_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfy_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfy_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfy_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfy_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfy_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfy_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfy_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfy_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfy_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfy_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfy_2x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfx_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfx_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfx_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfx_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfx_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfx_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfx_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfx_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfx_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfx_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfx_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfx_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfy_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfy_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfy_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfy_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfy_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfy_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfy_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfy_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfy_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfy_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfy_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfy_2x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfx_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfx_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfx_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfx_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfx_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfx_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfx_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfx_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfx_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfx_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfx_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfx_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfy_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfy_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfy_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfy_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfy_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfy_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfy_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfy_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfy_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfy_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfy_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfy_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfz_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfz_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfz_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfz_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfz_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfz_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfz_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfz_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfz_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfz_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfz_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfz_3x_ser_p1(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfx_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfx_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfx_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfx_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfx_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfx_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfx_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfx_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfx_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfx_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfx_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfx_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfy_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfy_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfy_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfy_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfy_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfy_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfy_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfy_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfy_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfy_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfy_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfy_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfz_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfz_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfz_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfz_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfz_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfz_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfz_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfz_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfz_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfz_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfz_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfz_3x_ser_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfx_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfx_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfx_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfx_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfx_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfx_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfx_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfx_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfx_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfx_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfx_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfx_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfy_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfy_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfy_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfy_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfy_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfy_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfy_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfy_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfy_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfy_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfy_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfy_2x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfx_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfx_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfx_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfx_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfx_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfx_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfx_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfx_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfx_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfx_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfx_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfx_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfy_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfy_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfy_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfy_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfy_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfy_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfy_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfy_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfy_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfy_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfy_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfy_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_surfz_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_surfz_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_surfz_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_pkpm_surfz_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_pkpm_surfz_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_pkpm_surfz_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_iso_euler_surfz_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_iso_euler_surfz_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_iso_euler_surfz_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void dg_diffusion_euler_surfz_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion4_euler_surfz_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void dg_diffusion6_euler_surfz_3x_tensor_p2(const double *w, const double *dxv, 
      double D, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out); 

EXTERN_C_END 
