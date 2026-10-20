#pragma once 
#include <math.h> 
#include <stdbool.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double pkpm_fluid_vol_1x_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *vlasov_pkpm_moms, const double *bvar, const double *u, const double *fluid,
  double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_flux_accum_x_1x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT I); 
GKYL_CU_DH void pkpm_fluid_flux_x_1x_ser_p1(const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *fluid_l, const double *fluid_r, const double *I, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_fluid_surfx_1x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH double pkpm_fluid_vol_2x_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *vlasov_pkpm_moms, const double *bvar, const double *u, const double *fluid,
  double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_flux_accum_x_2x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT I); 
GKYL_CU_DH void pkpm_fluid_flux_x_2x_ser_p1(const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *fluid_l, const double *fluid_r, const double *I, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_fluid_surfx_2x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_flux_accum_y_2x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT I); 
GKYL_CU_DH void pkpm_fluid_flux_y_2x_ser_p1(const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *fluid_l, const double *fluid_r, const double *I, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_fluid_surfy_2x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

GKYL_CU_DH double pkpm_fluid_vol_3x_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *vlasov_pkpm_moms, const double *bvar, const double *u, const double *fluid,
  double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_flux_accum_x_3x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT I); 
GKYL_CU_DH void pkpm_fluid_flux_x_3x_ser_p1(const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *fluid_l, const double *fluid_r, const double *I, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_fluid_surfx_3x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_flux_accum_y_3x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT I); 
GKYL_CU_DH void pkpm_fluid_flux_y_3x_ser_p1(const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *fluid_l, const double *fluid_r, const double *I, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_fluid_surfy_3x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_fluid_flux_accum_z_3x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT I); 
GKYL_CU_DH void pkpm_fluid_flux_z_3x_ser_p1(const double *bvar_l, const double *bvar_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *fluid_l, const double *fluid_r, const double *I, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_fluid_surfz_3x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 

EXTERN_C_END 
