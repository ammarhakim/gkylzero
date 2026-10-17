#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double pkpm_kinetic_conf_vol_1x1v_ser_p1(const double *w, const double *dxv,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta,
  const double *nu_sum, const double *nu_vthsq_sum,
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_conf_flux_x_1x1v_ser_p1(const double *w, const double *dxv,
  const double *b_l, const double *b_r, const double *u_l, const double *u_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux); 
GKYL_CU_DH void pkpm_conf_flux_corr_x_1x1v_ser_p1(double mass, const double *moms_l, const double *moms_r,
  const double *roe_mass_flux, const double *kin_mass_flux, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_stream_surfx_1x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_accel_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_accel_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_vel_flux_vpar_1x1v_ser_p1(const double *w, const double *dxv,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta,
  const double *a_P, const double *f_l, const double *f_r, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_kinetic_vel_vol_1x1v_ser_p1(const double *w, const double *dxv,
  const double *a_P, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH double pkpm_kinetic_conf_vol_2x1v_ser_p1(const double *w, const double *dxv,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta,
  const double *nu_sum, const double *nu_vthsq_sum,
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_conf_flux_x_2x1v_ser_p1(const double *w, const double *dxv,
  const double *b_l, const double *b_r, const double *u_l, const double *u_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux); 
GKYL_CU_DH void pkpm_conf_flux_corr_x_2x1v_ser_p1(double mass, const double *moms_l, const double *moms_r,
  const double *roe_mass_flux, const double *kin_mass_flux, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_stream_surfx_2x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_conf_flux_y_2x1v_ser_p1(const double *w, const double *dxv,
  const double *b_l, const double *b_r, const double *u_l, const double *u_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux); 
GKYL_CU_DH void pkpm_conf_flux_corr_y_2x1v_ser_p1(double mass, const double *moms_l, const double *moms_r,
  const double *roe_mass_flux, const double *kin_mass_flux, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_stream_surfy_2x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_accel_surfvpar_2x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_accel_boundary_surfvpar_2x1v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_vel_flux_vpar_2x1v_ser_p1(const double *w, const double *dxv,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta,
  const double *a_P, const double *f_l, const double *f_r, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_kinetic_vel_vol_2x1v_ser_p1(const double *w, const double *dxv,
  const double *a_P, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH double pkpm_kinetic_conf_vol_3x1v_ser_p1(const double *w, const double *dxv,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta,
  const double *nu_sum, const double *nu_vthsq_sum,
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_conf_flux_x_3x1v_ser_p1(const double *w, const double *dxv,
  const double *b_l, const double *b_r, const double *u_l, const double *u_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux); 
GKYL_CU_DH void pkpm_conf_flux_corr_x_3x1v_ser_p1(double mass, const double *moms_l, const double *moms_r,
  const double *roe_mass_flux, const double *kin_mass_flux, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_stream_surfx_3x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_conf_flux_y_3x1v_ser_p1(const double *w, const double *dxv,
  const double *b_l, const double *b_r, const double *u_l, const double *u_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux); 
GKYL_CU_DH void pkpm_conf_flux_corr_y_3x1v_ser_p1(double mass, const double *moms_l, const double *moms_r,
  const double *roe_mass_flux, const double *kin_mass_flux, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_stream_surfy_3x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_conf_flux_z_3x1v_ser_p1(const double *w, const double *dxv,
  const double *b_l, const double *b_r, const double *u_l, const double *u_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux, double* GKYL_RESTRICT mass_flux); 
GKYL_CU_DH void pkpm_conf_flux_corr_z_3x1v_ser_p1(double mass, const double *moms_l, const double *moms_r,
  const double *roe_mass_flux, const double *kin_mass_flux, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_stream_surfz_3x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_accel_surfvpar_3x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_accel_boundary_surfvpar_3x1v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_vel_flux_vpar_3x1v_ser_p1(const double *w, const double *dxv,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta,
  const double *a_P, const double *f_l, const double *f_r, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_kinetic_vel_vol_3x1v_ser_p1(const double *w, const double *dxv,
  const double *a_P, const double *f, double* GKYL_RESTRICT out); 

EXTERN_C_END 
