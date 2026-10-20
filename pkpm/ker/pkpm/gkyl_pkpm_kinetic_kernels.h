#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double pkpm_kinetic_conf_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *bvar, const double *u, const double *xi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH int pkpm_dist_normalize_1x1v_ser_p1(const double *nmoms, const double *f, double* GKYL_RESTRICT fbar); 
GKYL_CU_DH void pkpm_dist_transverse_1x1v_ser_p1(const double *bvar, double* GKYL_RESTRICT f); 
GKYL_CU_DH double pkpm_conf_flux_x_1x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_stream_surfx_1x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_accel_surfvpar_1x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_vel_flux_vpar_1x1v_ser_p1(const double *w, const double *dxv, const double *nv,
  const double *f_l, const double *f_r, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_kinetic_vel_vol_1x1v_ser_p1(const double *w, const double *dxv,
  const double *c_0, const double *alpha_0, const double *c_h, const double *a_perp,
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_harm_build_1x_ser_p1(double qbym, const double *bvar, const double *grad_b, const double *u,
  const double *grad_u, const double *xi, const double *grad_xi, const double *grad_ln_n, const double *theta,
  const double *dtb, const double *em, double* GKYL_RESTRICT ho); 
GKYL_CU_DH void pkpm_harm_mat_1x_ser_p1(const double *dxv, const double *jacob_pos, const double *ho,
  double* GKYL_RESTRICT mat); 
GKYL_CU_DH void pkpm_harm_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *mat, const double *f,
  double* GKYL_RESTRICT out); 

GKYL_CU_DH double pkpm_kinetic_conf_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *bvar, const double *u, const double *xi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH int pkpm_dist_normalize_2x1v_ser_p1(const double *nmoms, const double *f, double* GKYL_RESTRICT fbar); 
GKYL_CU_DH void pkpm_dist_transverse_2x1v_ser_p1(const double *bvar, double* GKYL_RESTRICT f); 
GKYL_CU_DH double pkpm_conf_flux_x_2x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_stream_surfx_2x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_conf_flux_y_2x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_stream_surfy_2x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_accel_surfvpar_2x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_vel_flux_vpar_2x1v_ser_p1(const double *w, const double *dxv, const double *nv,
  const double *f_l, const double *f_r, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_kinetic_vel_vol_2x1v_ser_p1(const double *w, const double *dxv,
  const double *c_0, const double *alpha_0, const double *c_h, const double *a_perp,
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_harm_build_2x_ser_p1(double qbym, const double *bvar, const double *grad_b, const double *u,
  const double *grad_u, const double *xi, const double *grad_xi, const double *grad_ln_n, const double *theta,
  const double *dtb, const double *em, double* GKYL_RESTRICT ho); 
GKYL_CU_DH void pkpm_harm_mat_2x_ser_p1(const double *dxv, const double *jacob_pos, const double *ho,
  double* GKYL_RESTRICT mat); 
GKYL_CU_DH void pkpm_harm_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *mat, const double *f,
  double* GKYL_RESTRICT out); 

GKYL_CU_DH double pkpm_kinetic_conf_vol_3x1v_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *bvar, const double *u, const double *xi, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH int pkpm_dist_normalize_3x1v_ser_p1(const double *nmoms, const double *f, double* GKYL_RESTRICT fbar); 
GKYL_CU_DH void pkpm_dist_transverse_3x1v_ser_p1(const double *bvar, double* GKYL_RESTRICT f); 
GKYL_CU_DH double pkpm_conf_flux_x_3x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_stream_surfx_3x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_conf_flux_y_3x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_stream_surfy_3x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_conf_flux_z_3x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *b_l, const double *b_r,
  const double *u_l, const double *u_r, const double *xi_l, const double *xi_r,
  const double *theta_l, const double *theta_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_stream_surfz_3x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_accel_surfvpar_3x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH double pkpm_vel_flux_vpar_3x1v_ser_p1(const double *w, const double *dxv, const double *nv,
  const double *f_l, const double *f_r, double* GKYL_RESTRICT flux); 
GKYL_CU_DH double pkpm_kinetic_vel_vol_3x1v_ser_p1(const double *w, const double *dxv,
  const double *c_0, const double *alpha_0, const double *c_h, const double *a_perp,
  const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_harm_build_3x_ser_p1(double qbym, const double *bvar, const double *grad_b, const double *u,
  const double *grad_u, const double *xi, const double *grad_xi, const double *grad_ln_n, const double *theta,
  const double *dtb, const double *em, double* GKYL_RESTRICT ho); 
GKYL_CU_DH void pkpm_harm_mat_3x_ser_p1(const double *dxv, const double *jacob_pos, const double *ho,
  double* GKYL_RESTRICT mat); 
GKYL_CU_DH void pkpm_harm_vol_3x1v_ser_p1(const double *w, const double *dxv, const double *mat, const double *f,
  double* GKYL_RESTRICT out); 

EXTERN_C_END 
