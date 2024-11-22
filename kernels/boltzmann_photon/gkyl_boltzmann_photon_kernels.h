#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double boltzmann_photon_vol_1x2v_ser_p1(const double *w, const double *dxv, double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_surfx_1x2v_ser_p1(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_surfvx_1x2v_ser_p1(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_surfvy_1x2v_ser_p1(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double boltzmann_photon_vol_1x2v_ser_p2(const double *w, const double *dxv, double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_surfx_1x2v_ser_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_surfvx_1x2v_ser_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_surfvy_1x2v_ser_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double boltzmann_photon_vol_1x2v_tensor_p2(const double *w, const double *dxv, double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_surfx_1x2v_tensor_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_boundary_surfx_1x2v_tensor_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_surfvx_1x2v_tensor_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_boundary_surfvx_1x2v_tensor_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_surfvy_1x2v_tensor_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double boltzmann_photon_boundary_surfvy_1x2v_tensor_p2(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 

EXTERN_C_END 
