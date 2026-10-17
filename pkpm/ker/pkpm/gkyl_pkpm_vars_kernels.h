#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH int pkpm_vars_theta_1x_ser_p1(const double *vlasov_pkpm_moms, double* GKYL_RESTRICT theta); 
GKYL_CU_DH void pkpm_vars_grad_u_x_1x_ser_p1(const double *dxv, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u); 
GKYL_CU_DH int pkpm_vars_p_force_closure_1x_ser_p1(const double *vlasov_pkpm_moms, const double *rhs_moms, double* GKYL_RESTRICT a_P); 
GKYL_CU_DH void pkpm_vars_pressure_proj_1x_ser_p1(const double *vlasov_pkpm_moms, const double *bvar, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_diag_1x_ser_p1(const double *dxv, const double *vlasov_pkpm_moms, const double *euler_pkpm,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta, const double *a_P,
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 

GKYL_CU_DH int pkpm_vars_theta_2x_ser_p1(const double *vlasov_pkpm_moms, double* GKYL_RESTRICT theta); 
GKYL_CU_DH void pkpm_vars_grad_u_x_2x_ser_p1(const double *dxv, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u); 
GKYL_CU_DH void pkpm_vars_grad_u_y_2x_ser_p1(const double *dxv, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u); 
GKYL_CU_DH int pkpm_vars_p_force_closure_2x_ser_p1(const double *vlasov_pkpm_moms, const double *rhs_moms, double* GKYL_RESTRICT a_P); 
GKYL_CU_DH void pkpm_vars_pressure_proj_2x_ser_p1(const double *vlasov_pkpm_moms, const double *bvar, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_diag_2x_ser_p1(const double *dxv, const double *vlasov_pkpm_moms, const double *euler_pkpm,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta, const double *a_P,
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 

GKYL_CU_DH int pkpm_vars_theta_3x_ser_p1(const double *vlasov_pkpm_moms, double* GKYL_RESTRICT theta); 
GKYL_CU_DH void pkpm_vars_grad_u_x_3x_ser_p1(const double *dxv, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u); 
GKYL_CU_DH void pkpm_vars_grad_u_y_3x_ser_p1(const double *dxv, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u); 
GKYL_CU_DH void pkpm_vars_grad_u_z_3x_ser_p1(const double *dxv, const double *u_l, const double *u_c, const double *u_r, double* GKYL_RESTRICT grad_u); 
GKYL_CU_DH int pkpm_vars_p_force_closure_3x_ser_p1(const double *vlasov_pkpm_moms, const double *rhs_moms, double* GKYL_RESTRICT a_P); 
GKYL_CU_DH void pkpm_vars_pressure_proj_3x_ser_p1(const double *vlasov_pkpm_moms, const double *bvar, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH void pkpm_vars_diag_3x_ser_p1(const double *dxv, const double *vlasov_pkpm_moms, const double *euler_pkpm,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta, const double *a_P,
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 

EXTERN_C_END 
