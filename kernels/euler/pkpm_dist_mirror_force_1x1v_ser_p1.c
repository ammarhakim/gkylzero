#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void pkpm_dist_mirror_force_1x1v_ser_p1(const double* T_perp_over_m, const double* f, const double* F_k_p_1, double* g_dist_source, double* F_k_m_1) 
{ 
  // T_perp_over_m: Input T_perp/m = p_perp/rho.
  // f:             Input distribution function [F_0, T_perp/m G = T_perp/m (F_0 - F_1)].
  // F_k_p_1:       Input k+1 distribution function. F_2 expansion is the first NP coefficients. 
  // g_dist_source: Output [T_perp/m G, -(T_perp/m)^2 F_1]].
  // F_k_m_1:       Output k-1 distribution function. F_1 expansion is the first NP coefficients. 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[6]; 
  const double *F_2 = &F_k_p_1[0]; 
  double *out_F_0_source = &g_dist_source[0]; 
  double *out_G_1_source = &g_dist_source[6]; 
  double *out_F_1 = &F_k_m_1[0]; 
  double tmp_T_perp_F_0[6] = {0.0}; 

  double tmp_T_perp_F_1[6] = {0.0}; 

  tmp_T_perp_F_0[0] = 0.7071067811865475*F_0[1]*T_perp_over_m[1]+0.7071067811865475*F_0[0]*T_perp_over_m[0]; 
  tmp_T_perp_F_0[1] = 0.7071067811865475*F_0[0]*T_perp_over_m[1]+0.7071067811865475*T_perp_over_m[0]*F_0[1]; 
  tmp_T_perp_F_0[2] = 0.7071067811865475*T_perp_over_m[1]*F_0[3]+0.7071067811865475*T_perp_over_m[0]*F_0[2]; 
  tmp_T_perp_F_0[3] = 0.7071067811865475*T_perp_over_m[0]*F_0[3]+0.7071067811865475*T_perp_over_m[1]*F_0[2]; 
  tmp_T_perp_F_0[4] = 0.7071067811865475*T_perp_over_m[1]*F_0[5]+0.7071067811865475*T_perp_over_m[0]*F_0[4]; 
  tmp_T_perp_F_0[5] = 0.7071067811865475*T_perp_over_m[0]*F_0[5]+0.7071067811865475*T_perp_over_m[1]*F_0[4]; 

  tmp_T_perp_F_1[0] = tmp_T_perp_F_0[0]-1.0*G_1[0]; 
  tmp_T_perp_F_1[1] = tmp_T_perp_F_0[1]-1.0*G_1[1]; 
  tmp_T_perp_F_1[2] = tmp_T_perp_F_0[2]-1.0*G_1[2]; 
  tmp_T_perp_F_1[3] = tmp_T_perp_F_0[3]-1.0*G_1[3]; 
  tmp_T_perp_F_1[4] = tmp_T_perp_F_0[4]-1.0*G_1[4]; 
  tmp_T_perp_F_1[5] = tmp_T_perp_F_0[5]-1.0*G_1[5]; 

  out_F_0_source[0] = G_1[0]; 
  out_F_0_source[1] = G_1[1]; 
  out_F_0_source[2] = G_1[2]; 
  out_F_0_source[3] = G_1[3]; 
  out_F_0_source[4] = G_1[4]; 
  out_F_0_source[5] = G_1[5]; 
  out_G_1_source[0] = (-1.414213562373095*T_perp_over_m[1]*tmp_T_perp_F_1[1])-1.414213562373095*T_perp_over_m[0]*tmp_T_perp_F_1[0]; 
  out_G_1_source[1] = (-1.414213562373095*T_perp_over_m[0]*tmp_T_perp_F_1[1])-1.414213562373095*tmp_T_perp_F_1[0]*T_perp_over_m[1]; 
  out_G_1_source[2] = (-1.414213562373095*T_perp_over_m[1]*tmp_T_perp_F_1[3])-1.414213562373095*T_perp_over_m[0]*tmp_T_perp_F_1[2]; 
  out_G_1_source[3] = (-1.414213562373095*T_perp_over_m[0]*tmp_T_perp_F_1[3])-1.414213562373095*T_perp_over_m[1]*tmp_T_perp_F_1[2]; 
  out_G_1_source[4] = (-1.414213562373095*T_perp_over_m[1]*tmp_T_perp_F_1[5])-1.414213562373095*T_perp_over_m[0]*tmp_T_perp_F_1[4]; 
  out_G_1_source[5] = (-1.414213562373095*T_perp_over_m[0]*tmp_T_perp_F_1[5])-1.414213562373095*T_perp_over_m[1]*tmp_T_perp_F_1[4]; 
} 
