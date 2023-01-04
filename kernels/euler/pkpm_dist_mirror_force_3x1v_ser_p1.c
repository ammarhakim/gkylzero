#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
GKYL_CU_DH void pkpm_dist_mirror_force_3x1v_ser_p1(const double* T_perp_over_m, const double* f, const double* F_k_p_1, double* g_dist_source, double* F_k_m_1) 
{ 
  // T_perp_over_m: Input T_perp/m = p_perp/rho.
  // f:             Input distribution function [F_0, T_perp/m G = T_perp/m (F_0 - F_1)].
  // F_k_p_1:       Input k+1 distribution function. F_2 expansion is the first NP coefficients. 
  // g_dist_source: Output [T_perp/m G, -(T_perp/m)^2 F_1]].
  // F_k_m_1:       Output k-1 distribution function. F_1 expansion is the first NP coefficients. 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[24]; 
  const double *F_2 = &F_k_p_1[0]; 
  double *out_F_0_source = &g_dist_source[0]; 
  double *out_G_1_source = &g_dist_source[24]; 
  double *out_F_1 = &F_k_m_1[0]; 
  double tmp_T_perp_F_0[24] = {0.0}; 

  double tmp_T_perp_F_1[24] = {0.0}; 

  tmp_T_perp_F_0[0] = 0.3535533905932737*T_perp_over_m[7]*F_0[11]+0.3535533905932737*T_perp_over_m[6]*F_0[7]+0.3535533905932737*T_perp_over_m[5]*F_0[6]+0.3535533905932737*T_perp_over_m[4]*F_0[5]+0.3535533905932737*F_0[3]*T_perp_over_m[3]+0.3535533905932737*F_0[2]*T_perp_over_m[2]+0.3535533905932737*F_0[1]*T_perp_over_m[1]+0.3535533905932737*F_0[0]*T_perp_over_m[0]; 
  tmp_T_perp_F_0[1] = 0.3535533905932737*T_perp_over_m[6]*F_0[11]+0.3535533905932737*F_0[7]*T_perp_over_m[7]+0.3535533905932737*T_perp_over_m[3]*F_0[6]+0.3535533905932737*F_0[3]*T_perp_over_m[5]+0.3535533905932737*T_perp_over_m[2]*F_0[5]+0.3535533905932737*F_0[2]*T_perp_over_m[4]+0.3535533905932737*F_0[0]*T_perp_over_m[1]+0.3535533905932737*T_perp_over_m[0]*F_0[1]; 
  tmp_T_perp_F_0[2] = 0.3535533905932737*T_perp_over_m[5]*F_0[11]+0.3535533905932737*F_0[6]*T_perp_over_m[7]+0.3535533905932737*T_perp_over_m[3]*F_0[7]+0.3535533905932737*F_0[3]*T_perp_over_m[6]+0.3535533905932737*T_perp_over_m[1]*F_0[5]+0.3535533905932737*F_0[1]*T_perp_over_m[4]+0.3535533905932737*F_0[0]*T_perp_over_m[2]+0.3535533905932737*T_perp_over_m[0]*F_0[2]; 
  tmp_T_perp_F_0[3] = 0.3535533905932737*T_perp_over_m[4]*F_0[11]+0.3535533905932737*F_0[5]*T_perp_over_m[7]+0.3535533905932737*T_perp_over_m[2]*F_0[7]+0.3535533905932737*F_0[2]*T_perp_over_m[6]+0.3535533905932737*T_perp_over_m[1]*F_0[6]+0.3535533905932737*F_0[1]*T_perp_over_m[5]+0.3535533905932737*F_0[0]*T_perp_over_m[3]+0.3535533905932737*T_perp_over_m[0]*F_0[3]; 
  tmp_T_perp_F_0[4] = 0.3535533905932737*T_perp_over_m[7]*F_0[15]+0.3535533905932737*T_perp_over_m[6]*F_0[14]+0.3535533905932737*T_perp_over_m[5]*F_0[13]+0.3535533905932737*T_perp_over_m[4]*F_0[12]+0.3535533905932737*T_perp_over_m[3]*F_0[10]+0.3535533905932737*T_perp_over_m[2]*F_0[9]+0.3535533905932737*T_perp_over_m[1]*F_0[8]+0.3535533905932737*T_perp_over_m[0]*F_0[4]; 
  tmp_T_perp_F_0[5] = 0.3535533905932737*T_perp_over_m[3]*F_0[11]+0.3535533905932737*F_0[3]*T_perp_over_m[7]+0.3535533905932737*T_perp_over_m[5]*F_0[7]+0.3535533905932737*F_0[6]*T_perp_over_m[6]+0.3535533905932737*T_perp_over_m[0]*F_0[5]+0.3535533905932737*F_0[0]*T_perp_over_m[4]+0.3535533905932737*F_0[1]*T_perp_over_m[2]+0.3535533905932737*T_perp_over_m[1]*F_0[2]; 
  tmp_T_perp_F_0[6] = 0.3535533905932737*T_perp_over_m[2]*F_0[11]+0.3535533905932737*F_0[2]*T_perp_over_m[7]+0.3535533905932737*T_perp_over_m[4]*F_0[7]+0.3535533905932737*F_0[5]*T_perp_over_m[6]+0.3535533905932737*T_perp_over_m[0]*F_0[6]+0.3535533905932737*F_0[0]*T_perp_over_m[5]+0.3535533905932737*F_0[1]*T_perp_over_m[3]+0.3535533905932737*T_perp_over_m[1]*F_0[3]; 
  tmp_T_perp_F_0[7] = 0.3535533905932737*T_perp_over_m[1]*F_0[11]+0.3535533905932737*F_0[1]*T_perp_over_m[7]+0.3535533905932737*T_perp_over_m[0]*F_0[7]+0.3535533905932737*F_0[0]*T_perp_over_m[6]+0.3535533905932737*T_perp_over_m[4]*F_0[6]+0.3535533905932737*F_0[5]*T_perp_over_m[5]+0.3535533905932737*F_0[2]*T_perp_over_m[3]+0.3535533905932737*T_perp_over_m[2]*F_0[3]; 
  tmp_T_perp_F_0[8] = 0.3535533905932737*T_perp_over_m[6]*F_0[15]+0.3535533905932737*T_perp_over_m[7]*F_0[14]+0.3535533905932737*T_perp_over_m[3]*F_0[13]+0.3535533905932737*T_perp_over_m[2]*F_0[12]+0.3535533905932737*T_perp_over_m[5]*F_0[10]+0.3535533905932737*T_perp_over_m[4]*F_0[9]+0.3535533905932737*T_perp_over_m[0]*F_0[8]+0.3535533905932737*T_perp_over_m[1]*F_0[4]; 
  tmp_T_perp_F_0[9] = 0.3535533905932737*T_perp_over_m[5]*F_0[15]+0.3535533905932737*T_perp_over_m[3]*F_0[14]+0.3535533905932737*T_perp_over_m[7]*F_0[13]+0.3535533905932737*T_perp_over_m[1]*F_0[12]+0.3535533905932737*T_perp_over_m[6]*F_0[10]+0.3535533905932737*T_perp_over_m[0]*F_0[9]+0.3535533905932737*T_perp_over_m[4]*F_0[8]+0.3535533905932737*T_perp_over_m[2]*F_0[4]; 
  tmp_T_perp_F_0[10] = 0.3535533905932737*T_perp_over_m[4]*F_0[15]+0.3535533905932737*T_perp_over_m[2]*F_0[14]+0.3535533905932737*T_perp_over_m[1]*F_0[13]+0.3535533905932737*T_perp_over_m[7]*F_0[12]+0.3535533905932737*T_perp_over_m[0]*F_0[10]+0.3535533905932737*T_perp_over_m[6]*F_0[9]+0.3535533905932737*T_perp_over_m[5]*F_0[8]+0.3535533905932737*T_perp_over_m[3]*F_0[4]; 
  tmp_T_perp_F_0[11] = 0.3535533905932737*T_perp_over_m[0]*F_0[11]+0.3535533905932737*F_0[0]*T_perp_over_m[7]+0.3535533905932737*T_perp_over_m[1]*F_0[7]+0.3535533905932737*F_0[1]*T_perp_over_m[6]+0.3535533905932737*T_perp_over_m[2]*F_0[6]+0.3535533905932737*F_0[2]*T_perp_over_m[5]+0.3535533905932737*T_perp_over_m[3]*F_0[5]+0.3535533905932737*F_0[3]*T_perp_over_m[4]; 
  tmp_T_perp_F_0[12] = 0.3535533905932737*T_perp_over_m[3]*F_0[15]+0.3535533905932737*T_perp_over_m[5]*F_0[14]+0.3535533905932737*T_perp_over_m[6]*F_0[13]+0.3535533905932737*T_perp_over_m[0]*F_0[12]+0.3535533905932737*T_perp_over_m[7]*F_0[10]+0.3535533905932737*T_perp_over_m[1]*F_0[9]+0.3535533905932737*T_perp_over_m[2]*F_0[8]+0.3535533905932737*F_0[4]*T_perp_over_m[4]; 
  tmp_T_perp_F_0[13] = 0.3535533905932737*T_perp_over_m[2]*F_0[15]+0.3535533905932737*T_perp_over_m[4]*F_0[14]+0.3535533905932737*T_perp_over_m[0]*F_0[13]+0.3535533905932737*T_perp_over_m[6]*F_0[12]+0.3535533905932737*T_perp_over_m[1]*F_0[10]+0.3535533905932737*T_perp_over_m[7]*F_0[9]+0.3535533905932737*T_perp_over_m[3]*F_0[8]+0.3535533905932737*F_0[4]*T_perp_over_m[5]; 
  tmp_T_perp_F_0[14] = 0.3535533905932737*T_perp_over_m[1]*F_0[15]+0.3535533905932737*T_perp_over_m[0]*F_0[14]+0.3535533905932737*T_perp_over_m[4]*F_0[13]+0.3535533905932737*T_perp_over_m[5]*F_0[12]+0.3535533905932737*T_perp_over_m[2]*F_0[10]+0.3535533905932737*T_perp_over_m[3]*F_0[9]+0.3535533905932737*T_perp_over_m[7]*F_0[8]+0.3535533905932737*F_0[4]*T_perp_over_m[6]; 
  tmp_T_perp_F_0[15] = 0.3535533905932737*T_perp_over_m[0]*F_0[15]+0.3535533905932737*T_perp_over_m[1]*F_0[14]+0.3535533905932737*T_perp_over_m[2]*F_0[13]+0.3535533905932737*T_perp_over_m[3]*F_0[12]+0.3535533905932737*T_perp_over_m[4]*F_0[10]+0.3535533905932737*T_perp_over_m[5]*F_0[9]+0.3535533905932737*T_perp_over_m[6]*F_0[8]+0.3535533905932737*F_0[4]*T_perp_over_m[7]; 
  tmp_T_perp_F_0[16] = 0.3535533905932737*T_perp_over_m[7]*F_0[23]+0.3535533905932737*T_perp_over_m[6]*F_0[22]+0.3535533905932737*T_perp_over_m[5]*F_0[21]+0.3535533905932737*T_perp_over_m[4]*F_0[20]+0.3535533905932737*T_perp_over_m[3]*F_0[19]+0.3535533905932737*T_perp_over_m[2]*F_0[18]+0.3535533905932737*T_perp_over_m[1]*F_0[17]+0.3535533905932737*T_perp_over_m[0]*F_0[16]; 
  tmp_T_perp_F_0[17] = 0.3535533905932737*T_perp_over_m[6]*F_0[23]+0.3535533905932737*T_perp_over_m[7]*F_0[22]+0.3535533905932737*T_perp_over_m[3]*F_0[21]+0.3535533905932737*T_perp_over_m[2]*F_0[20]+0.3535533905932737*T_perp_over_m[5]*F_0[19]+0.3535533905932737*T_perp_over_m[4]*F_0[18]+0.3535533905932737*T_perp_over_m[0]*F_0[17]+0.3535533905932737*T_perp_over_m[1]*F_0[16]; 
  tmp_T_perp_F_0[18] = 0.3535533905932737*T_perp_over_m[5]*F_0[23]+0.3535533905932737*T_perp_over_m[3]*F_0[22]+0.3535533905932737*T_perp_over_m[7]*F_0[21]+0.3535533905932737*T_perp_over_m[1]*F_0[20]+0.3535533905932737*T_perp_over_m[6]*F_0[19]+0.3535533905932737*T_perp_over_m[0]*F_0[18]+0.3535533905932737*T_perp_over_m[4]*F_0[17]+0.3535533905932737*T_perp_over_m[2]*F_0[16]; 
  tmp_T_perp_F_0[19] = 0.3535533905932737*T_perp_over_m[4]*F_0[23]+0.3535533905932737*T_perp_over_m[2]*F_0[22]+0.3535533905932737*T_perp_over_m[1]*F_0[21]+0.3535533905932737*T_perp_over_m[7]*F_0[20]+0.3535533905932737*T_perp_over_m[0]*F_0[19]+0.3535533905932737*T_perp_over_m[6]*F_0[18]+0.3535533905932737*T_perp_over_m[5]*F_0[17]+0.3535533905932737*T_perp_over_m[3]*F_0[16]; 
  tmp_T_perp_F_0[20] = 0.3535533905932737*T_perp_over_m[3]*F_0[23]+0.3535533905932737*T_perp_over_m[5]*F_0[22]+0.3535533905932737*T_perp_over_m[6]*F_0[21]+0.3535533905932737*T_perp_over_m[0]*F_0[20]+0.3535533905932737*T_perp_over_m[7]*F_0[19]+0.3535533905932737*T_perp_over_m[1]*F_0[18]+0.3535533905932737*T_perp_over_m[2]*F_0[17]+0.3535533905932737*T_perp_over_m[4]*F_0[16]; 
  tmp_T_perp_F_0[21] = 0.3535533905932737*T_perp_over_m[2]*F_0[23]+0.3535533905932737*T_perp_over_m[4]*F_0[22]+0.3535533905932737*T_perp_over_m[0]*F_0[21]+0.3535533905932737*T_perp_over_m[6]*F_0[20]+0.3535533905932737*T_perp_over_m[1]*F_0[19]+0.3535533905932737*T_perp_over_m[7]*F_0[18]+0.3535533905932737*T_perp_over_m[3]*F_0[17]+0.3535533905932737*T_perp_over_m[5]*F_0[16]; 
  tmp_T_perp_F_0[22] = 0.3535533905932737*T_perp_over_m[1]*F_0[23]+0.3535533905932737*T_perp_over_m[0]*F_0[22]+0.3535533905932737*T_perp_over_m[4]*F_0[21]+0.3535533905932737*T_perp_over_m[5]*F_0[20]+0.3535533905932737*T_perp_over_m[2]*F_0[19]+0.3535533905932737*T_perp_over_m[3]*F_0[18]+0.3535533905932737*T_perp_over_m[7]*F_0[17]+0.3535533905932737*T_perp_over_m[6]*F_0[16]; 
  tmp_T_perp_F_0[23] = 0.3535533905932737*T_perp_over_m[0]*F_0[23]+0.3535533905932737*T_perp_over_m[1]*F_0[22]+0.3535533905932737*T_perp_over_m[2]*F_0[21]+0.3535533905932737*T_perp_over_m[3]*F_0[20]+0.3535533905932737*T_perp_over_m[4]*F_0[19]+0.3535533905932737*T_perp_over_m[5]*F_0[18]+0.3535533905932737*T_perp_over_m[6]*F_0[17]+0.3535533905932737*T_perp_over_m[7]*F_0[16]; 

  tmp_T_perp_F_1[0] = tmp_T_perp_F_0[0]-1.0*G_1[0]; 
  tmp_T_perp_F_1[1] = tmp_T_perp_F_0[1]-1.0*G_1[1]; 
  tmp_T_perp_F_1[2] = tmp_T_perp_F_0[2]-1.0*G_1[2]; 
  tmp_T_perp_F_1[3] = tmp_T_perp_F_0[3]-1.0*G_1[3]; 
  tmp_T_perp_F_1[4] = tmp_T_perp_F_0[4]-1.0*G_1[4]; 
  tmp_T_perp_F_1[5] = tmp_T_perp_F_0[5]-1.0*G_1[5]; 
  tmp_T_perp_F_1[6] = tmp_T_perp_F_0[6]-1.0*G_1[6]; 
  tmp_T_perp_F_1[7] = tmp_T_perp_F_0[7]-1.0*G_1[7]; 
  tmp_T_perp_F_1[8] = tmp_T_perp_F_0[8]-1.0*G_1[8]; 
  tmp_T_perp_F_1[9] = tmp_T_perp_F_0[9]-1.0*G_1[9]; 
  tmp_T_perp_F_1[10] = tmp_T_perp_F_0[10]-1.0*G_1[10]; 
  tmp_T_perp_F_1[11] = tmp_T_perp_F_0[11]-1.0*G_1[11]; 
  tmp_T_perp_F_1[12] = tmp_T_perp_F_0[12]-1.0*G_1[12]; 
  tmp_T_perp_F_1[13] = tmp_T_perp_F_0[13]-1.0*G_1[13]; 
  tmp_T_perp_F_1[14] = tmp_T_perp_F_0[14]-1.0*G_1[14]; 
  tmp_T_perp_F_1[15] = tmp_T_perp_F_0[15]-1.0*G_1[15]; 
  tmp_T_perp_F_1[16] = tmp_T_perp_F_0[16]-1.0*G_1[16]; 
  tmp_T_perp_F_1[17] = tmp_T_perp_F_0[17]-1.0*G_1[17]; 
  tmp_T_perp_F_1[18] = tmp_T_perp_F_0[18]-1.0*G_1[18]; 
  tmp_T_perp_F_1[19] = tmp_T_perp_F_0[19]-1.0*G_1[19]; 
  tmp_T_perp_F_1[20] = tmp_T_perp_F_0[20]-1.0*G_1[20]; 
  tmp_T_perp_F_1[21] = tmp_T_perp_F_0[21]-1.0*G_1[21]; 
  tmp_T_perp_F_1[22] = tmp_T_perp_F_0[22]-1.0*G_1[22]; 
  tmp_T_perp_F_1[23] = tmp_T_perp_F_0[23]-1.0*G_1[23]; 

  out_F_0_source[0] = G_1[0]; 
  out_F_0_source[1] = G_1[1]; 
  out_F_0_source[2] = G_1[2]; 
  out_F_0_source[3] = G_1[3]; 
  out_F_0_source[4] = G_1[4]; 
  out_F_0_source[5] = G_1[5]; 
  out_F_0_source[6] = G_1[6]; 
  out_F_0_source[7] = G_1[7]; 
  out_F_0_source[8] = G_1[8]; 
  out_F_0_source[9] = G_1[9]; 
  out_F_0_source[10] = G_1[10]; 
  out_F_0_source[11] = G_1[11]; 
  out_F_0_source[12] = G_1[12]; 
  out_F_0_source[13] = G_1[13]; 
  out_F_0_source[14] = G_1[14]; 
  out_F_0_source[15] = G_1[15]; 
  out_F_0_source[16] = G_1[16]; 
  out_F_0_source[17] = G_1[17]; 
  out_F_0_source[18] = G_1[18]; 
  out_F_0_source[19] = G_1[19]; 
  out_F_0_source[20] = G_1[20]; 
  out_F_0_source[21] = G_1[21]; 
  out_F_0_source[22] = G_1[22]; 
  out_F_0_source[23] = G_1[23]; 
  out_G_1_source[0] = (-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[11])-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[7]-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[6]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[5]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[3]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[2]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[1]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[0]; 
  out_G_1_source[1] = (-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[11])-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[7]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[6]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[5]-0.7071067811865475*tmp_T_perp_F_1[3]*T_perp_over_m[5]-0.7071067811865475*tmp_T_perp_F_1[2]*T_perp_over_m[4]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[1]-0.7071067811865475*tmp_T_perp_F_1[0]*T_perp_over_m[1]; 
  out_G_1_source[2] = (-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[11])-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[7]-0.7071067811865475*tmp_T_perp_F_1[6]*T_perp_over_m[7]-0.7071067811865475*tmp_T_perp_F_1[3]*T_perp_over_m[6]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[5]-0.7071067811865475*tmp_T_perp_F_1[1]*T_perp_over_m[4]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[2]-0.7071067811865475*tmp_T_perp_F_1[0]*T_perp_over_m[2]; 
  out_G_1_source[3] = (-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[11])-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[7]-0.7071067811865475*tmp_T_perp_F_1[5]*T_perp_over_m[7]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[6]-0.7071067811865475*tmp_T_perp_F_1[2]*T_perp_over_m[6]-0.7071067811865475*tmp_T_perp_F_1[1]*T_perp_over_m[5]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[3]-0.7071067811865475*tmp_T_perp_F_1[0]*T_perp_over_m[3]; 
  out_G_1_source[4] = (-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[15])-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[14]-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[13]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[12]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[10]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[9]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[8]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[4]; 
  out_G_1_source[5] = (-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[11])-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[7]-0.7071067811865475*tmp_T_perp_F_1[3]*T_perp_over_m[7]-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[6]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[5]-0.7071067811865475*tmp_T_perp_F_1[0]*T_perp_over_m[4]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[2]-0.7071067811865475*tmp_T_perp_F_1[1]*T_perp_over_m[2]; 
  out_G_1_source[6] = (-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[11])-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[7]-0.7071067811865475*tmp_T_perp_F_1[2]*T_perp_over_m[7]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[6]-0.7071067811865475*tmp_T_perp_F_1[5]*T_perp_over_m[6]-0.7071067811865475*tmp_T_perp_F_1[0]*T_perp_over_m[5]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[3]-0.7071067811865475*tmp_T_perp_F_1[1]*T_perp_over_m[3]; 
  out_G_1_source[7] = (-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[11])-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[7]-0.7071067811865475*tmp_T_perp_F_1[1]*T_perp_over_m[7]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[6]-0.7071067811865475*tmp_T_perp_F_1[0]*T_perp_over_m[6]-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[5]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[3]-0.7071067811865475*tmp_T_perp_F_1[2]*T_perp_over_m[3]; 
  out_G_1_source[8] = (-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[15])-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[14]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[13]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[12]-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[10]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[9]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[8]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[4]; 
  out_G_1_source[9] = (-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[15])-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[14]-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[13]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[12]-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[10]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[9]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[8]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[4]; 
  out_G_1_source[10] = (-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[15])-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[14]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[13]-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[12]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[10]-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[9]-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[8]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[4]; 
  out_G_1_source[11] = (-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[11])-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[7]-0.7071067811865475*tmp_T_perp_F_1[0]*T_perp_over_m[7]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[6]-0.7071067811865475*tmp_T_perp_F_1[1]*T_perp_over_m[6]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[5]-0.7071067811865475*tmp_T_perp_F_1[2]*T_perp_over_m[5]-0.7071067811865475*tmp_T_perp_F_1[3]*T_perp_over_m[4]; 
  out_G_1_source[12] = (-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[15])-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[14]-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[13]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[12]-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[10]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[9]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[8]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[4]; 
  out_G_1_source[13] = (-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[15])-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[14]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[13]-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[12]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[10]-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[9]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[8]-0.7071067811865475*tmp_T_perp_F_1[4]*T_perp_over_m[5]; 
  out_G_1_source[14] = (-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[15])-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[14]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[13]-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[12]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[10]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[9]-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[8]-0.7071067811865475*tmp_T_perp_F_1[4]*T_perp_over_m[6]; 
  out_G_1_source[15] = (-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[15])-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[14]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[13]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[12]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[10]-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[9]-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[8]-0.7071067811865475*tmp_T_perp_F_1[4]*T_perp_over_m[7]; 
  out_G_1_source[16] = (-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[23])-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[22]-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[21]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[20]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[19]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[18]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[17]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[16]; 
  out_G_1_source[17] = (-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[23])-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[22]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[21]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[20]-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[19]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[18]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[17]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[16]; 
  out_G_1_source[18] = (-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[23])-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[22]-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[21]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[20]-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[19]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[18]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[17]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[16]; 
  out_G_1_source[19] = (-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[23])-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[22]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[21]-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[20]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[19]-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[18]-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[17]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[16]; 
  out_G_1_source[20] = (-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[23])-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[22]-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[21]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[20]-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[19]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[18]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[17]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[16]; 
  out_G_1_source[21] = (-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[23])-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[22]-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[21]-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[20]-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[19]-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[18]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[17]-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[16]; 
  out_G_1_source[22] = (-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[23])-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[22]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[21]-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[20]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[19]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[18]-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[17]-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[16]; 
  out_G_1_source[23] = (-0.7071067811865475*T_perp_over_m[0]*tmp_T_perp_F_1[23])-0.7071067811865475*T_perp_over_m[1]*tmp_T_perp_F_1[22]-0.7071067811865475*T_perp_over_m[2]*tmp_T_perp_F_1[21]-0.7071067811865475*T_perp_over_m[3]*tmp_T_perp_F_1[20]-0.7071067811865475*T_perp_over_m[4]*tmp_T_perp_F_1[19]-0.7071067811865475*T_perp_over_m[5]*tmp_T_perp_F_1[18]-0.7071067811865475*T_perp_over_m[6]*tmp_T_perp_F_1[17]-0.7071067811865475*T_perp_over_m[7]*tmp_T_perp_F_1[16]; 
} 
