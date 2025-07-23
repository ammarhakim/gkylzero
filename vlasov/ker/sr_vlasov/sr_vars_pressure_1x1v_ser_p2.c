#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_pressure_1x1v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure) 
{ 
  // gamma:       Particle Lorentz boost factor sqrt(1 + p^2).
  // gamma_inv:   Inverse particle Lorentz boost factor 1/sqrt(1 + p^2).
  // u_i:         Spatial components of bulk four-velocity = GammaV*V_drift. 
  // u_i_sq:      Squared spatial components of bulk four-velocity = u_i^2. 
  // GammaV:      Bulk four-velocity Lorentz factor = sqrt(1 + |u_i|^2). 
  // GammaV_sq:   Squared bulk four-velocity Lorentz factor = 1 + |u_i|^2. 
  // f:           Input distribution function.
  // sr_pressure: Output relativistic pressure.
  const double volFact = dxv[1]/2; 
 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double *V_0 = &u_i[0]; 
  const double *V_0_sq = &u_i_sq[0]; 
 
  double temp[8] = {0.0}; 
  double temp_sq[8] = {0.0}; 
  double p_fac[8] = {0.0}; 
  temp[0] = 1.414213562373095*V_0[0]*wx1; 
  temp[1] = 1.414213562373095*V_0[1]*wx1; 
  temp[2] = 0.408248290463863*V_0[0]*dv1; 
  temp[3] = 0.408248290463863*V_0[1]*dv1; 
  temp[4] = 1.414213562373095*V_0[2]*wx1; 
  temp[6] = 0.408248290463863*V_0[2]*dv1; 

  temp_sq[0] = 1.414213562373095*V_0_sq[0]*wx1_sq+0.1178511301977579*V_0_sq[0]*dv1_sq; 
  temp_sq[1] = 1.414213562373095*V_0_sq[1]*wx1_sq+0.1178511301977579*V_0_sq[1]*dv1_sq; 
  temp_sq[2] = 0.8164965809277261*V_0_sq[0]*dv1*wx1; 
  temp_sq[3] = 0.8164965809277261*V_0_sq[1]*dv1*wx1; 
  temp_sq[4] = 1.414213562373095*V_0_sq[2]*wx1_sq+0.1178511301977579*V_0_sq[2]*dv1_sq; 
  temp_sq[5] = 0.105409255338946*V_0_sq[0]*dv1_sq; 
  temp_sq[6] = 0.816496580927726*V_0_sq[2]*dv1*wx1; 
  temp_sq[7] = 0.105409255338946*V_0_sq[1]*dv1_sq; 

  p_fac[0] = 0.7071067811865475*gamma_inv[2]*temp_sq[5]-1.414213562373095*GammaV[2]*temp[4]+0.7071067811865475*gamma_inv[1]*temp_sq[2]-1.414213562373095*GammaV[1]*temp[1]+GammaV_sq[0]*gamma[0]+0.7071067811865475*gamma_inv[0]*temp_sq[0]-1.414213562373095*GammaV[0]*temp[0]-1.414213562373095*gamma_inv[0]; 
  p_fac[1] = 0.7071067811865475*gamma_inv[2]*temp_sq[7]-1.264911064067352*GammaV[1]*temp[4]+0.7071067811865475*gamma_inv[1]*temp_sq[3]-1.264911064067352*temp[1]*GammaV[2]+0.7071067811865475*gamma_inv[0]*temp_sq[1]-1.414213562373095*GammaV[0]*temp[1]+gamma[0]*GammaV_sq[1]-1.414213562373095*temp[0]*GammaV[1]; 
  p_fac[2] = (-1.414213562373095*GammaV[2]*temp[6])+0.6324555320336759*gamma_inv[1]*temp_sq[5]-1.414213562373095*GammaV[1]*temp[3]+0.6324555320336759*gamma_inv[2]*temp_sq[2]+0.7071067811865475*gamma_inv[0]*temp_sq[2]-1.414213562373095*GammaV[0]*temp[2]+GammaV_sq[0]*gamma[1]+0.7071067811865475*temp_sq[0]*gamma_inv[1]-1.414213562373095*gamma_inv[1]; 
  p_fac[3] = 0.632455532033676*gamma_inv[1]*temp_sq[7]-1.264911064067352*GammaV[1]*temp[6]+0.6324555320336759*gamma_inv[2]*temp_sq[3]+0.7071067811865475*gamma_inv[0]*temp_sq[3]-1.264911064067352*GammaV[2]*temp[3]-1.414213562373095*GammaV[0]*temp[3]-1.414213562373095*GammaV[1]*temp[2]+GammaV_sq[1]*gamma[1]+0.7071067811865475*gamma_inv[1]*temp_sq[1]; 
  p_fac[4] = 0.7071067811865475*gamma_inv[1]*temp_sq[6]+0.7071067811865475*gamma_inv[0]*temp_sq[4]-0.9035079029052515*GammaV[2]*temp[4]-1.414213562373095*GammaV[0]*temp[4]+gamma[0]*GammaV_sq[2]-1.414213562373095*temp[0]*GammaV[2]-1.264911064067352*GammaV[1]*temp[1]; 
  p_fac[5] = 0.4517539514526256*gamma_inv[2]*temp_sq[5]+0.7071067811865475*gamma_inv[0]*temp_sq[5]+GammaV_sq[0]*gamma[2]+0.6324555320336759*gamma_inv[1]*temp_sq[2]+0.7071067811865475*temp_sq[0]*gamma_inv[2]-1.414213562373095*gamma_inv[2]; 
  p_fac[6] = 0.6324555320336759*gamma_inv[2]*temp_sq[6]+0.7071067811865475*gamma_inv[0]*temp_sq[6]-0.9035079029052515*GammaV[2]*temp[6]-1.414213562373095*GammaV[0]*temp[6]+0.7071067811865475*gamma_inv[1]*temp_sq[4]-1.264911064067352*GammaV[1]*temp[3]-1.414213562373095*GammaV[2]*temp[2]+1.0*gamma[1]*GammaV_sq[2]; 
  p_fac[7] = 0.4517539514526256*gamma_inv[2]*temp_sq[7]+0.7071067811865475*gamma_inv[0]*temp_sq[7]+0.632455532033676*gamma_inv[1]*temp_sq[3]+1.0*GammaV_sq[1]*gamma[2]+0.7071067811865475*temp_sq[1]*gamma_inv[2]; 

  sr_pressure[0] += (0.7071067811865475*f[7]*p_fac[7]+0.7071067811865475*f[6]*p_fac[6]+0.7071067811865475*f[5]*p_fac[5]+0.7071067811865475*f[4]*p_fac[4]+0.7071067811865475*f[3]*p_fac[3]+0.7071067811865475*f[2]*p_fac[2]+0.7071067811865475*f[1]*p_fac[1]+0.7071067811865475*f[0]*p_fac[0])*volFact; 
  sr_pressure[1] += (0.7071067811865475*f[5]*p_fac[7]+0.7071067811865475*p_fac[5]*f[7]+0.632455532033676*f[3]*p_fac[6]+0.632455532033676*p_fac[3]*f[6]+0.6324555320336759*f[1]*p_fac[4]+0.6324555320336759*p_fac[1]*f[4]+0.7071067811865475*f[2]*p_fac[3]+0.7071067811865475*p_fac[2]*f[3]+0.7071067811865475*f[0]*p_fac[1]+0.7071067811865475*p_fac[0]*f[1])*volFact; 
  sr_pressure[2] += (0.6324555320336759*f[7]*p_fac[7]+0.4517539514526256*f[6]*p_fac[6]+0.7071067811865475*f[2]*p_fac[6]+0.7071067811865475*p_fac[2]*f[6]+0.4517539514526256*f[4]*p_fac[4]+0.7071067811865475*f[0]*p_fac[4]+0.7071067811865475*p_fac[0]*f[4]+0.6324555320336759*f[3]*p_fac[3]+0.6324555320336759*f[1]*p_fac[1])*volFact; 
} 
