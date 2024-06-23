#include <gkyl_mom_vlasov_sr_kernels.h> 
#include <gkyl_basis_ser_1x_p1_exp_sq.h> 
GKYL_CU_DH void vlasov_sr_M0_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double dv10 = 2.0/dxv[1]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[2]*dv10; 

  out[0] += (p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (p0_over_gamma[1]*f[3]+p0_over_gamma[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_Ni_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double dv10 = 2.0/dxv[1]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[2]*dv10; 

  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += (p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[3] += (p0_over_gamma[1]*f[3]+p0_over_gamma[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_Energy_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += (gamma[2]*f[4]+gamma[1]*f[2]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0*gamma[2]*f[5]+gamma[1]*f[3]+gamma[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_Pressure_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma_inv, const double *gamma, const double *GammaV2, const double *V_drift, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double *V_drift_0 = &V_drift[0]; 
  double V_drift_0_sq[2] = {0.0}; 
  ser_1x_p1_exp_sq(V_drift_0, V_drift_0_sq); 
 
  double temp[6] = {0.0}; 
  double temp_sq[6] = {0.0}; 
  double p_fac[6] = {0.0}; 
  temp[0] = 1.414213562373095*V_drift_0[0]*wx1; 
  temp[1] = 1.414213562373095*V_drift_0[1]*wx1; 
  temp[2] = 0.408248290463863*V_drift_0[0]*dv1; 
  temp[3] = 0.408248290463863*V_drift_0[1]*dv1; 

  temp_sq[0] = GammaV2[1]*V_drift_0_sq[1]*wx1_sq+GammaV2[0]*V_drift_0_sq[0]*wx1_sq+0.08333333333333333*GammaV2[1]*V_drift_0_sq[1]*dv1_sq+0.08333333333333333*GammaV2[0]*V_drift_0_sq[0]*dv1_sq; 
  temp_sq[1] = GammaV2[0]*V_drift_0_sq[1]*wx1_sq+V_drift_0_sq[0]*GammaV2[1]*wx1_sq+0.08333333333333333*GammaV2[0]*V_drift_0_sq[1]*dv1_sq+0.08333333333333333*V_drift_0_sq[0]*GammaV2[1]*dv1_sq; 
  temp_sq[2] = 0.5773502691896258*GammaV2[1]*V_drift_0_sq[1]*dv1*wx1+0.5773502691896258*GammaV2[0]*V_drift_0_sq[0]*dv1*wx1; 
  temp_sq[3] = 0.5773502691896258*GammaV2[0]*V_drift_0_sq[1]*dv1*wx1+0.5773502691896258*V_drift_0_sq[0]*GammaV2[1]*dv1*wx1; 
  temp_sq[4] = 0.07453559924999298*GammaV2[1]*V_drift_0_sq[1]*dv1_sq+0.07453559924999298*GammaV2[0]*V_drift_0_sq[0]*dv1_sq; 
  temp_sq[5] = 0.07453559924999302*GammaV2[0]*V_drift_0_sq[1]*dv1_sq+0.07453559924999302*V_drift_0_sq[0]*GammaV2[1]*dv1_sq; 

  p_fac[0] = 0.7071067811865475*gamma_inv[2]*temp_sq[4]+0.7071067811865475*gamma_inv[1]*temp_sq[2]-1.414213562373095*GammaV2[1]*temp[1]+GammaV2[0]*gamma[0]+0.7071067811865475*gamma_inv[0]*temp_sq[0]-1.414213562373095*GammaV2[0]*temp[0]-1.414213562373095*gamma_inv[0]; 
  p_fac[1] = 0.7071067811865475*gamma_inv[2]*temp_sq[5]+0.7071067811865475*gamma_inv[1]*temp_sq[3]+0.7071067811865475*gamma_inv[0]*temp_sq[1]-1.414213562373095*GammaV2[0]*temp[1]+gamma[0]*GammaV2[1]-1.414213562373095*temp[0]*GammaV2[1]; 
  p_fac[2] = 0.6324555320336759*gamma_inv[1]*temp_sq[4]-1.414213562373095*GammaV2[1]*temp[3]+0.6324555320336759*gamma_inv[2]*temp_sq[2]+0.7071067811865475*gamma_inv[0]*temp_sq[2]-1.414213562373095*GammaV2[0]*temp[2]+GammaV2[0]*gamma[1]+0.7071067811865475*temp_sq[0]*gamma_inv[1]-1.414213562373095*gamma_inv[1]; 
  p_fac[3] = 0.632455532033676*gamma_inv[1]*temp_sq[5]+0.6324555320336759*gamma_inv[2]*temp_sq[3]+0.7071067811865475*gamma_inv[0]*temp_sq[3]-1.414213562373095*GammaV2[0]*temp[3]-1.414213562373095*GammaV2[1]*temp[2]+GammaV2[1]*gamma[1]+0.7071067811865475*gamma_inv[1]*temp_sq[1]; 
  p_fac[4] = 0.4517539514526256*gamma_inv[2]*temp_sq[4]+0.7071067811865475*gamma_inv[0]*temp_sq[4]+GammaV2[0]*gamma[2]+0.6324555320336759*gamma_inv[1]*temp_sq[2]+0.7071067811865475*temp_sq[0]*gamma_inv[2]-1.414213562373095*gamma_inv[2]; 
  p_fac[5] = 0.4517539514526256*gamma_inv[2]*temp_sq[5]+0.7071067811865475*gamma_inv[0]*temp_sq[5]+0.632455532033676*gamma_inv[1]*temp_sq[3]+1.0*GammaV2[1]*gamma[2]+0.7071067811865475*temp_sq[1]*gamma_inv[2]; 

  out[0] += (0.7071067811865475*f[5]*p_fac[5]+0.7071067811865475*f[4]*p_fac[4]+0.7071067811865475*f[3]*p_fac[3]+0.7071067811865475*f[2]*p_fac[2]+0.7071067811865475*f[1]*p_fac[1]+0.7071067811865475*f[0]*p_fac[0])*volFact; 
  out[1] += (0.7071067811865475*f[4]*p_fac[5]+0.7071067811865475*p_fac[4]*f[5]+0.7071067811865475*f[2]*p_fac[3]+0.7071067811865475*p_fac[2]*f[3]+0.7071067811865475*f[0]*p_fac[1]+0.7071067811865475*p_fac[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_Tij_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double dv10 = 2.0/dxv[1]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[2]*dv10; 

  out[0] += (gamma[2]*f[4]+gamma[1]*f[2]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0*gamma[2]*f[5]+gamma[1]*f[3]+gamma[0]*f[1])*volFact; 
  out[2] += volFact*(1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1); 
  out[3] += volFact*(1.414213562373095*f[1]*wx1+0.408248290463863*f[3]*dv1); 
  out[4] += volFact*(p0_over_gamma[1]*f[2]*wx1+f[0]*p0_over_gamma[0]*wx1+0.2581988897471612*p0_over_gamma[1]*f[4]*dv1+0.2886751345948129*p0_over_gamma[0]*f[2]*dv1+0.2886751345948129*f[0]*p0_over_gamma[1]*dv1); 
  out[5] += volFact*(p0_over_gamma[1]*f[3]*wx1+p0_over_gamma[0]*f[1]*wx1+0.2581988897471611*p0_over_gamma[1]*f[5]*dv1+0.2886751345948129*p0_over_gamma[0]*f[3]*dv1+0.2886751345948129*f[1]*p0_over_gamma[1]*dv1); 
} 
GKYL_CU_DH void vlasov_sr_int_mom_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*0.25; 
  const double wx1 = w[1], dv1 = dxv[1]; 
 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += (1.414213562373095*gamma[2]*f[4]+1.414213562373095*gamma[1]*f[2]+1.414213562373095*f[0]*gamma[0])*volFact; 
  out[2] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[2]*dv1); 
} 
