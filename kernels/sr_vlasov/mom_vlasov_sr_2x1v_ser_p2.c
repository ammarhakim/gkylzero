#include <gkyl_mom_vlasov_sr_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[2]*volFact; 
  out[3] += 1.414213562373095*f[4]*volFact; 
  out[4] += 1.414213562373095*f[7]*volFact; 
  out[5] += 1.414213562373095*f[8]*volFact; 
  out[6] += 1.414213562373095*f[11]*volFact; 
  out[7] += 1.414213562373095*f[12]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  const double dv10 = 2.0/dxv[2]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[2]*dv10; 

  out[0] += (p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (p0_over_gamma[1]*f[5]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (p0_over_gamma[1]*f[6]+p0_over_gamma[0]*f[2])*volFact; 
  out[3] += (p0_over_gamma[1]*f[10]+p0_over_gamma[0]*f[4])*volFact; 
  out[4] += (1.0*p0_over_gamma[1]*f[13]+p0_over_gamma[0]*f[7])*volFact; 
  out[5] += (1.0*p0_over_gamma[1]*f[14]+p0_over_gamma[0]*f[8])*volFact; 
  out[6] += (1.0*p0_over_gamma[1]*f[17]+p0_over_gamma[0]*f[11])*volFact; 
  out[7] += (1.0*p0_over_gamma[1]*f[18]+p0_over_gamma[0]*f[12])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M2_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  out[0] += (gamma[2]*f[9]+gamma[1]*f[3]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0*gamma[2]*f[15]+gamma[1]*f[5]+gamma[0]*f[1])*volFact; 
  out[2] += (1.0*gamma[2]*f[16]+gamma[1]*f[6]+gamma[0]*f[2])*volFact; 
  out[3] += (gamma[2]*f[19]+gamma[1]*f[10]+gamma[0]*f[4])*volFact; 
  out[4] += (1.0*gamma[1]*f[13]+gamma[0]*f[7])*volFact; 
  out[5] += (1.0*gamma[1]*f[14]+gamma[0]*f[8])*volFact; 
  out[6] += (1.0*gamma[1]*f[17]+gamma[0]*f[11])*volFact; 
  out[7] += (1.0*gamma[1]*f[18]+gamma[0]*f[12])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M3i_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1+0.408248290463863*f[3]*dv1); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1+0.408248290463863*f[5]*dv1); 
  out[2] += volFact*(1.414213562373095*f[2]*wx1+0.408248290463863*f[6]*dv1); 
  out[3] += volFact*(1.414213562373095*f[4]*wx1+0.408248290463863*f[10]*dv1); 
  out[4] += volFact*(1.414213562373095*f[7]*wx1+0.408248290463863*f[13]*dv1); 
  out[5] += volFact*(1.414213562373095*f[8]*wx1+0.408248290463863*f[14]*dv1); 
  out[6] += volFact*(1.414213562373095*f[11]*wx1+0.408248290463863*f[17]*dv1); 
  out[7] += volFact*(1.414213562373095*f[12]*wx1+0.408248290463863*f[18]*dv1); 
} 
GKYL_CU_DH void vlasov_sr_Ni_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  const double dv10 = 2.0/dxv[2]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[2]*dv10; 

  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[2]*volFact; 
  out[3] += 1.414213562373095*f[4]*volFact; 
  out[4] += 1.414213562373095*f[7]*volFact; 
  out[5] += 1.414213562373095*f[8]*volFact; 
  out[6] += 1.414213562373095*f[11]*volFact; 
  out[7] += 1.414213562373095*f[12]*volFact; 
  out[8] += (p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*volFact; 
  out[9] += (p0_over_gamma[1]*f[5]+p0_over_gamma[0]*f[1])*volFact; 
  out[10] += (p0_over_gamma[1]*f[6]+p0_over_gamma[0]*f[2])*volFact; 
  out[11] += (p0_over_gamma[1]*f[10]+p0_over_gamma[0]*f[4])*volFact; 
  out[12] += (1.0*p0_over_gamma[1]*f[13]+p0_over_gamma[0]*f[7])*volFact; 
  out[13] += (1.0*p0_over_gamma[1]*f[14]+p0_over_gamma[0]*f[8])*volFact; 
  out[14] += (1.0*p0_over_gamma[1]*f[17]+p0_over_gamma[0]*f[11])*volFact; 
  out[15] += (1.0*p0_over_gamma[1]*f[18]+p0_over_gamma[0]*f[12])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_Tij_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double dv10 = 2.0/dxv[2]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[2]*dv10; 

  out[0] += (gamma[2]*f[9]+gamma[1]*f[3]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0*gamma[2]*f[15]+gamma[1]*f[5]+gamma[0]*f[1])*volFact; 
  out[2] += (1.0*gamma[2]*f[16]+gamma[1]*f[6]+gamma[0]*f[2])*volFact; 
  out[3] += (gamma[2]*f[19]+gamma[1]*f[10]+gamma[0]*f[4])*volFact; 
  out[4] += (1.0*gamma[1]*f[13]+gamma[0]*f[7])*volFact; 
  out[5] += (1.0*gamma[1]*f[14]+gamma[0]*f[8])*volFact; 
  out[6] += (1.0*gamma[1]*f[17]+gamma[0]*f[11])*volFact; 
  out[7] += (1.0*gamma[1]*f[18]+gamma[0]*f[12])*volFact; 
  out[8] += volFact*(1.414213562373095*f[0]*wx1+0.408248290463863*f[3]*dv1); 
  out[9] += volFact*(1.414213562373095*f[1]*wx1+0.408248290463863*f[5]*dv1); 
  out[10] += volFact*(1.414213562373095*f[2]*wx1+0.408248290463863*f[6]*dv1); 
  out[11] += volFact*(1.414213562373095*f[4]*wx1+0.408248290463863*f[10]*dv1); 
  out[12] += volFact*(1.414213562373095*f[7]*wx1+0.408248290463863*f[13]*dv1); 
  out[13] += volFact*(1.414213562373095*f[8]*wx1+0.408248290463863*f[14]*dv1); 
  out[14] += volFact*(1.414213562373095*f[11]*wx1+0.408248290463863*f[17]*dv1); 
  out[15] += volFact*(1.414213562373095*f[12]*wx1+0.408248290463863*f[18]*dv1); 
  out[16] += volFact*(p0_over_gamma[1]*f[3]*wx1+f[0]*p0_over_gamma[0]*wx1+0.2581988897471612*p0_over_gamma[1]*f[9]*dv1+0.2886751345948129*p0_over_gamma[0]*f[3]*dv1+0.2886751345948129*f[0]*p0_over_gamma[1]*dv1); 
  out[17] += volFact*(p0_over_gamma[1]*f[5]*wx1+p0_over_gamma[0]*f[1]*wx1+0.2581988897471611*p0_over_gamma[1]*f[15]*dv1+0.2886751345948129*p0_over_gamma[0]*f[5]*dv1+0.2886751345948129*f[1]*p0_over_gamma[1]*dv1); 
  out[18] += volFact*(p0_over_gamma[1]*f[6]*wx1+p0_over_gamma[0]*f[2]*wx1+0.2581988897471611*p0_over_gamma[1]*f[16]*dv1+0.2886751345948129*p0_over_gamma[0]*f[6]*dv1+0.2886751345948129*p0_over_gamma[1]*f[2]*dv1); 
  out[19] += volFact*(p0_over_gamma[1]*f[10]*wx1+p0_over_gamma[0]*f[4]*wx1+0.2581988897471612*p0_over_gamma[1]*f[19]*dv1+0.2886751345948129*p0_over_gamma[0]*f[10]*dv1+0.2886751345948129*p0_over_gamma[1]*f[4]*dv1); 
  out[20] += volFact*(1.0*p0_over_gamma[1]*f[13]*wx1+p0_over_gamma[0]*f[7]*wx1+0.2886751345948129*p0_over_gamma[0]*f[13]*dv1+0.2886751345948129*p0_over_gamma[1]*f[7]*dv1); 
  out[21] += volFact*(1.0*p0_over_gamma[1]*f[14]*wx1+p0_over_gamma[0]*f[8]*wx1+0.2886751345948129*p0_over_gamma[0]*f[14]*dv1+0.2886751345948129*p0_over_gamma[1]*f[8]*dv1); 
  out[22] += volFact*(1.0*p0_over_gamma[1]*f[17]*wx1+p0_over_gamma[0]*f[11]*wx1+0.2886751345948129*p0_over_gamma[0]*f[17]*dv1+0.2886751345948129*p0_over_gamma[1]*f[11]*dv1); 
  out[23] += volFact*(1.0*p0_over_gamma[1]*f[18]*wx1+p0_over_gamma[0]*f[12]*wx1+0.2886751345948129*p0_over_gamma[0]*f[18]*dv1+0.2886751345948129*p0_over_gamma[1]*f[12]*dv1); 
} 
GKYL_CU_DH void vlasov_sr_int_mom_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += (2.0*gamma[2]*f[9]+2.0*gamma[1]*f[3]+2.0*f[0]*gamma[0])*volFact; 
  out[2] += volFact*(2.828427124746191*f[0]*wx1+0.8164965809277261*f[3]*dv1); 
} 
