#include <gkyl_mom_vlasov_sr_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[2]*volFact; 
  out[3] += 1.414213562373095*f[4]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  const double dv10 = 2.0/dxv[2]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[2]*dv10; 

  out[0] += (p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (p0_over_gamma[1]*f[5]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (p0_over_gamma[1]*f[6]+p0_over_gamma[0]*f[2])*volFact; 
  out[3] += (p0_over_gamma[1]*f[7]+p0_over_gamma[0]*f[4])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M2_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  out[0] += (gamma[2]*f[8]+gamma[1]*f[3]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0*gamma[2]*f[9]+gamma[1]*f[5]+gamma[0]*f[1])*volFact; 
  out[2] += (1.0*gamma[2]*f[10]+gamma[1]*f[6]+gamma[0]*f[2])*volFact; 
  out[3] += (gamma[2]*f[11]+gamma[1]*f[7]+gamma[0]*f[4])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M3i_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1+0.408248290463863*f[3]*dv1); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1+0.408248290463863*f[5]*dv1); 
  out[2] += volFact*(1.414213562373095*f[2]*wx1+0.408248290463863*f[6]*dv1); 
  out[3] += volFact*(1.414213562373095*f[4]*wx1+0.408248290463863*f[7]*dv1); 
} 
GKYL_CU_DH void vlasov_sr_Ni_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
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
  out[4] += (p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*volFact; 
  out[5] += (p0_over_gamma[1]*f[5]+p0_over_gamma[0]*f[1])*volFact; 
  out[6] += (p0_over_gamma[1]*f[6]+p0_over_gamma[0]*f[2])*volFact; 
  out[7] += (p0_over_gamma[1]*f[7]+p0_over_gamma[0]*f[4])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_Tij_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double dv10 = 2.0/dxv[2]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[2]*dv10; 

  out[0] += (gamma[2]*f[8]+gamma[1]*f[3]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0*gamma[2]*f[9]+gamma[1]*f[5]+gamma[0]*f[1])*volFact; 
  out[2] += (1.0*gamma[2]*f[10]+gamma[1]*f[6]+gamma[0]*f[2])*volFact; 
  out[3] += (gamma[2]*f[11]+gamma[1]*f[7]+gamma[0]*f[4])*volFact; 
  out[4] += volFact*(1.414213562373095*f[0]*wx1+0.408248290463863*f[3]*dv1); 
  out[5] += volFact*(1.414213562373095*f[1]*wx1+0.408248290463863*f[5]*dv1); 
  out[6] += volFact*(1.414213562373095*f[2]*wx1+0.408248290463863*f[6]*dv1); 
  out[7] += volFact*(1.414213562373095*f[4]*wx1+0.408248290463863*f[7]*dv1); 
  out[8] += volFact*(p0_over_gamma[1]*f[3]*wx1+f[0]*p0_over_gamma[0]*wx1+0.2581988897471612*p0_over_gamma[1]*f[8]*dv1+0.2886751345948129*p0_over_gamma[0]*f[3]*dv1+0.2886751345948129*f[0]*p0_over_gamma[1]*dv1); 
  out[9] += volFact*(p0_over_gamma[1]*f[5]*wx1+p0_over_gamma[0]*f[1]*wx1+0.2581988897471611*p0_over_gamma[1]*f[9]*dv1+0.2886751345948129*p0_over_gamma[0]*f[5]*dv1+0.2886751345948129*f[1]*p0_over_gamma[1]*dv1); 
  out[10] += volFact*(p0_over_gamma[1]*f[6]*wx1+p0_over_gamma[0]*f[2]*wx1+0.2581988897471611*p0_over_gamma[1]*f[10]*dv1+0.2886751345948129*p0_over_gamma[0]*f[6]*dv1+0.2886751345948129*p0_over_gamma[1]*f[2]*dv1); 
  out[11] += volFact*(p0_over_gamma[1]*f[7]*wx1+p0_over_gamma[0]*f[4]*wx1+0.2581988897471612*p0_over_gamma[1]*f[11]*dv1+0.2886751345948129*p0_over_gamma[0]*f[7]*dv1+0.2886751345948129*p0_over_gamma[1]*f[4]*dv1); 
} 
GKYL_CU_DH void vlasov_sr_int_mom_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += (2.0*gamma[2]*f[8]+2.0*gamma[1]*f[3]+2.0*f[0]*gamma[0])*volFact; 
  out[2] += volFact*(2.828427124746191*f[0]*wx1+0.8164965809277261*f[3]*dv1); 
} 
