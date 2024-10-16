#include <gkyl_mom_vlasov_sr_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  double p0_over_gamma[8] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[4]*dv10; 
  p0_over_gamma[2] = 1.732050807568877*gamma[3]*dv10; 
  p0_over_gamma[3] = 3.872983346207417*gamma[6]*dv10; 
  p0_over_gamma[5] = 1.732050807568877*gamma[7]*dv10; 

  double p1_over_gamma[8] = {0.0}; 
  p1_over_gamma[0] = 1.732050807568877*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.732050807568877*gamma[3]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[5]*dv11; 
  p1_over_gamma[3] = 3.872983346207417*gamma[7]*dv11; 
  p1_over_gamma[4] = 1.732050807568877*gamma[6]*dv11; 

  out[0] += (p0_over_gamma[5]*f[24]+p0_over_gamma[3]*f[10]+p0_over_gamma[2]*f[4]+p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (1.0*p0_over_gamma[5]*f[25]+p0_over_gamma[3]*f[13]+p0_over_gamma[2]*f[8]+p0_over_gamma[1]*f[6]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (1.0*p0_over_gamma[5]*f[26]+p0_over_gamma[3]*f[14]+p0_over_gamma[2]*f[9]+p0_over_gamma[1]*f[7]+p0_over_gamma[0]*f[2])*volFact; 
  out[3] += (p0_over_gamma[5]*f[28]+p0_over_gamma[3]*f[15]+p0_over_gamma[2]*f[12]+p0_over_gamma[1]*f[11]+p0_over_gamma[0]*f[5])*volFact; 
  out[4] += (p1_over_gamma[4]*f[16]+p1_over_gamma[3]*f[10]+p1_over_gamma[2]*f[4]+p1_over_gamma[1]*f[3]+f[0]*p1_over_gamma[0])*volFact; 
  out[5] += (1.0*p1_over_gamma[4]*f[17]+p1_over_gamma[3]*f[13]+p1_over_gamma[2]*f[8]+p1_over_gamma[1]*f[6]+p1_over_gamma[0]*f[1])*volFact; 
  out[6] += (1.0*p1_over_gamma[4]*f[18]+p1_over_gamma[3]*f[14]+p1_over_gamma[2]*f[9]+p1_over_gamma[1]*f[7]+p1_over_gamma[0]*f[2])*volFact; 
  out[7] += (p1_over_gamma[4]*f[20]+p1_over_gamma[3]*f[15]+p1_over_gamma[2]*f[12]+p1_over_gamma[1]*f[11]+p1_over_gamma[0]*f[5])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  out[0] += (gamma[7]*f[27]+gamma[5]*f[24]+gamma[6]*f[19]+gamma[4]*f[16]+gamma[3]*f[10]+gamma[2]*f[4]+gamma[1]*f[3]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0*gamma[7]*f[29]+1.0*gamma[5]*f[25]+1.0*gamma[6]*f[21]+1.0*gamma[4]*f[17]+gamma[3]*f[13]+gamma[2]*f[8]+gamma[1]*f[6]+gamma[0]*f[1])*volFact; 
  out[2] += (1.0*gamma[7]*f[30]+1.0*gamma[5]*f[26]+1.0*gamma[6]*f[22]+1.0*gamma[4]*f[18]+gamma[3]*f[14]+gamma[2]*f[9]+gamma[1]*f[7]+gamma[0]*f[2])*volFact; 
  out[3] += (gamma[7]*f[31]+gamma[5]*f[28]+gamma[6]*f[23]+gamma[4]*f[20]+gamma[3]*f[15]+gamma[2]*f[12]+gamma[1]*f[11]+gamma[0]*f[5])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M3i_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1); 
  out[2] += volFact*(2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1); 
  out[3] += volFact*(2.0*f[5]*wx1+0.5773502691896258*f[11]*dv1); 
  out[4] += volFact*(2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2); 
  out[5] += volFact*(2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2); 
  out[6] += volFact*(2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2); 
  out[7] += volFact*(2.0*f[5]*wx2+0.5773502691896258*f[12]*dv2); 
} 
GKYL_CU_DH void vlasov_sr_Ni_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  double p0_over_gamma[8] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[4]*dv10; 
  p0_over_gamma[2] = 1.732050807568877*gamma[3]*dv10; 
  p0_over_gamma[3] = 3.872983346207417*gamma[6]*dv10; 
  p0_over_gamma[5] = 1.732050807568877*gamma[7]*dv10; 

  double p1_over_gamma[8] = {0.0}; 
  p1_over_gamma[0] = 1.732050807568877*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.732050807568877*gamma[3]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[5]*dv11; 
  p1_over_gamma[3] = 3.872983346207417*gamma[7]*dv11; 
  p1_over_gamma[4] = 1.732050807568877*gamma[6]*dv11; 

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += (p0_over_gamma[5]*f[24]+p0_over_gamma[3]*f[10]+p0_over_gamma[2]*f[4]+p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*volFact; 
  out[5] += (1.0*p0_over_gamma[5]*f[25]+p0_over_gamma[3]*f[13]+p0_over_gamma[2]*f[8]+p0_over_gamma[1]*f[6]+p0_over_gamma[0]*f[1])*volFact; 
  out[6] += (1.0*p0_over_gamma[5]*f[26]+p0_over_gamma[3]*f[14]+p0_over_gamma[2]*f[9]+p0_over_gamma[1]*f[7]+p0_over_gamma[0]*f[2])*volFact; 
  out[7] += (p0_over_gamma[5]*f[28]+p0_over_gamma[3]*f[15]+p0_over_gamma[2]*f[12]+p0_over_gamma[1]*f[11]+p0_over_gamma[0]*f[5])*volFact; 
  out[8] += (p1_over_gamma[4]*f[16]+p1_over_gamma[3]*f[10]+p1_over_gamma[2]*f[4]+p1_over_gamma[1]*f[3]+f[0]*p1_over_gamma[0])*volFact; 
  out[9] += (1.0*p1_over_gamma[4]*f[17]+p1_over_gamma[3]*f[13]+p1_over_gamma[2]*f[8]+p1_over_gamma[1]*f[6]+p1_over_gamma[0]*f[1])*volFact; 
  out[10] += (1.0*p1_over_gamma[4]*f[18]+p1_over_gamma[3]*f[14]+p1_over_gamma[2]*f[9]+p1_over_gamma[1]*f[7]+p1_over_gamma[0]*f[2])*volFact; 
  out[11] += (p1_over_gamma[4]*f[20]+p1_over_gamma[3]*f[15]+p1_over_gamma[2]*f[12]+p1_over_gamma[1]*f[11]+p1_over_gamma[0]*f[5])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_Tij_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double dv10 = 2.0/dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double dv11 = 2.0/dxv[3]; 
  double p0_over_gamma[8] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[4]*dv10; 
  p0_over_gamma[2] = 1.732050807568877*gamma[3]*dv10; 
  p0_over_gamma[3] = 3.872983346207417*gamma[6]*dv10; 
  p0_over_gamma[5] = 1.732050807568877*gamma[7]*dv10; 

  double p1_over_gamma[8] = {0.0}; 
  p1_over_gamma[0] = 1.732050807568877*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.732050807568877*gamma[3]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[5]*dv11; 
  p1_over_gamma[3] = 3.872983346207417*gamma[7]*dv11; 
  p1_over_gamma[4] = 1.732050807568877*gamma[6]*dv11; 

  out[0] += (gamma[7]*f[27]+gamma[5]*f[24]+gamma[6]*f[19]+gamma[4]*f[16]+gamma[3]*f[10]+gamma[2]*f[4]+gamma[1]*f[3]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0*gamma[7]*f[29]+1.0*gamma[5]*f[25]+1.0*gamma[6]*f[21]+1.0*gamma[4]*f[17]+gamma[3]*f[13]+gamma[2]*f[8]+gamma[1]*f[6]+gamma[0]*f[1])*volFact; 
  out[2] += (1.0*gamma[7]*f[30]+1.0*gamma[5]*f[26]+1.0*gamma[6]*f[22]+1.0*gamma[4]*f[18]+gamma[3]*f[14]+gamma[2]*f[9]+gamma[1]*f[7]+gamma[0]*f[2])*volFact; 
  out[3] += (gamma[7]*f[31]+gamma[5]*f[28]+gamma[6]*f[23]+gamma[4]*f[20]+gamma[3]*f[15]+gamma[2]*f[12]+gamma[1]*f[11]+gamma[0]*f[5])*volFact; 
  out[4] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1); 
  out[5] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1); 
  out[6] += volFact*(2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1); 
  out[7] += volFact*(2.0*f[5]*wx1+0.5773502691896258*f[11]*dv1); 
  out[8] += volFact*(2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2); 
  out[9] += volFact*(2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2); 
  out[10] += volFact*(2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2); 
  out[11] += volFact*(2.0*f[5]*wx2+0.5773502691896258*f[12]*dv2); 
  out[12] += volFact*(p0_over_gamma[5]*f[24]*wx1+p0_over_gamma[3]*f[10]*wx1+p0_over_gamma[2]*f[4]*wx1+p0_over_gamma[1]*f[3]*wx1+f[0]*p0_over_gamma[0]*wx1+0.2886751345948129*p0_over_gamma[5]*f[27]*dv1+0.2581988897471611*p0_over_gamma[3]*f[19]*dv1+0.2581988897471612*p0_over_gamma[1]*f[16]*dv1+0.2886751345948129*p0_over_gamma[2]*f[10]*dv1+0.2886751345948129*p0_over_gamma[3]*f[4]*dv1+0.2886751345948129*p0_over_gamma[0]*f[3]*dv1+0.2886751345948129*f[0]*p0_over_gamma[1]*dv1); 
  out[13] += volFact*(1.0*p0_over_gamma[5]*f[25]*wx1+p0_over_gamma[3]*f[13]*wx1+p0_over_gamma[2]*f[8]*wx1+p0_over_gamma[1]*f[6]*wx1+p0_over_gamma[0]*f[1]*wx1+0.2886751345948129*p0_over_gamma[5]*f[29]*dv1+0.2581988897471612*p0_over_gamma[3]*f[21]*dv1+0.2581988897471611*p0_over_gamma[1]*f[17]*dv1+0.2886751345948129*p0_over_gamma[2]*f[13]*dv1+0.2886751345948129*p0_over_gamma[3]*f[8]*dv1+0.2886751345948129*p0_over_gamma[0]*f[6]*dv1+0.2886751345948129*f[1]*p0_over_gamma[1]*dv1); 
  out[14] += volFact*(1.0*p0_over_gamma[5]*f[26]*wx1+p0_over_gamma[3]*f[14]*wx1+p0_over_gamma[2]*f[9]*wx1+p0_over_gamma[1]*f[7]*wx1+p0_over_gamma[0]*f[2]*wx1+0.2886751345948129*p0_over_gamma[5]*f[30]*dv1+0.2581988897471612*p0_over_gamma[3]*f[22]*dv1+0.2581988897471611*p0_over_gamma[1]*f[18]*dv1+0.2886751345948129*p0_over_gamma[2]*f[14]*dv1+0.2886751345948129*p0_over_gamma[3]*f[9]*dv1+0.2886751345948129*p0_over_gamma[0]*f[7]*dv1+0.2886751345948129*p0_over_gamma[1]*f[2]*dv1); 
  out[15] += volFact*(p0_over_gamma[5]*f[28]*wx1+p0_over_gamma[3]*f[15]*wx1+p0_over_gamma[2]*f[12]*wx1+p0_over_gamma[1]*f[11]*wx1+p0_over_gamma[0]*f[5]*wx1+0.2886751345948129*p0_over_gamma[5]*f[31]*dv1+0.2581988897471611*p0_over_gamma[3]*f[23]*dv1+0.2581988897471612*p0_over_gamma[1]*f[20]*dv1+0.2886751345948129*p0_over_gamma[2]*f[15]*dv1+0.2886751345948129*p0_over_gamma[3]*f[12]*dv1+0.2886751345948129*p0_over_gamma[0]*f[11]*dv1+0.2886751345948129*p0_over_gamma[1]*f[5]*dv1); 
  out[16] += volFact*(p0_over_gamma[5]*f[24]*wx2+p0_over_gamma[3]*f[10]*wx2+p0_over_gamma[2]*f[4]*wx2+p0_over_gamma[1]*f[3]*wx2+f[0]*p0_over_gamma[0]*wx2+0.2581988897471611*p0_over_gamma[3]*f[27]*dv2+0.2581988897471612*p0_over_gamma[2]*f[24]*dv2+0.2886751345948129*p0_over_gamma[1]*f[10]*dv2+0.2581988897471612*f[4]*p0_over_gamma[5]*dv2+0.2886751345948129*p0_over_gamma[0]*f[4]*dv2+0.2886751345948129*f[3]*p0_over_gamma[3]*dv2+0.2886751345948129*f[0]*p0_over_gamma[2]*dv2); 
  out[17] += volFact*(1.0*p0_over_gamma[5]*f[25]*wx2+p0_over_gamma[3]*f[13]*wx2+p0_over_gamma[2]*f[8]*wx2+p0_over_gamma[1]*f[6]*wx2+p0_over_gamma[0]*f[1]*wx2+0.2581988897471612*p0_over_gamma[3]*f[29]*dv2+0.2581988897471611*p0_over_gamma[2]*f[25]*dv2+0.2886751345948129*p0_over_gamma[1]*f[13]*dv2+0.2581988897471612*p0_over_gamma[5]*f[8]*dv2+0.2886751345948129*p0_over_gamma[0]*f[8]*dv2+0.2886751345948129*p0_over_gamma[3]*f[6]*dv2+0.2886751345948129*f[1]*p0_over_gamma[2]*dv2); 
  out[18] += volFact*(1.0*p0_over_gamma[5]*f[26]*wx2+p0_over_gamma[3]*f[14]*wx2+p0_over_gamma[2]*f[9]*wx2+p0_over_gamma[1]*f[7]*wx2+p0_over_gamma[0]*f[2]*wx2+0.2581988897471612*p0_over_gamma[3]*f[30]*dv2+0.2581988897471611*p0_over_gamma[2]*f[26]*dv2+0.2886751345948129*p0_over_gamma[1]*f[14]*dv2+0.2581988897471612*p0_over_gamma[5]*f[9]*dv2+0.2886751345948129*p0_over_gamma[0]*f[9]*dv2+0.2886751345948129*p0_over_gamma[3]*f[7]*dv2+0.2886751345948129*f[2]*p0_over_gamma[2]*dv2); 
  out[19] += volFact*(p0_over_gamma[5]*f[28]*wx2+p0_over_gamma[3]*f[15]*wx2+p0_over_gamma[2]*f[12]*wx2+p0_over_gamma[1]*f[11]*wx2+p0_over_gamma[0]*f[5]*wx2+0.2581988897471611*p0_over_gamma[3]*f[31]*dv2+0.2581988897471612*p0_over_gamma[2]*f[28]*dv2+0.2886751345948129*p0_over_gamma[1]*f[15]*dv2+0.2581988897471612*p0_over_gamma[5]*f[12]*dv2+0.2886751345948129*p0_over_gamma[0]*f[12]*dv2+0.2886751345948129*p0_over_gamma[3]*f[11]*dv2+0.2886751345948129*p0_over_gamma[2]*f[5]*dv2); 
  out[20] += volFact*(p1_over_gamma[4]*f[16]*wx2+p1_over_gamma[3]*f[10]*wx2+p1_over_gamma[2]*f[4]*wx2+p1_over_gamma[1]*f[3]*wx2+f[0]*p1_over_gamma[0]*wx2+0.2581988897471611*p1_over_gamma[3]*f[27]*dv2+0.2581988897471612*p1_over_gamma[2]*f[24]*dv2+0.2886751345948129*p1_over_gamma[4]*f[19]*dv2+0.2886751345948129*p1_over_gamma[1]*f[10]*dv2+0.2886751345948129*p1_over_gamma[0]*f[4]*dv2+0.2886751345948129*f[3]*p1_over_gamma[3]*dv2+0.2886751345948129*f[0]*p1_over_gamma[2]*dv2); 
  out[21] += volFact*(1.0*p1_over_gamma[4]*f[17]*wx2+p1_over_gamma[3]*f[13]*wx2+p1_over_gamma[2]*f[8]*wx2+p1_over_gamma[1]*f[6]*wx2+p1_over_gamma[0]*f[1]*wx2+0.2581988897471612*p1_over_gamma[3]*f[29]*dv2+0.2581988897471611*p1_over_gamma[2]*f[25]*dv2+0.2886751345948129*p1_over_gamma[4]*f[21]*dv2+0.2886751345948129*p1_over_gamma[1]*f[13]*dv2+0.2886751345948129*p1_over_gamma[0]*f[8]*dv2+0.2886751345948129*p1_over_gamma[3]*f[6]*dv2+0.2886751345948129*f[1]*p1_over_gamma[2]*dv2); 
  out[22] += volFact*(1.0*p1_over_gamma[4]*f[18]*wx2+p1_over_gamma[3]*f[14]*wx2+p1_over_gamma[2]*f[9]*wx2+p1_over_gamma[1]*f[7]*wx2+p1_over_gamma[0]*f[2]*wx2+0.2581988897471612*p1_over_gamma[3]*f[30]*dv2+0.2581988897471611*p1_over_gamma[2]*f[26]*dv2+0.2886751345948129*p1_over_gamma[4]*f[22]*dv2+0.2886751345948129*p1_over_gamma[1]*f[14]*dv2+0.2886751345948129*p1_over_gamma[0]*f[9]*dv2+0.2886751345948129*p1_over_gamma[3]*f[7]*dv2+0.2886751345948129*f[2]*p1_over_gamma[2]*dv2); 
  out[23] += volFact*(p1_over_gamma[4]*f[20]*wx2+p1_over_gamma[3]*f[15]*wx2+p1_over_gamma[2]*f[12]*wx2+p1_over_gamma[1]*f[11]*wx2+p1_over_gamma[0]*f[5]*wx2+0.2581988897471611*p1_over_gamma[3]*f[31]*dv2+0.2581988897471612*p1_over_gamma[2]*f[28]*dv2+0.2886751345948129*p1_over_gamma[4]*f[23]*dv2+0.2886751345948129*p1_over_gamma[1]*f[15]*dv2+0.2886751345948129*p1_over_gamma[0]*f[12]*dv2+0.2886751345948129*p1_over_gamma[3]*f[11]*dv2+0.2886751345948129*p1_over_gamma[2]*f[5]*dv2); 
} 
GKYL_CU_DH void vlasov_sr_int_mom_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (2.0*gamma[7]*f[27]+2.0*gamma[5]*f[24]+2.0*gamma[6]*f[19]+2.0*gamma[4]*f[16]+2.0*gamma[3]*f[10]+2.0*gamma[2]*f[4]+2.0*gamma[1]*f[3]+2.0*f[0]*gamma[0])*volFact; 
  out[2] += volFact*(4.0*f[0]*wx1+1.154700538379252*f[3]*dv1); 
  out[3] += volFact*(4.0*f[0]*wx2+1.154700538379252*f[4]*dv2); 
} 
