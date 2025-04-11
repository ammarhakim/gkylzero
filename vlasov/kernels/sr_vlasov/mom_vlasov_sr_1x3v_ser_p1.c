#include <gkyl_mom_vlasov_sr_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  double p0_over_gamma[20] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[7]*dv10; 
  p0_over_gamma[2] = 1.732050807568877*gamma[4]*dv10; 
  p0_over_gamma[3] = 1.732050807568877*gamma[5]*dv10; 
  p0_over_gamma[4] = 3.872983346207417*gamma[11]*dv10; 
  p0_over_gamma[5] = 3.872983346207417*gamma[13]*dv10; 
  p0_over_gamma[6] = 1.732050807568877*gamma[10]*dv10; 
  p0_over_gamma[8] = 1.732050807568877*gamma[12]*dv10; 
  p0_over_gamma[9] = 1.732050807568877*gamma[15]*dv10; 
  p0_over_gamma[10] = 3.872983346207417*gamma[17]*dv10; 
  p0_over_gamma[14] = 1.732050807568877*gamma[18]*dv10; 
  p0_over_gamma[16] = 1.732050807568877*gamma[19]*dv10; 

  double p1_over_gamma[20] = {0.0}; 
  p1_over_gamma[0] = 1.732050807568877*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.732050807568877*gamma[4]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[8]*dv11; 
  p1_over_gamma[3] = 1.732050807568877*gamma[6]*dv11; 
  p1_over_gamma[4] = 3.872983346207417*gamma[12]*dv11; 
  p1_over_gamma[5] = 1.732050807568877*gamma[10]*dv11; 
  p1_over_gamma[6] = 3.872983346207417*gamma[14]*dv11; 
  p1_over_gamma[7] = 1.732050807568877*gamma[11]*dv11; 
  p1_over_gamma[9] = 1.732050807568877*gamma[16]*dv11; 
  p1_over_gamma[10] = 3.872983346207417*gamma[18]*dv11; 
  p1_over_gamma[13] = 1.732050807568877*gamma[17]*dv11; 
  p1_over_gamma[15] = 1.732050807568877*gamma[19]*dv11; 

  double p2_over_gamma[20] = {0.0}; 
  p2_over_gamma[0] = 1.732050807568877*gamma[3]*dv12; 
  p2_over_gamma[1] = 1.732050807568877*gamma[5]*dv12; 
  p2_over_gamma[2] = 1.732050807568877*gamma[6]*dv12; 
  p2_over_gamma[3] = 3.872983346207417*gamma[9]*dv12; 
  p2_over_gamma[4] = 1.732050807568877*gamma[10]*dv12; 
  p2_over_gamma[5] = 3.872983346207417*gamma[15]*dv12; 
  p2_over_gamma[6] = 3.872983346207417*gamma[16]*dv12; 
  p2_over_gamma[7] = 1.732050807568877*gamma[13]*dv12; 
  p2_over_gamma[8] = 1.732050807568877*gamma[14]*dv12; 
  p2_over_gamma[10] = 3.872983346207417*gamma[19]*dv12; 
  p2_over_gamma[11] = 1.732050807568877*gamma[17]*dv12; 
  p2_over_gamma[12] = 1.732050807568877*gamma[18]*dv12; 

  out[0] += (p0_over_gamma[16]*f[35]+p0_over_gamma[9]*f[32]+p0_over_gamma[14]*f[27]+p0_over_gamma[8]*f[24]+p0_over_gamma[10]*f[14]+p0_over_gamma[6]*f[10]+p0_over_gamma[5]*f[9]+p0_over_gamma[4]*f[7]+p0_over_gamma[3]*f[4]+p0_over_gamma[2]*f[3]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (1.0*p0_over_gamma[16]*f[37]+1.0*p0_over_gamma[9]*f[33]+1.0*p0_over_gamma[14]*f[29]+1.0*p0_over_gamma[8]*f[25]+p0_over_gamma[10]*f[15]+p0_over_gamma[6]*f[13]+p0_over_gamma[5]*f[12]+p0_over_gamma[4]*f[11]+p0_over_gamma[3]*f[8]+p0_over_gamma[2]*f[6]+p0_over_gamma[1]*f[5]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (p1_over_gamma[15]*f[34]+p1_over_gamma[9]*f[32]+p1_over_gamma[13]*f[19]+p1_over_gamma[7]*f[16]+p1_over_gamma[10]*f[14]+p1_over_gamma[6]*f[10]+p1_over_gamma[5]*f[9]+p1_over_gamma[4]*f[7]+p1_over_gamma[3]*f[4]+p1_over_gamma[2]*f[3]+p1_over_gamma[1]*f[2]+f[0]*p1_over_gamma[0])*volFact; 
  out[3] += (1.0*p1_over_gamma[15]*f[36]+1.0*p1_over_gamma[9]*f[33]+1.0*p1_over_gamma[13]*f[21]+1.0*p1_over_gamma[7]*f[17]+p1_over_gamma[10]*f[15]+p1_over_gamma[6]*f[13]+p1_over_gamma[5]*f[12]+p1_over_gamma[4]*f[11]+p1_over_gamma[3]*f[8]+p1_over_gamma[2]*f[6]+p1_over_gamma[1]*f[5]+p1_over_gamma[0]*f[1])*volFact; 
  out[4] += (p2_over_gamma[12]*f[26]+p2_over_gamma[8]*f[24]+p2_over_gamma[11]*f[18]+p2_over_gamma[7]*f[16]+p2_over_gamma[10]*f[14]+p2_over_gamma[6]*f[10]+p2_over_gamma[5]*f[9]+p2_over_gamma[4]*f[7]+p2_over_gamma[3]*f[4]+p2_over_gamma[2]*f[3]+p2_over_gamma[1]*f[2]+f[0]*p2_over_gamma[0])*volFact; 
  out[5] += (1.0*p2_over_gamma[12]*f[28]+1.0*p2_over_gamma[8]*f[25]+1.0*p2_over_gamma[11]*f[20]+1.0*p2_over_gamma[7]*f[17]+p2_over_gamma[10]*f[15]+p2_over_gamma[6]*f[13]+p2_over_gamma[5]*f[12]+p2_over_gamma[4]*f[11]+p2_over_gamma[3]*f[8]+p2_over_gamma[2]*f[6]+p2_over_gamma[1]*f[5]+p2_over_gamma[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M2_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += (gamma[19]*f[38]+gamma[16]*f[35]+gamma[15]*f[34]+gamma[9]*f[32]+gamma[18]*f[30]+gamma[14]*f[27]+gamma[12]*f[26]+gamma[8]*f[24]+gamma[17]*f[22]+gamma[13]*f[19]+gamma[11]*f[18]+gamma[7]*f[16]+gamma[10]*f[14]+gamma[6]*f[10]+gamma[5]*f[9]+gamma[4]*f[7]+gamma[3]*f[4]+gamma[2]*f[3]+gamma[1]*f[2]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0*gamma[19]*f[39]+1.0*gamma[16]*f[37]+1.0*gamma[15]*f[36]+1.0*gamma[9]*f[33]+1.0*gamma[18]*f[31]+1.0*gamma[14]*f[29]+1.0*gamma[12]*f[28]+1.0*gamma[8]*f[25]+1.0*gamma[17]*f[23]+1.0*gamma[13]*f[21]+1.0*gamma[11]*f[20]+1.0*gamma[7]*f[17]+gamma[10]*f[15]+gamma[6]*f[13]+gamma[5]*f[12]+gamma[4]*f[11]+gamma[3]*f[8]+gamma[2]*f[6]+gamma[1]*f[5]+gamma[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M3i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  out[0] += volFact*(2.828427124746191*f[0]*wx1+0.8164965809277261*f[2]*dv1); 
  out[1] += volFact*(2.828427124746191*f[1]*wx1+0.8164965809277261*f[5]*dv1); 
  out[2] += volFact*(2.828427124746191*f[0]*wx2+0.8164965809277261*f[3]*dv2); 
  out[3] += volFact*(2.828427124746191*f[1]*wx2+0.8164965809277261*f[6]*dv2); 
  out[4] += volFact*(2.828427124746191*f[0]*wx3+0.8164965809277261*f[4]*dv3); 
  out[5] += volFact*(2.828427124746191*f[1]*wx3+0.8164965809277261*f[8]*dv3); 
} 
GKYL_CU_DH void vlasov_sr_Ni_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  double p0_over_gamma[20] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[7]*dv10; 
  p0_over_gamma[2] = 1.732050807568877*gamma[4]*dv10; 
  p0_over_gamma[3] = 1.732050807568877*gamma[5]*dv10; 
  p0_over_gamma[4] = 3.872983346207417*gamma[11]*dv10; 
  p0_over_gamma[5] = 3.872983346207417*gamma[13]*dv10; 
  p0_over_gamma[6] = 1.732050807568877*gamma[10]*dv10; 
  p0_over_gamma[8] = 1.732050807568877*gamma[12]*dv10; 
  p0_over_gamma[9] = 1.732050807568877*gamma[15]*dv10; 
  p0_over_gamma[10] = 3.872983346207417*gamma[17]*dv10; 
  p0_over_gamma[14] = 1.732050807568877*gamma[18]*dv10; 
  p0_over_gamma[16] = 1.732050807568877*gamma[19]*dv10; 

  double p1_over_gamma[20] = {0.0}; 
  p1_over_gamma[0] = 1.732050807568877*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.732050807568877*gamma[4]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[8]*dv11; 
  p1_over_gamma[3] = 1.732050807568877*gamma[6]*dv11; 
  p1_over_gamma[4] = 3.872983346207417*gamma[12]*dv11; 
  p1_over_gamma[5] = 1.732050807568877*gamma[10]*dv11; 
  p1_over_gamma[6] = 3.872983346207417*gamma[14]*dv11; 
  p1_over_gamma[7] = 1.732050807568877*gamma[11]*dv11; 
  p1_over_gamma[9] = 1.732050807568877*gamma[16]*dv11; 
  p1_over_gamma[10] = 3.872983346207417*gamma[18]*dv11; 
  p1_over_gamma[13] = 1.732050807568877*gamma[17]*dv11; 
  p1_over_gamma[15] = 1.732050807568877*gamma[19]*dv11; 

  double p2_over_gamma[20] = {0.0}; 
  p2_over_gamma[0] = 1.732050807568877*gamma[3]*dv12; 
  p2_over_gamma[1] = 1.732050807568877*gamma[5]*dv12; 
  p2_over_gamma[2] = 1.732050807568877*gamma[6]*dv12; 
  p2_over_gamma[3] = 3.872983346207417*gamma[9]*dv12; 
  p2_over_gamma[4] = 1.732050807568877*gamma[10]*dv12; 
  p2_over_gamma[5] = 3.872983346207417*gamma[15]*dv12; 
  p2_over_gamma[6] = 3.872983346207417*gamma[16]*dv12; 
  p2_over_gamma[7] = 1.732050807568877*gamma[13]*dv12; 
  p2_over_gamma[8] = 1.732050807568877*gamma[14]*dv12; 
  p2_over_gamma[10] = 3.872983346207417*gamma[19]*dv12; 
  p2_over_gamma[11] = 1.732050807568877*gamma[17]*dv12; 
  p2_over_gamma[12] = 1.732050807568877*gamma[18]*dv12; 

  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += (p0_over_gamma[16]*f[35]+p0_over_gamma[9]*f[32]+p0_over_gamma[14]*f[27]+p0_over_gamma[8]*f[24]+p0_over_gamma[10]*f[14]+p0_over_gamma[6]*f[10]+p0_over_gamma[5]*f[9]+p0_over_gamma[4]*f[7]+p0_over_gamma[3]*f[4]+p0_over_gamma[2]*f[3]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[3] += (1.0*p0_over_gamma[16]*f[37]+1.0*p0_over_gamma[9]*f[33]+1.0*p0_over_gamma[14]*f[29]+1.0*p0_over_gamma[8]*f[25]+p0_over_gamma[10]*f[15]+p0_over_gamma[6]*f[13]+p0_over_gamma[5]*f[12]+p0_over_gamma[4]*f[11]+p0_over_gamma[3]*f[8]+p0_over_gamma[2]*f[6]+p0_over_gamma[1]*f[5]+p0_over_gamma[0]*f[1])*volFact; 
  out[4] += (p1_over_gamma[15]*f[34]+p1_over_gamma[9]*f[32]+p1_over_gamma[13]*f[19]+p1_over_gamma[7]*f[16]+p1_over_gamma[10]*f[14]+p1_over_gamma[6]*f[10]+p1_over_gamma[5]*f[9]+p1_over_gamma[4]*f[7]+p1_over_gamma[3]*f[4]+p1_over_gamma[2]*f[3]+p1_over_gamma[1]*f[2]+f[0]*p1_over_gamma[0])*volFact; 
  out[5] += (1.0*p1_over_gamma[15]*f[36]+1.0*p1_over_gamma[9]*f[33]+1.0*p1_over_gamma[13]*f[21]+1.0*p1_over_gamma[7]*f[17]+p1_over_gamma[10]*f[15]+p1_over_gamma[6]*f[13]+p1_over_gamma[5]*f[12]+p1_over_gamma[4]*f[11]+p1_over_gamma[3]*f[8]+p1_over_gamma[2]*f[6]+p1_over_gamma[1]*f[5]+p1_over_gamma[0]*f[1])*volFact; 
  out[6] += (p2_over_gamma[12]*f[26]+p2_over_gamma[8]*f[24]+p2_over_gamma[11]*f[18]+p2_over_gamma[7]*f[16]+p2_over_gamma[10]*f[14]+p2_over_gamma[6]*f[10]+p2_over_gamma[5]*f[9]+p2_over_gamma[4]*f[7]+p2_over_gamma[3]*f[4]+p2_over_gamma[2]*f[3]+p2_over_gamma[1]*f[2]+f[0]*p2_over_gamma[0])*volFact; 
  out[7] += (1.0*p2_over_gamma[12]*f[28]+1.0*p2_over_gamma[8]*f[25]+1.0*p2_over_gamma[11]*f[20]+1.0*p2_over_gamma[7]*f[17]+p2_over_gamma[10]*f[15]+p2_over_gamma[6]*f[13]+p2_over_gamma[5]*f[12]+p2_over_gamma[4]*f[11]+p2_over_gamma[3]*f[8]+p2_over_gamma[2]*f[6]+p2_over_gamma[1]*f[5]+p2_over_gamma[0]*f[1])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_Tij_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double dv10 = 2.0/dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double dv11 = 2.0/dxv[2]; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double dv12 = 2.0/dxv[3]; 
  double p0_over_gamma[20] = {0.0}; 
  p0_over_gamma[0] = 1.732050807568877*gamma[1]*dv10; 
  p0_over_gamma[1] = 3.872983346207417*gamma[7]*dv10; 
  p0_over_gamma[2] = 1.732050807568877*gamma[4]*dv10; 
  p0_over_gamma[3] = 1.732050807568877*gamma[5]*dv10; 
  p0_over_gamma[4] = 3.872983346207417*gamma[11]*dv10; 
  p0_over_gamma[5] = 3.872983346207417*gamma[13]*dv10; 
  p0_over_gamma[6] = 1.732050807568877*gamma[10]*dv10; 
  p0_over_gamma[8] = 1.732050807568877*gamma[12]*dv10; 
  p0_over_gamma[9] = 1.732050807568877*gamma[15]*dv10; 
  p0_over_gamma[10] = 3.872983346207417*gamma[17]*dv10; 
  p0_over_gamma[14] = 1.732050807568877*gamma[18]*dv10; 
  p0_over_gamma[16] = 1.732050807568877*gamma[19]*dv10; 

  double p1_over_gamma[20] = {0.0}; 
  p1_over_gamma[0] = 1.732050807568877*gamma[2]*dv11; 
  p1_over_gamma[1] = 1.732050807568877*gamma[4]*dv11; 
  p1_over_gamma[2] = 3.872983346207417*gamma[8]*dv11; 
  p1_over_gamma[3] = 1.732050807568877*gamma[6]*dv11; 
  p1_over_gamma[4] = 3.872983346207417*gamma[12]*dv11; 
  p1_over_gamma[5] = 1.732050807568877*gamma[10]*dv11; 
  p1_over_gamma[6] = 3.872983346207417*gamma[14]*dv11; 
  p1_over_gamma[7] = 1.732050807568877*gamma[11]*dv11; 
  p1_over_gamma[9] = 1.732050807568877*gamma[16]*dv11; 
  p1_over_gamma[10] = 3.872983346207417*gamma[18]*dv11; 
  p1_over_gamma[13] = 1.732050807568877*gamma[17]*dv11; 
  p1_over_gamma[15] = 1.732050807568877*gamma[19]*dv11; 

  double p2_over_gamma[20] = {0.0}; 
  p2_over_gamma[0] = 1.732050807568877*gamma[3]*dv12; 
  p2_over_gamma[1] = 1.732050807568877*gamma[5]*dv12; 
  p2_over_gamma[2] = 1.732050807568877*gamma[6]*dv12; 
  p2_over_gamma[3] = 3.872983346207417*gamma[9]*dv12; 
  p2_over_gamma[4] = 1.732050807568877*gamma[10]*dv12; 
  p2_over_gamma[5] = 3.872983346207417*gamma[15]*dv12; 
  p2_over_gamma[6] = 3.872983346207417*gamma[16]*dv12; 
  p2_over_gamma[7] = 1.732050807568877*gamma[13]*dv12; 
  p2_over_gamma[8] = 1.732050807568877*gamma[14]*dv12; 
  p2_over_gamma[10] = 3.872983346207417*gamma[19]*dv12; 
  p2_over_gamma[11] = 1.732050807568877*gamma[17]*dv12; 
  p2_over_gamma[12] = 1.732050807568877*gamma[18]*dv12; 

  out[0] += (gamma[19]*f[38]+gamma[16]*f[35]+gamma[15]*f[34]+gamma[9]*f[32]+gamma[18]*f[30]+gamma[14]*f[27]+gamma[12]*f[26]+gamma[8]*f[24]+gamma[17]*f[22]+gamma[13]*f[19]+gamma[11]*f[18]+gamma[7]*f[16]+gamma[10]*f[14]+gamma[6]*f[10]+gamma[5]*f[9]+gamma[4]*f[7]+gamma[3]*f[4]+gamma[2]*f[3]+gamma[1]*f[2]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0*gamma[19]*f[39]+1.0*gamma[16]*f[37]+1.0*gamma[15]*f[36]+1.0*gamma[9]*f[33]+1.0*gamma[18]*f[31]+1.0*gamma[14]*f[29]+1.0*gamma[12]*f[28]+1.0*gamma[8]*f[25]+1.0*gamma[17]*f[23]+1.0*gamma[13]*f[21]+1.0*gamma[11]*f[20]+1.0*gamma[7]*f[17]+gamma[10]*f[15]+gamma[6]*f[13]+gamma[5]*f[12]+gamma[4]*f[11]+gamma[3]*f[8]+gamma[2]*f[6]+gamma[1]*f[5]+gamma[0]*f[1])*volFact; 
  out[2] += volFact*(2.828427124746191*f[0]*wx1+0.8164965809277261*f[2]*dv1); 
  out[3] += volFact*(2.828427124746191*f[1]*wx1+0.8164965809277261*f[5]*dv1); 
  out[4] += volFact*(2.828427124746191*f[0]*wx2+0.8164965809277261*f[3]*dv2); 
  out[5] += volFact*(2.828427124746191*f[1]*wx2+0.8164965809277261*f[6]*dv2); 
  out[6] += volFact*(2.828427124746191*f[0]*wx3+0.8164965809277261*f[4]*dv3); 
  out[7] += volFact*(2.828427124746191*f[1]*wx3+0.8164965809277261*f[8]*dv3); 
  out[8] += volFact*(p0_over_gamma[16]*f[35]*wx1+p0_over_gamma[9]*f[32]*wx1+p0_over_gamma[14]*f[27]*wx1+p0_over_gamma[8]*f[24]*wx1+p0_over_gamma[10]*f[14]*wx1+p0_over_gamma[6]*f[10]*wx1+p0_over_gamma[5]*f[9]*wx1+p0_over_gamma[4]*f[7]*wx1+p0_over_gamma[3]*f[4]*wx1+p0_over_gamma[2]*f[3]*wx1+p0_over_gamma[1]*f[2]*wx1+f[0]*p0_over_gamma[0]*wx1+0.2886751345948129*p0_over_gamma[16]*f[38]*dv1+0.2886751345948129*p0_over_gamma[9]*f[34]*dv1+0.2886751345948129*p0_over_gamma[14]*f[30]*dv1+0.2886751345948129*p0_over_gamma[8]*f[26]*dv1+0.2581988897471612*p0_over_gamma[10]*f[22]*dv1+0.2581988897471611*p0_over_gamma[5]*f[19]*dv1+0.2581988897471611*p0_over_gamma[4]*f[18]*dv1+0.2581988897471612*p0_over_gamma[1]*f[16]*dv1+0.2886751345948129*p0_over_gamma[6]*f[14]*dv1+0.2886751345948129*f[10]*p0_over_gamma[10]*dv1+0.2886751345948129*p0_over_gamma[3]*f[9]*dv1+0.2886751345948129*p0_over_gamma[2]*f[7]*dv1+0.2886751345948129*f[4]*p0_over_gamma[5]*dv1+0.2886751345948129*f[3]*p0_over_gamma[4]*dv1+0.2886751345948129*p0_over_gamma[0]*f[2]*dv1+0.2886751345948129*f[0]*p0_over_gamma[1]*dv1); 
  out[9] += volFact*(1.0*p0_over_gamma[16]*f[37]*wx1+1.0*p0_over_gamma[9]*f[33]*wx1+1.0*p0_over_gamma[14]*f[29]*wx1+1.0*p0_over_gamma[8]*f[25]*wx1+p0_over_gamma[10]*f[15]*wx1+p0_over_gamma[6]*f[13]*wx1+p0_over_gamma[5]*f[12]*wx1+p0_over_gamma[4]*f[11]*wx1+p0_over_gamma[3]*f[8]*wx1+p0_over_gamma[2]*f[6]*wx1+p0_over_gamma[1]*f[5]*wx1+p0_over_gamma[0]*f[1]*wx1+0.2886751345948129*p0_over_gamma[16]*f[39]*dv1+0.2886751345948129*p0_over_gamma[9]*f[36]*dv1+0.2886751345948129*p0_over_gamma[14]*f[31]*dv1+0.2886751345948129*p0_over_gamma[8]*f[28]*dv1+0.2581988897471611*p0_over_gamma[10]*f[23]*dv1+0.2581988897471612*p0_over_gamma[5]*f[21]*dv1+0.2581988897471612*p0_over_gamma[4]*f[20]*dv1+0.2581988897471611*p0_over_gamma[1]*f[17]*dv1+0.2886751345948129*p0_over_gamma[6]*f[15]*dv1+0.2886751345948129*p0_over_gamma[10]*f[13]*dv1+0.2886751345948129*p0_over_gamma[3]*f[12]*dv1+0.2886751345948129*p0_over_gamma[2]*f[11]*dv1+0.2886751345948129*p0_over_gamma[5]*f[8]*dv1+0.2886751345948129*p0_over_gamma[4]*f[6]*dv1+0.2886751345948129*p0_over_gamma[0]*f[5]*dv1+0.2886751345948129*f[1]*p0_over_gamma[1]*dv1); 
  out[10] += volFact*(p0_over_gamma[16]*f[35]*wx2+p0_over_gamma[9]*f[32]*wx2+p0_over_gamma[14]*f[27]*wx2+p0_over_gamma[8]*f[24]*wx2+p0_over_gamma[10]*f[14]*wx2+p0_over_gamma[6]*f[10]*wx2+p0_over_gamma[5]*f[9]*wx2+p0_over_gamma[4]*f[7]*wx2+p0_over_gamma[3]*f[4]*wx2+p0_over_gamma[2]*f[3]*wx2+p0_over_gamma[1]*f[2]*wx2+f[0]*p0_over_gamma[0]*wx2+0.2886751345948129*p0_over_gamma[9]*f[35]*dv2+0.2886751345948129*p0_over_gamma[16]*f[32]*dv2+0.2581988897471612*p0_over_gamma[10]*f[30]*dv2+0.2581988897471611*p0_over_gamma[6]*f[27]*dv2+0.2581988897471611*p0_over_gamma[4]*f[26]*dv2+0.2581988897471612*p0_over_gamma[2]*f[24]*dv2+0.2581988897471611*f[10]*p0_over_gamma[14]*dv2+0.2886751345948129*p0_over_gamma[5]*f[14]*dv2+0.2886751345948129*f[9]*p0_over_gamma[10]*dv2+0.2886751345948129*p0_over_gamma[3]*f[10]*dv2+0.2581988897471612*f[3]*p0_over_gamma[8]*dv2+0.2886751345948129*p0_over_gamma[1]*f[7]*dv2+0.2886751345948129*f[4]*p0_over_gamma[6]*dv2+0.2886751345948129*f[2]*p0_over_gamma[4]*dv2+0.2886751345948129*p0_over_gamma[0]*f[3]*dv2+0.2886751345948129*f[0]*p0_over_gamma[2]*dv2); 
  out[11] += volFact*(1.0*p0_over_gamma[16]*f[37]*wx2+1.0*p0_over_gamma[9]*f[33]*wx2+1.0*p0_over_gamma[14]*f[29]*wx2+1.0*p0_over_gamma[8]*f[25]*wx2+p0_over_gamma[10]*f[15]*wx2+p0_over_gamma[6]*f[13]*wx2+p0_over_gamma[5]*f[12]*wx2+p0_over_gamma[4]*f[11]*wx2+p0_over_gamma[3]*f[8]*wx2+p0_over_gamma[2]*f[6]*wx2+p0_over_gamma[1]*f[5]*wx2+p0_over_gamma[0]*f[1]*wx2+0.2886751345948129*p0_over_gamma[9]*f[37]*dv2+0.2886751345948129*p0_over_gamma[16]*f[33]*dv2+0.2581988897471611*p0_over_gamma[10]*f[31]*dv2+0.2581988897471612*p0_over_gamma[6]*f[29]*dv2+0.2581988897471612*p0_over_gamma[4]*f[28]*dv2+0.2581988897471611*p0_over_gamma[2]*f[25]*dv2+0.2886751345948129*p0_over_gamma[5]*f[15]*dv2+0.2581988897471611*f[13]*p0_over_gamma[14]*dv2+0.2886751345948129*p0_over_gamma[3]*f[13]*dv2+0.2886751345948129*p0_over_gamma[10]*f[12]*dv2+0.2886751345948129*p0_over_gamma[1]*f[11]*dv2+0.2581988897471612*f[6]*p0_over_gamma[8]*dv2+0.2886751345948129*p0_over_gamma[6]*f[8]*dv2+0.2886751345948129*p0_over_gamma[0]*f[6]*dv2+0.2886751345948129*p0_over_gamma[4]*f[5]*dv2+0.2886751345948129*f[1]*p0_over_gamma[2]*dv2); 
  out[12] += volFact*(p0_over_gamma[16]*f[35]*wx3+p0_over_gamma[9]*f[32]*wx3+p0_over_gamma[14]*f[27]*wx3+p0_over_gamma[8]*f[24]*wx3+p0_over_gamma[10]*f[14]*wx3+p0_over_gamma[6]*f[10]*wx3+p0_over_gamma[5]*f[9]*wx3+p0_over_gamma[4]*f[7]*wx3+p0_over_gamma[3]*f[4]*wx3+p0_over_gamma[2]*f[3]*wx3+p0_over_gamma[1]*f[2]*wx3+f[0]*p0_over_gamma[0]*wx3+0.2581988897471612*p0_over_gamma[10]*f[38]*dv3+0.2581988897471611*p0_over_gamma[6]*f[35]*dv3+0.2581988897471611*p0_over_gamma[5]*f[34]*dv3+0.2581988897471612*p0_over_gamma[3]*f[32]*dv3+0.2886751345948129*p0_over_gamma[8]*f[27]*dv3+0.2886751345948129*p0_over_gamma[14]*f[24]*dv3+0.2581988897471611*f[10]*p0_over_gamma[16]*dv3+0.2886751345948129*p0_over_gamma[4]*f[14]*dv3+0.2886751345948129*f[7]*p0_over_gamma[10]*dv3+0.2886751345948129*p0_over_gamma[2]*f[10]*dv3+0.2581988897471612*f[4]*p0_over_gamma[9]*dv3+0.2886751345948129*p0_over_gamma[1]*f[9]*dv3+0.2886751345948129*f[3]*p0_over_gamma[6]*dv3+0.2886751345948129*f[2]*p0_over_gamma[5]*dv3+0.2886751345948129*p0_over_gamma[0]*f[4]*dv3+0.2886751345948129*f[0]*p0_over_gamma[3]*dv3); 
  out[13] += volFact*(1.0*p0_over_gamma[16]*f[37]*wx3+1.0*p0_over_gamma[9]*f[33]*wx3+1.0*p0_over_gamma[14]*f[29]*wx3+1.0*p0_over_gamma[8]*f[25]*wx3+p0_over_gamma[10]*f[15]*wx3+p0_over_gamma[6]*f[13]*wx3+p0_over_gamma[5]*f[12]*wx3+p0_over_gamma[4]*f[11]*wx3+p0_over_gamma[3]*f[8]*wx3+p0_over_gamma[2]*f[6]*wx3+p0_over_gamma[1]*f[5]*wx3+p0_over_gamma[0]*f[1]*wx3+0.2581988897471611*p0_over_gamma[10]*f[39]*dv3+0.2581988897471612*p0_over_gamma[6]*f[37]*dv3+0.2581988897471612*p0_over_gamma[5]*f[36]*dv3+0.2581988897471611*p0_over_gamma[3]*f[33]*dv3+0.2886751345948129*p0_over_gamma[8]*f[29]*dv3+0.2886751345948129*p0_over_gamma[14]*f[25]*dv3+0.2581988897471611*f[13]*p0_over_gamma[16]*dv3+0.2886751345948129*p0_over_gamma[4]*f[15]*dv3+0.2886751345948129*p0_over_gamma[2]*f[13]*dv3+0.2886751345948129*p0_over_gamma[1]*f[12]*dv3+0.2886751345948129*p0_over_gamma[10]*f[11]*dv3+0.2581988897471612*f[8]*p0_over_gamma[9]*dv3+0.2886751345948129*p0_over_gamma[0]*f[8]*dv3+0.2886751345948129*f[6]*p0_over_gamma[6]*dv3+0.2886751345948129*f[5]*p0_over_gamma[5]*dv3+0.2886751345948129*f[1]*p0_over_gamma[3]*dv3); 
  out[14] += volFact*(p1_over_gamma[15]*f[34]*wx2+p1_over_gamma[9]*f[32]*wx2+p1_over_gamma[13]*f[19]*wx2+p1_over_gamma[7]*f[16]*wx2+p1_over_gamma[10]*f[14]*wx2+p1_over_gamma[6]*f[10]*wx2+p1_over_gamma[5]*f[9]*wx2+p1_over_gamma[4]*f[7]*wx2+p1_over_gamma[3]*f[4]*wx2+p1_over_gamma[2]*f[3]*wx2+p1_over_gamma[1]*f[2]*wx2+f[0]*p1_over_gamma[0]*wx2+0.2886751345948129*p1_over_gamma[15]*f[38]*dv2+0.2886751345948129*p1_over_gamma[9]*f[35]*dv2+0.2581988897471612*p1_over_gamma[10]*f[30]*dv2+0.2581988897471611*p1_over_gamma[6]*f[27]*dv2+0.2581988897471611*p1_over_gamma[4]*f[26]*dv2+0.2581988897471612*p1_over_gamma[2]*f[24]*dv2+0.2886751345948129*p1_over_gamma[13]*f[22]*dv2+0.2886751345948129*p1_over_gamma[7]*f[18]*dv2+0.2886751345948129*p1_over_gamma[5]*f[14]*dv2+0.2886751345948129*f[9]*p1_over_gamma[10]*dv2+0.2886751345948129*p1_over_gamma[3]*f[10]*dv2+0.2886751345948129*p1_over_gamma[1]*f[7]*dv2+0.2886751345948129*f[4]*p1_over_gamma[6]*dv2+0.2886751345948129*f[2]*p1_over_gamma[4]*dv2+0.2886751345948129*p1_over_gamma[0]*f[3]*dv2+0.2886751345948129*f[0]*p1_over_gamma[2]*dv2); 
  out[15] += volFact*(1.0*p1_over_gamma[15]*f[36]*wx2+1.0*p1_over_gamma[9]*f[33]*wx2+1.0*p1_over_gamma[13]*f[21]*wx2+1.0*p1_over_gamma[7]*f[17]*wx2+p1_over_gamma[10]*f[15]*wx2+p1_over_gamma[6]*f[13]*wx2+p1_over_gamma[5]*f[12]*wx2+p1_over_gamma[4]*f[11]*wx2+p1_over_gamma[3]*f[8]*wx2+p1_over_gamma[2]*f[6]*wx2+p1_over_gamma[1]*f[5]*wx2+p1_over_gamma[0]*f[1]*wx2+0.2886751345948129*p1_over_gamma[15]*f[39]*dv2+0.2886751345948129*p1_over_gamma[9]*f[37]*dv2+0.2581988897471611*p1_over_gamma[10]*f[31]*dv2+0.2581988897471612*p1_over_gamma[6]*f[29]*dv2+0.2581988897471612*p1_over_gamma[4]*f[28]*dv2+0.2581988897471611*p1_over_gamma[2]*f[25]*dv2+0.2886751345948129*p1_over_gamma[13]*f[23]*dv2+0.2886751345948129*p1_over_gamma[7]*f[20]*dv2+0.2886751345948129*p1_over_gamma[5]*f[15]*dv2+0.2886751345948129*p1_over_gamma[3]*f[13]*dv2+0.2886751345948129*p1_over_gamma[10]*f[12]*dv2+0.2886751345948129*p1_over_gamma[1]*f[11]*dv2+0.2886751345948129*p1_over_gamma[6]*f[8]*dv2+0.2886751345948129*p1_over_gamma[0]*f[6]*dv2+0.2886751345948129*p1_over_gamma[4]*f[5]*dv2+0.2886751345948129*f[1]*p1_over_gamma[2]*dv2); 
  out[16] += volFact*(p1_over_gamma[15]*f[34]*wx3+p1_over_gamma[9]*f[32]*wx3+p1_over_gamma[13]*f[19]*wx3+p1_over_gamma[7]*f[16]*wx3+p1_over_gamma[10]*f[14]*wx3+p1_over_gamma[6]*f[10]*wx3+p1_over_gamma[5]*f[9]*wx3+p1_over_gamma[4]*f[7]*wx3+p1_over_gamma[3]*f[4]*wx3+p1_over_gamma[2]*f[3]*wx3+p1_over_gamma[1]*f[2]*wx3+f[0]*p1_over_gamma[0]*wx3+0.2581988897471612*p1_over_gamma[10]*f[38]*dv3+0.2581988897471611*p1_over_gamma[6]*f[35]*dv3+0.2581988897471611*p1_over_gamma[5]*f[34]*dv3+0.2581988897471612*p1_over_gamma[3]*f[32]*dv3+0.2886751345948129*p1_over_gamma[7]*f[19]*dv3+0.2886751345948129*p1_over_gamma[13]*f[16]*dv3+0.2581988897471611*f[9]*p1_over_gamma[15]*dv3+0.2886751345948129*p1_over_gamma[4]*f[14]*dv3+0.2886751345948129*f[7]*p1_over_gamma[10]*dv3+0.2886751345948129*p1_over_gamma[2]*f[10]*dv3+0.2581988897471612*f[4]*p1_over_gamma[9]*dv3+0.2886751345948129*p1_over_gamma[1]*f[9]*dv3+0.2886751345948129*f[3]*p1_over_gamma[6]*dv3+0.2886751345948129*f[2]*p1_over_gamma[5]*dv3+0.2886751345948129*p1_over_gamma[0]*f[4]*dv3+0.2886751345948129*f[0]*p1_over_gamma[3]*dv3); 
  out[17] += volFact*(1.0*p1_over_gamma[15]*f[36]*wx3+1.0*p1_over_gamma[9]*f[33]*wx3+1.0*p1_over_gamma[13]*f[21]*wx3+1.0*p1_over_gamma[7]*f[17]*wx3+p1_over_gamma[10]*f[15]*wx3+p1_over_gamma[6]*f[13]*wx3+p1_over_gamma[5]*f[12]*wx3+p1_over_gamma[4]*f[11]*wx3+p1_over_gamma[3]*f[8]*wx3+p1_over_gamma[2]*f[6]*wx3+p1_over_gamma[1]*f[5]*wx3+p1_over_gamma[0]*f[1]*wx3+0.2581988897471611*p1_over_gamma[10]*f[39]*dv3+0.2581988897471612*p1_over_gamma[6]*f[37]*dv3+0.2581988897471612*p1_over_gamma[5]*f[36]*dv3+0.2581988897471611*p1_over_gamma[3]*f[33]*dv3+0.2886751345948129*p1_over_gamma[7]*f[21]*dv3+0.2886751345948129*p1_over_gamma[13]*f[17]*dv3+0.2581988897471611*f[12]*p1_over_gamma[15]*dv3+0.2886751345948129*p1_over_gamma[4]*f[15]*dv3+0.2886751345948129*p1_over_gamma[2]*f[13]*dv3+0.2886751345948129*p1_over_gamma[1]*f[12]*dv3+0.2886751345948129*p1_over_gamma[10]*f[11]*dv3+0.2581988897471612*f[8]*p1_over_gamma[9]*dv3+0.2886751345948129*p1_over_gamma[0]*f[8]*dv3+0.2886751345948129*f[6]*p1_over_gamma[6]*dv3+0.2886751345948129*f[5]*p1_over_gamma[5]*dv3+0.2886751345948129*f[1]*p1_over_gamma[3]*dv3); 
  out[18] += volFact*(p2_over_gamma[12]*f[26]*wx3+p2_over_gamma[8]*f[24]*wx3+p2_over_gamma[11]*f[18]*wx3+p2_over_gamma[7]*f[16]*wx3+p2_over_gamma[10]*f[14]*wx3+p2_over_gamma[6]*f[10]*wx3+p2_over_gamma[5]*f[9]*wx3+p2_over_gamma[4]*f[7]*wx3+p2_over_gamma[3]*f[4]*wx3+p2_over_gamma[2]*f[3]*wx3+p2_over_gamma[1]*f[2]*wx3+f[0]*p2_over_gamma[0]*wx3+0.2581988897471612*p2_over_gamma[10]*f[38]*dv3+0.2581988897471611*p2_over_gamma[6]*f[35]*dv3+0.2581988897471611*p2_over_gamma[5]*f[34]*dv3+0.2581988897471612*p2_over_gamma[3]*f[32]*dv3+0.2886751345948129*p2_over_gamma[12]*f[30]*dv3+0.2886751345948129*p2_over_gamma[8]*f[27]*dv3+0.2886751345948129*p2_over_gamma[11]*f[22]*dv3+0.2886751345948129*p2_over_gamma[7]*f[19]*dv3+0.2886751345948129*p2_over_gamma[4]*f[14]*dv3+0.2886751345948129*f[7]*p2_over_gamma[10]*dv3+0.2886751345948129*p2_over_gamma[2]*f[10]*dv3+0.2886751345948129*p2_over_gamma[1]*f[9]*dv3+0.2886751345948129*f[3]*p2_over_gamma[6]*dv3+0.2886751345948129*f[2]*p2_over_gamma[5]*dv3+0.2886751345948129*p2_over_gamma[0]*f[4]*dv3+0.2886751345948129*f[0]*p2_over_gamma[3]*dv3); 
  out[19] += volFact*(1.0*p2_over_gamma[12]*f[28]*wx3+1.0*p2_over_gamma[8]*f[25]*wx3+1.0*p2_over_gamma[11]*f[20]*wx3+1.0*p2_over_gamma[7]*f[17]*wx3+p2_over_gamma[10]*f[15]*wx3+p2_over_gamma[6]*f[13]*wx3+p2_over_gamma[5]*f[12]*wx3+p2_over_gamma[4]*f[11]*wx3+p2_over_gamma[3]*f[8]*wx3+p2_over_gamma[2]*f[6]*wx3+p2_over_gamma[1]*f[5]*wx3+p2_over_gamma[0]*f[1]*wx3+0.2581988897471611*p2_over_gamma[10]*f[39]*dv3+0.2581988897471612*p2_over_gamma[6]*f[37]*dv3+0.2581988897471612*p2_over_gamma[5]*f[36]*dv3+0.2581988897471611*p2_over_gamma[3]*f[33]*dv3+0.2886751345948129*p2_over_gamma[12]*f[31]*dv3+0.2886751345948129*p2_over_gamma[8]*f[29]*dv3+0.2886751345948129*p2_over_gamma[11]*f[23]*dv3+0.2886751345948129*p2_over_gamma[7]*f[21]*dv3+0.2886751345948129*p2_over_gamma[4]*f[15]*dv3+0.2886751345948129*p2_over_gamma[2]*f[13]*dv3+0.2886751345948129*p2_over_gamma[1]*f[12]*dv3+0.2886751345948129*p2_over_gamma[10]*f[11]*dv3+0.2886751345948129*p2_over_gamma[0]*f[8]*dv3+0.2886751345948129*f[6]*p2_over_gamma[6]*dv3+0.2886751345948129*f[5]*p2_over_gamma[5]*dv3+0.2886751345948129*f[1]*p2_over_gamma[3]*dv3); 
} 
GKYL_CU_DH void vlasov_sr_int_mom_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (1.414213562373095*gamma[19]*f[38]+1.414213562373095*gamma[16]*f[35]+1.414213562373095*gamma[15]*f[34]+1.414213562373095*gamma[9]*f[32]+1.414213562373095*gamma[18]*f[30]+1.414213562373095*gamma[14]*f[27]+1.414213562373095*gamma[12]*f[26]+1.414213562373095*gamma[8]*f[24]+1.414213562373095*gamma[17]*f[22]+1.414213562373095*gamma[13]*f[19]+1.414213562373095*gamma[11]*f[18]+1.414213562373095*gamma[7]*f[16]+1.414213562373095*gamma[10]*f[14]+1.414213562373095*gamma[6]*f[10]+1.414213562373095*gamma[5]*f[9]+1.414213562373095*gamma[4]*f[7]+1.414213562373095*gamma[3]*f[4]+1.414213562373095*gamma[2]*f[3]+1.414213562373095*gamma[1]*f[2]+1.414213562373095*f[0]*gamma[0])*volFact; 
  out[2] += volFact*(4.0*f[0]*wx1+1.154700538379252*f[2]*dv1); 
  out[3] += volFact*(4.0*f[0]*wx2+1.154700538379252*f[3]*dv2); 
  out[4] += volFact*(4.0*f[0]*wx3+1.154700538379252*f[4]*dv3); 
} 
