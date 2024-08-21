#include <gkyl_mom_vlasov_sr_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[6]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 
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

  out[0] += (p0_over_gamma[16]*f[68]+p0_over_gamma[9]*f[64]+p0_over_gamma[14]*f[52]+p0_over_gamma[8]*f[48]+p0_over_gamma[10]*f[25]+p0_over_gamma[6]*f[15]+p0_over_gamma[5]*f[14]+p0_over_gamma[4]*f[11]+p0_over_gamma[3]*f[5]+p0_over_gamma[2]*f[4]+p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (1.0*p0_over_gamma[16]*f[72]+1.0*p0_over_gamma[9]*f[65]+1.0*p0_over_gamma[14]*f[56]+1.0*p0_over_gamma[8]*f[49]+p0_over_gamma[10]*f[29]+p0_over_gamma[6]*f[23]+p0_over_gamma[5]*f[21]+p0_over_gamma[4]*f[18]+p0_over_gamma[3]*f[12]+p0_over_gamma[2]*f[9]+p0_over_gamma[1]*f[7]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (1.0*p0_over_gamma[16]*f[73]+1.0*p0_over_gamma[9]*f[66]+1.0*p0_over_gamma[14]*f[57]+1.0*p0_over_gamma[8]*f[50]+p0_over_gamma[10]*f[30]+p0_over_gamma[6]*f[24]+p0_over_gamma[5]*f[22]+p0_over_gamma[4]*f[19]+p0_over_gamma[3]*f[13]+p0_over_gamma[2]*f[10]+p0_over_gamma[1]*f[8]+p0_over_gamma[0]*f[2])*volFact; 
  out[3] += (p0_over_gamma[16]*f[76]+p0_over_gamma[9]*f[69]+p0_over_gamma[14]*f[60]+p0_over_gamma[8]*f[53]+p0_over_gamma[10]*f[31]+p0_over_gamma[6]*f[28]+p0_over_gamma[5]*f[27]+p0_over_gamma[4]*f[26]+p0_over_gamma[3]*f[20]+p0_over_gamma[2]*f[17]+p0_over_gamma[1]*f[16]+p0_over_gamma[0]*f[6])*volFact; 
  out[4] += (p1_over_gamma[15]*f[67]+p1_over_gamma[9]*f[64]+p1_over_gamma[13]*f[36]+p1_over_gamma[7]*f[32]+p1_over_gamma[10]*f[25]+p1_over_gamma[6]*f[15]+p1_over_gamma[5]*f[14]+p1_over_gamma[4]*f[11]+p1_over_gamma[3]*f[5]+p1_over_gamma[2]*f[4]+p1_over_gamma[1]*f[3]+f[0]*p1_over_gamma[0])*volFact; 
  out[5] += (1.0*p1_over_gamma[15]*f[70]+1.0*p1_over_gamma[9]*f[65]+1.0*p1_over_gamma[13]*f[40]+1.0*p1_over_gamma[7]*f[33]+p1_over_gamma[10]*f[29]+p1_over_gamma[6]*f[23]+p1_over_gamma[5]*f[21]+p1_over_gamma[4]*f[18]+p1_over_gamma[3]*f[12]+p1_over_gamma[2]*f[9]+p1_over_gamma[1]*f[7]+p1_over_gamma[0]*f[1])*volFact; 
  out[6] += (1.0*p1_over_gamma[15]*f[71]+1.0*p1_over_gamma[9]*f[66]+1.0*p1_over_gamma[13]*f[41]+1.0*p1_over_gamma[7]*f[34]+p1_over_gamma[10]*f[30]+p1_over_gamma[6]*f[24]+p1_over_gamma[5]*f[22]+p1_over_gamma[4]*f[19]+p1_over_gamma[3]*f[13]+p1_over_gamma[2]*f[10]+p1_over_gamma[1]*f[8]+p1_over_gamma[0]*f[2])*volFact; 
  out[7] += (p1_over_gamma[15]*f[75]+p1_over_gamma[9]*f[69]+p1_over_gamma[13]*f[44]+p1_over_gamma[7]*f[37]+p1_over_gamma[10]*f[31]+p1_over_gamma[6]*f[28]+p1_over_gamma[5]*f[27]+p1_over_gamma[4]*f[26]+p1_over_gamma[3]*f[20]+p1_over_gamma[2]*f[17]+p1_over_gamma[1]*f[16]+p1_over_gamma[0]*f[6])*volFact; 
  out[8] += (p2_over_gamma[12]*f[51]+p2_over_gamma[8]*f[48]+p2_over_gamma[11]*f[35]+p2_over_gamma[7]*f[32]+p2_over_gamma[10]*f[25]+p2_over_gamma[6]*f[15]+p2_over_gamma[5]*f[14]+p2_over_gamma[4]*f[11]+p2_over_gamma[3]*f[5]+p2_over_gamma[2]*f[4]+p2_over_gamma[1]*f[3]+f[0]*p2_over_gamma[0])*volFact; 
  out[9] += (1.0*p2_over_gamma[12]*f[54]+1.0*p2_over_gamma[8]*f[49]+1.0*p2_over_gamma[11]*f[38]+1.0*p2_over_gamma[7]*f[33]+p2_over_gamma[10]*f[29]+p2_over_gamma[6]*f[23]+p2_over_gamma[5]*f[21]+p2_over_gamma[4]*f[18]+p2_over_gamma[3]*f[12]+p2_over_gamma[2]*f[9]+p2_over_gamma[1]*f[7]+p2_over_gamma[0]*f[1])*volFact; 
  out[10] += (1.0*p2_over_gamma[12]*f[55]+1.0*p2_over_gamma[8]*f[50]+1.0*p2_over_gamma[11]*f[39]+1.0*p2_over_gamma[7]*f[34]+p2_over_gamma[10]*f[30]+p2_over_gamma[6]*f[24]+p2_over_gamma[5]*f[22]+p2_over_gamma[4]*f[19]+p2_over_gamma[3]*f[13]+p2_over_gamma[2]*f[10]+p2_over_gamma[1]*f[8]+p2_over_gamma[0]*f[2])*volFact; 
  out[11] += (p2_over_gamma[12]*f[59]+p2_over_gamma[8]*f[53]+p2_over_gamma[11]*f[43]+p2_over_gamma[7]*f[37]+p2_over_gamma[10]*f[31]+p2_over_gamma[6]*f[28]+p2_over_gamma[5]*f[27]+p2_over_gamma[4]*f[26]+p2_over_gamma[3]*f[20]+p2_over_gamma[2]*f[17]+p2_over_gamma[1]*f[16]+p2_over_gamma[0]*f[6])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M2_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  out[0] += (gamma[19]*f[74]+gamma[16]*f[68]+gamma[15]*f[67]+gamma[9]*f[64]+gamma[18]*f[58]+gamma[14]*f[52]+gamma[12]*f[51]+gamma[8]*f[48]+gamma[17]*f[42]+gamma[13]*f[36]+gamma[11]*f[35]+gamma[7]*f[32]+gamma[10]*f[25]+gamma[6]*f[15]+gamma[5]*f[14]+gamma[4]*f[11]+gamma[3]*f[5]+gamma[2]*f[4]+gamma[1]*f[3]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0*gamma[19]*f[77]+1.0*gamma[16]*f[72]+1.0*gamma[15]*f[70]+1.0*gamma[9]*f[65]+1.0*gamma[18]*f[61]+1.0*gamma[14]*f[56]+1.0*gamma[12]*f[54]+1.0*gamma[8]*f[49]+1.0*gamma[17]*f[45]+1.0*gamma[13]*f[40]+1.0*gamma[11]*f[38]+1.0*gamma[7]*f[33]+gamma[10]*f[29]+gamma[6]*f[23]+gamma[5]*f[21]+gamma[4]*f[18]+gamma[3]*f[12]+gamma[2]*f[9]+gamma[1]*f[7]+gamma[0]*f[1])*volFact; 
  out[2] += (1.0*gamma[19]*f[78]+1.0*gamma[16]*f[73]+1.0*gamma[15]*f[71]+1.0*gamma[9]*f[66]+1.0*gamma[18]*f[62]+1.0*gamma[14]*f[57]+1.0*gamma[12]*f[55]+1.0*gamma[8]*f[50]+1.0*gamma[17]*f[46]+1.0*gamma[13]*f[41]+1.0*gamma[11]*f[39]+1.0*gamma[7]*f[34]+gamma[10]*f[30]+gamma[6]*f[24]+gamma[5]*f[22]+gamma[4]*f[19]+gamma[3]*f[13]+gamma[2]*f[10]+gamma[1]*f[8]+gamma[0]*f[2])*volFact; 
  out[3] += (gamma[19]*f[79]+gamma[16]*f[76]+gamma[15]*f[75]+gamma[9]*f[69]+gamma[18]*f[63]+gamma[14]*f[60]+gamma[12]*f[59]+gamma[8]*f[53]+gamma[17]*f[47]+gamma[13]*f[44]+gamma[11]*f[43]+gamma[7]*f[37]+gamma[10]*f[31]+gamma[6]*f[28]+gamma[5]*f[27]+gamma[4]*f[26]+gamma[3]*f[20]+gamma[2]*f[17]+gamma[1]*f[16]+gamma[0]*f[6])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M3i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  out[0] += volFact*(2.828427124746191*f[0]*wx1+0.8164965809277261*f[3]*dv1); 
  out[1] += volFact*(2.828427124746191*f[1]*wx1+0.8164965809277261*f[7]*dv1); 
  out[2] += volFact*(2.828427124746191*f[2]*wx1+0.8164965809277261*f[8]*dv1); 
  out[3] += volFact*(2.828427124746191*f[6]*wx1+0.8164965809277261*f[16]*dv1); 
  out[4] += volFact*(2.828427124746191*f[0]*wx2+0.8164965809277261*f[4]*dv2); 
  out[5] += volFact*(2.828427124746191*f[1]*wx2+0.8164965809277261*f[9]*dv2); 
  out[6] += volFact*(2.828427124746191*f[2]*wx2+0.8164965809277261*f[10]*dv2); 
  out[7] += volFact*(2.828427124746191*f[6]*wx2+0.8164965809277261*f[17]*dv2); 
  out[8] += volFact*(2.828427124746191*f[0]*wx3+0.8164965809277261*f[5]*dv3); 
  out[9] += volFact*(2.828427124746191*f[1]*wx3+0.8164965809277261*f[12]*dv3); 
  out[10] += volFact*(2.828427124746191*f[2]*wx3+0.8164965809277261*f[13]*dv3); 
  out[11] += volFact*(2.828427124746191*f[6]*wx3+0.8164965809277261*f[20]*dv3); 
} 
GKYL_CU_DH void vlasov_sr_Ni_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 
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
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[6]*volFact; 
  out[4] += (p0_over_gamma[16]*f[68]+p0_over_gamma[9]*f[64]+p0_over_gamma[14]*f[52]+p0_over_gamma[8]*f[48]+p0_over_gamma[10]*f[25]+p0_over_gamma[6]*f[15]+p0_over_gamma[5]*f[14]+p0_over_gamma[4]*f[11]+p0_over_gamma[3]*f[5]+p0_over_gamma[2]*f[4]+p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*volFact; 
  out[5] += (1.0*p0_over_gamma[16]*f[72]+1.0*p0_over_gamma[9]*f[65]+1.0*p0_over_gamma[14]*f[56]+1.0*p0_over_gamma[8]*f[49]+p0_over_gamma[10]*f[29]+p0_over_gamma[6]*f[23]+p0_over_gamma[5]*f[21]+p0_over_gamma[4]*f[18]+p0_over_gamma[3]*f[12]+p0_over_gamma[2]*f[9]+p0_over_gamma[1]*f[7]+p0_over_gamma[0]*f[1])*volFact; 
  out[6] += (1.0*p0_over_gamma[16]*f[73]+1.0*p0_over_gamma[9]*f[66]+1.0*p0_over_gamma[14]*f[57]+1.0*p0_over_gamma[8]*f[50]+p0_over_gamma[10]*f[30]+p0_over_gamma[6]*f[24]+p0_over_gamma[5]*f[22]+p0_over_gamma[4]*f[19]+p0_over_gamma[3]*f[13]+p0_over_gamma[2]*f[10]+p0_over_gamma[1]*f[8]+p0_over_gamma[0]*f[2])*volFact; 
  out[7] += (p0_over_gamma[16]*f[76]+p0_over_gamma[9]*f[69]+p0_over_gamma[14]*f[60]+p0_over_gamma[8]*f[53]+p0_over_gamma[10]*f[31]+p0_over_gamma[6]*f[28]+p0_over_gamma[5]*f[27]+p0_over_gamma[4]*f[26]+p0_over_gamma[3]*f[20]+p0_over_gamma[2]*f[17]+p0_over_gamma[1]*f[16]+p0_over_gamma[0]*f[6])*volFact; 
  out[8] += (p1_over_gamma[15]*f[67]+p1_over_gamma[9]*f[64]+p1_over_gamma[13]*f[36]+p1_over_gamma[7]*f[32]+p1_over_gamma[10]*f[25]+p1_over_gamma[6]*f[15]+p1_over_gamma[5]*f[14]+p1_over_gamma[4]*f[11]+p1_over_gamma[3]*f[5]+p1_over_gamma[2]*f[4]+p1_over_gamma[1]*f[3]+f[0]*p1_over_gamma[0])*volFact; 
  out[9] += (1.0*p1_over_gamma[15]*f[70]+1.0*p1_over_gamma[9]*f[65]+1.0*p1_over_gamma[13]*f[40]+1.0*p1_over_gamma[7]*f[33]+p1_over_gamma[10]*f[29]+p1_over_gamma[6]*f[23]+p1_over_gamma[5]*f[21]+p1_over_gamma[4]*f[18]+p1_over_gamma[3]*f[12]+p1_over_gamma[2]*f[9]+p1_over_gamma[1]*f[7]+p1_over_gamma[0]*f[1])*volFact; 
  out[10] += (1.0*p1_over_gamma[15]*f[71]+1.0*p1_over_gamma[9]*f[66]+1.0*p1_over_gamma[13]*f[41]+1.0*p1_over_gamma[7]*f[34]+p1_over_gamma[10]*f[30]+p1_over_gamma[6]*f[24]+p1_over_gamma[5]*f[22]+p1_over_gamma[4]*f[19]+p1_over_gamma[3]*f[13]+p1_over_gamma[2]*f[10]+p1_over_gamma[1]*f[8]+p1_over_gamma[0]*f[2])*volFact; 
  out[11] += (p1_over_gamma[15]*f[75]+p1_over_gamma[9]*f[69]+p1_over_gamma[13]*f[44]+p1_over_gamma[7]*f[37]+p1_over_gamma[10]*f[31]+p1_over_gamma[6]*f[28]+p1_over_gamma[5]*f[27]+p1_over_gamma[4]*f[26]+p1_over_gamma[3]*f[20]+p1_over_gamma[2]*f[17]+p1_over_gamma[1]*f[16]+p1_over_gamma[0]*f[6])*volFact; 
  out[12] += (p2_over_gamma[12]*f[51]+p2_over_gamma[8]*f[48]+p2_over_gamma[11]*f[35]+p2_over_gamma[7]*f[32]+p2_over_gamma[10]*f[25]+p2_over_gamma[6]*f[15]+p2_over_gamma[5]*f[14]+p2_over_gamma[4]*f[11]+p2_over_gamma[3]*f[5]+p2_over_gamma[2]*f[4]+p2_over_gamma[1]*f[3]+f[0]*p2_over_gamma[0])*volFact; 
  out[13] += (1.0*p2_over_gamma[12]*f[54]+1.0*p2_over_gamma[8]*f[49]+1.0*p2_over_gamma[11]*f[38]+1.0*p2_over_gamma[7]*f[33]+p2_over_gamma[10]*f[29]+p2_over_gamma[6]*f[23]+p2_over_gamma[5]*f[21]+p2_over_gamma[4]*f[18]+p2_over_gamma[3]*f[12]+p2_over_gamma[2]*f[9]+p2_over_gamma[1]*f[7]+p2_over_gamma[0]*f[1])*volFact; 
  out[14] += (1.0*p2_over_gamma[12]*f[55]+1.0*p2_over_gamma[8]*f[50]+1.0*p2_over_gamma[11]*f[39]+1.0*p2_over_gamma[7]*f[34]+p2_over_gamma[10]*f[30]+p2_over_gamma[6]*f[24]+p2_over_gamma[5]*f[22]+p2_over_gamma[4]*f[19]+p2_over_gamma[3]*f[13]+p2_over_gamma[2]*f[10]+p2_over_gamma[1]*f[8]+p2_over_gamma[0]*f[2])*volFact; 
  out[15] += (p2_over_gamma[12]*f[59]+p2_over_gamma[8]*f[53]+p2_over_gamma[11]*f[43]+p2_over_gamma[7]*f[37]+p2_over_gamma[10]*f[31]+p2_over_gamma[6]*f[28]+p2_over_gamma[5]*f[27]+p2_over_gamma[4]*f[26]+p2_over_gamma[3]*f[20]+p2_over_gamma[2]*f[17]+p2_over_gamma[1]*f[16]+p2_over_gamma[0]*f[6])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_Tij_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double dv10 = 2.0/dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double dv11 = 2.0/dxv[3]; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double dv12 = 2.0/dxv[4]; 
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

  out[0] += (gamma[19]*f[74]+gamma[16]*f[68]+gamma[15]*f[67]+gamma[9]*f[64]+gamma[18]*f[58]+gamma[14]*f[52]+gamma[12]*f[51]+gamma[8]*f[48]+gamma[17]*f[42]+gamma[13]*f[36]+gamma[11]*f[35]+gamma[7]*f[32]+gamma[10]*f[25]+gamma[6]*f[15]+gamma[5]*f[14]+gamma[4]*f[11]+gamma[3]*f[5]+gamma[2]*f[4]+gamma[1]*f[3]+f[0]*gamma[0])*volFact; 
  out[1] += (1.0*gamma[19]*f[77]+1.0*gamma[16]*f[72]+1.0*gamma[15]*f[70]+1.0*gamma[9]*f[65]+1.0*gamma[18]*f[61]+1.0*gamma[14]*f[56]+1.0*gamma[12]*f[54]+1.0*gamma[8]*f[49]+1.0*gamma[17]*f[45]+1.0*gamma[13]*f[40]+1.0*gamma[11]*f[38]+1.0*gamma[7]*f[33]+gamma[10]*f[29]+gamma[6]*f[23]+gamma[5]*f[21]+gamma[4]*f[18]+gamma[3]*f[12]+gamma[2]*f[9]+gamma[1]*f[7]+gamma[0]*f[1])*volFact; 
  out[2] += (1.0*gamma[19]*f[78]+1.0*gamma[16]*f[73]+1.0*gamma[15]*f[71]+1.0*gamma[9]*f[66]+1.0*gamma[18]*f[62]+1.0*gamma[14]*f[57]+1.0*gamma[12]*f[55]+1.0*gamma[8]*f[50]+1.0*gamma[17]*f[46]+1.0*gamma[13]*f[41]+1.0*gamma[11]*f[39]+1.0*gamma[7]*f[34]+gamma[10]*f[30]+gamma[6]*f[24]+gamma[5]*f[22]+gamma[4]*f[19]+gamma[3]*f[13]+gamma[2]*f[10]+gamma[1]*f[8]+gamma[0]*f[2])*volFact; 
  out[3] += (gamma[19]*f[79]+gamma[16]*f[76]+gamma[15]*f[75]+gamma[9]*f[69]+gamma[18]*f[63]+gamma[14]*f[60]+gamma[12]*f[59]+gamma[8]*f[53]+gamma[17]*f[47]+gamma[13]*f[44]+gamma[11]*f[43]+gamma[7]*f[37]+gamma[10]*f[31]+gamma[6]*f[28]+gamma[5]*f[27]+gamma[4]*f[26]+gamma[3]*f[20]+gamma[2]*f[17]+gamma[1]*f[16]+gamma[0]*f[6])*volFact; 
  out[4] += volFact*(2.828427124746191*f[0]*wx1+0.8164965809277261*f[3]*dv1); 
  out[5] += volFact*(2.828427124746191*f[1]*wx1+0.8164965809277261*f[7]*dv1); 
  out[6] += volFact*(2.828427124746191*f[2]*wx1+0.8164965809277261*f[8]*dv1); 
  out[7] += volFact*(2.828427124746191*f[6]*wx1+0.8164965809277261*f[16]*dv1); 
  out[8] += volFact*(2.828427124746191*f[0]*wx2+0.8164965809277261*f[4]*dv2); 
  out[9] += volFact*(2.828427124746191*f[1]*wx2+0.8164965809277261*f[9]*dv2); 
  out[10] += volFact*(2.828427124746191*f[2]*wx2+0.8164965809277261*f[10]*dv2); 
  out[11] += volFact*(2.828427124746191*f[6]*wx2+0.8164965809277261*f[17]*dv2); 
  out[12] += volFact*(2.828427124746191*f[0]*wx3+0.8164965809277261*f[5]*dv3); 
  out[13] += volFact*(2.828427124746191*f[1]*wx3+0.8164965809277261*f[12]*dv3); 
  out[14] += volFact*(2.828427124746191*f[2]*wx3+0.8164965809277261*f[13]*dv3); 
  out[15] += volFact*(2.828427124746191*f[6]*wx3+0.8164965809277261*f[20]*dv3); 
  out[16] += volFact*(p0_over_gamma[16]*f[68]*wx1+p0_over_gamma[9]*f[64]*wx1+p0_over_gamma[14]*f[52]*wx1+p0_over_gamma[8]*f[48]*wx1+p0_over_gamma[10]*f[25]*wx1+p0_over_gamma[6]*f[15]*wx1+p0_over_gamma[5]*f[14]*wx1+p0_over_gamma[4]*f[11]*wx1+p0_over_gamma[3]*f[5]*wx1+p0_over_gamma[2]*f[4]*wx1+p0_over_gamma[1]*f[3]*wx1+f[0]*p0_over_gamma[0]*wx1+0.2886751345948129*p0_over_gamma[16]*f[74]*dv1+0.2886751345948129*p0_over_gamma[9]*f[67]*dv1+0.2886751345948129*p0_over_gamma[14]*f[58]*dv1+0.2886751345948129*p0_over_gamma[8]*f[51]*dv1+0.2581988897471612*p0_over_gamma[10]*f[42]*dv1+0.2581988897471611*p0_over_gamma[5]*f[36]*dv1+0.2581988897471611*p0_over_gamma[4]*f[35]*dv1+0.2581988897471612*p0_over_gamma[1]*f[32]*dv1+0.2886751345948129*p0_over_gamma[6]*f[25]*dv1+0.2886751345948129*p0_over_gamma[10]*f[15]*dv1+0.2886751345948129*p0_over_gamma[3]*f[14]*dv1+0.2886751345948129*p0_over_gamma[2]*f[11]*dv1+0.2886751345948129*f[5]*p0_over_gamma[5]*dv1+0.2886751345948129*f[4]*p0_over_gamma[4]*dv1+0.2886751345948129*p0_over_gamma[0]*f[3]*dv1+0.2886751345948129*f[0]*p0_over_gamma[1]*dv1); 
  out[17] += volFact*(1.0*p0_over_gamma[16]*f[72]*wx1+1.0*p0_over_gamma[9]*f[65]*wx1+1.0*p0_over_gamma[14]*f[56]*wx1+1.0*p0_over_gamma[8]*f[49]*wx1+p0_over_gamma[10]*f[29]*wx1+p0_over_gamma[6]*f[23]*wx1+p0_over_gamma[5]*f[21]*wx1+p0_over_gamma[4]*f[18]*wx1+p0_over_gamma[3]*f[12]*wx1+p0_over_gamma[2]*f[9]*wx1+p0_over_gamma[1]*f[7]*wx1+p0_over_gamma[0]*f[1]*wx1+0.2886751345948129*p0_over_gamma[16]*f[77]*dv1+0.2886751345948129*p0_over_gamma[9]*f[70]*dv1+0.2886751345948129*p0_over_gamma[14]*f[61]*dv1+0.2886751345948129*p0_over_gamma[8]*f[54]*dv1+0.2581988897471611*p0_over_gamma[10]*f[45]*dv1+0.2581988897471612*p0_over_gamma[5]*f[40]*dv1+0.2581988897471612*p0_over_gamma[4]*f[38]*dv1+0.2581988897471611*p0_over_gamma[1]*f[33]*dv1+0.2886751345948129*p0_over_gamma[6]*f[29]*dv1+0.2886751345948129*p0_over_gamma[10]*f[23]*dv1+0.2886751345948129*p0_over_gamma[3]*f[21]*dv1+0.2886751345948129*p0_over_gamma[2]*f[18]*dv1+0.2886751345948129*p0_over_gamma[5]*f[12]*dv1+0.2886751345948129*p0_over_gamma[4]*f[9]*dv1+0.2886751345948129*p0_over_gamma[0]*f[7]*dv1+0.2886751345948129*f[1]*p0_over_gamma[1]*dv1); 
  out[18] += volFact*(1.0*p0_over_gamma[16]*f[73]*wx1+1.0*p0_over_gamma[9]*f[66]*wx1+1.0*p0_over_gamma[14]*f[57]*wx1+1.0*p0_over_gamma[8]*f[50]*wx1+p0_over_gamma[10]*f[30]*wx1+p0_over_gamma[6]*f[24]*wx1+p0_over_gamma[5]*f[22]*wx1+p0_over_gamma[4]*f[19]*wx1+p0_over_gamma[3]*f[13]*wx1+p0_over_gamma[2]*f[10]*wx1+p0_over_gamma[1]*f[8]*wx1+p0_over_gamma[0]*f[2]*wx1+0.2886751345948129*p0_over_gamma[16]*f[78]*dv1+0.2886751345948129*p0_over_gamma[9]*f[71]*dv1+0.2886751345948129*p0_over_gamma[14]*f[62]*dv1+0.2886751345948129*p0_over_gamma[8]*f[55]*dv1+0.2581988897471611*p0_over_gamma[10]*f[46]*dv1+0.2581988897471612*p0_over_gamma[5]*f[41]*dv1+0.2581988897471612*p0_over_gamma[4]*f[39]*dv1+0.2581988897471611*p0_over_gamma[1]*f[34]*dv1+0.2886751345948129*p0_over_gamma[6]*f[30]*dv1+0.2886751345948129*p0_over_gamma[10]*f[24]*dv1+0.2886751345948129*p0_over_gamma[3]*f[22]*dv1+0.2886751345948129*p0_over_gamma[2]*f[19]*dv1+0.2886751345948129*p0_over_gamma[5]*f[13]*dv1+0.2886751345948129*p0_over_gamma[4]*f[10]*dv1+0.2886751345948129*p0_over_gamma[0]*f[8]*dv1+0.2886751345948129*p0_over_gamma[1]*f[2]*dv1); 
  out[19] += volFact*(p0_over_gamma[16]*f[76]*wx1+p0_over_gamma[9]*f[69]*wx1+p0_over_gamma[14]*f[60]*wx1+p0_over_gamma[8]*f[53]*wx1+p0_over_gamma[10]*f[31]*wx1+p0_over_gamma[6]*f[28]*wx1+p0_over_gamma[5]*f[27]*wx1+p0_over_gamma[4]*f[26]*wx1+p0_over_gamma[3]*f[20]*wx1+p0_over_gamma[2]*f[17]*wx1+p0_over_gamma[1]*f[16]*wx1+p0_over_gamma[0]*f[6]*wx1+0.2886751345948129*p0_over_gamma[16]*f[79]*dv1+0.2886751345948129*p0_over_gamma[9]*f[75]*dv1+0.2886751345948129*p0_over_gamma[14]*f[63]*dv1+0.2886751345948129*p0_over_gamma[8]*f[59]*dv1+0.2581988897471612*p0_over_gamma[10]*f[47]*dv1+0.2581988897471611*p0_over_gamma[5]*f[44]*dv1+0.2581988897471611*p0_over_gamma[4]*f[43]*dv1+0.2581988897471612*p0_over_gamma[1]*f[37]*dv1+0.2886751345948129*p0_over_gamma[6]*f[31]*dv1+0.2886751345948129*p0_over_gamma[10]*f[28]*dv1+0.2886751345948129*p0_over_gamma[3]*f[27]*dv1+0.2886751345948129*p0_over_gamma[2]*f[26]*dv1+0.2886751345948129*p0_over_gamma[5]*f[20]*dv1+0.2886751345948129*p0_over_gamma[4]*f[17]*dv1+0.2886751345948129*p0_over_gamma[0]*f[16]*dv1+0.2886751345948129*p0_over_gamma[1]*f[6]*dv1); 
  out[20] += volFact*(p0_over_gamma[16]*f[68]*wx2+p0_over_gamma[9]*f[64]*wx2+p0_over_gamma[14]*f[52]*wx2+p0_over_gamma[8]*f[48]*wx2+p0_over_gamma[10]*f[25]*wx2+p0_over_gamma[6]*f[15]*wx2+p0_over_gamma[5]*f[14]*wx2+p0_over_gamma[4]*f[11]*wx2+p0_over_gamma[3]*f[5]*wx2+p0_over_gamma[2]*f[4]*wx2+p0_over_gamma[1]*f[3]*wx2+f[0]*p0_over_gamma[0]*wx2+0.2886751345948129*p0_over_gamma[9]*f[68]*dv2+0.2886751345948129*p0_over_gamma[16]*f[64]*dv2+0.2581988897471612*p0_over_gamma[10]*f[58]*dv2+0.2581988897471611*p0_over_gamma[6]*f[52]*dv2+0.2581988897471611*p0_over_gamma[4]*f[51]*dv2+0.2581988897471612*p0_over_gamma[2]*f[48]*dv2+0.2886751345948129*p0_over_gamma[5]*f[25]*dv2+0.2581988897471611*p0_over_gamma[14]*f[15]*dv2+0.2886751345948129*p0_over_gamma[3]*f[15]*dv2+0.2886751345948129*p0_over_gamma[10]*f[14]*dv2+0.2886751345948129*p0_over_gamma[1]*f[11]*dv2+0.2581988897471612*f[4]*p0_over_gamma[8]*dv2+0.2886751345948129*f[5]*p0_over_gamma[6]*dv2+0.2886751345948129*f[3]*p0_over_gamma[4]*dv2+0.2886751345948129*p0_over_gamma[0]*f[4]*dv2+0.2886751345948129*f[0]*p0_over_gamma[2]*dv2); 
  out[21] += volFact*(1.0*p0_over_gamma[16]*f[72]*wx2+1.0*p0_over_gamma[9]*f[65]*wx2+1.0*p0_over_gamma[14]*f[56]*wx2+1.0*p0_over_gamma[8]*f[49]*wx2+p0_over_gamma[10]*f[29]*wx2+p0_over_gamma[6]*f[23]*wx2+p0_over_gamma[5]*f[21]*wx2+p0_over_gamma[4]*f[18]*wx2+p0_over_gamma[3]*f[12]*wx2+p0_over_gamma[2]*f[9]*wx2+p0_over_gamma[1]*f[7]*wx2+p0_over_gamma[0]*f[1]*wx2+0.2886751345948129*p0_over_gamma[9]*f[72]*dv2+0.2886751345948129*p0_over_gamma[16]*f[65]*dv2+0.2581988897471611*p0_over_gamma[10]*f[61]*dv2+0.2581988897471612*p0_over_gamma[6]*f[56]*dv2+0.2581988897471612*p0_over_gamma[4]*f[54]*dv2+0.2581988897471611*p0_over_gamma[2]*f[49]*dv2+0.2886751345948129*p0_over_gamma[5]*f[29]*dv2+0.2581988897471611*p0_over_gamma[14]*f[23]*dv2+0.2886751345948129*p0_over_gamma[3]*f[23]*dv2+0.2886751345948129*p0_over_gamma[10]*f[21]*dv2+0.2886751345948129*p0_over_gamma[1]*f[18]*dv2+0.2886751345948129*p0_over_gamma[6]*f[12]*dv2+0.2581988897471612*p0_over_gamma[8]*f[9]*dv2+0.2886751345948129*p0_over_gamma[0]*f[9]*dv2+0.2886751345948129*p0_over_gamma[4]*f[7]*dv2+0.2886751345948129*f[1]*p0_over_gamma[2]*dv2); 
  out[22] += volFact*(1.0*p0_over_gamma[16]*f[73]*wx2+1.0*p0_over_gamma[9]*f[66]*wx2+1.0*p0_over_gamma[14]*f[57]*wx2+1.0*p0_over_gamma[8]*f[50]*wx2+p0_over_gamma[10]*f[30]*wx2+p0_over_gamma[6]*f[24]*wx2+p0_over_gamma[5]*f[22]*wx2+p0_over_gamma[4]*f[19]*wx2+p0_over_gamma[3]*f[13]*wx2+p0_over_gamma[2]*f[10]*wx2+p0_over_gamma[1]*f[8]*wx2+p0_over_gamma[0]*f[2]*wx2+0.2886751345948129*p0_over_gamma[9]*f[73]*dv2+0.2886751345948129*p0_over_gamma[16]*f[66]*dv2+0.2581988897471611*p0_over_gamma[10]*f[62]*dv2+0.2581988897471612*p0_over_gamma[6]*f[57]*dv2+0.2581988897471612*p0_over_gamma[4]*f[55]*dv2+0.2581988897471611*p0_over_gamma[2]*f[50]*dv2+0.2886751345948129*p0_over_gamma[5]*f[30]*dv2+0.2581988897471611*p0_over_gamma[14]*f[24]*dv2+0.2886751345948129*p0_over_gamma[3]*f[24]*dv2+0.2886751345948129*p0_over_gamma[10]*f[22]*dv2+0.2886751345948129*p0_over_gamma[1]*f[19]*dv2+0.2886751345948129*p0_over_gamma[6]*f[13]*dv2+0.2581988897471612*p0_over_gamma[8]*f[10]*dv2+0.2886751345948129*p0_over_gamma[0]*f[10]*dv2+0.2886751345948129*p0_over_gamma[4]*f[8]*dv2+0.2886751345948129*f[2]*p0_over_gamma[2]*dv2); 
  out[23] += volFact*(p0_over_gamma[16]*f[76]*wx2+p0_over_gamma[9]*f[69]*wx2+p0_over_gamma[14]*f[60]*wx2+p0_over_gamma[8]*f[53]*wx2+p0_over_gamma[10]*f[31]*wx2+p0_over_gamma[6]*f[28]*wx2+p0_over_gamma[5]*f[27]*wx2+p0_over_gamma[4]*f[26]*wx2+p0_over_gamma[3]*f[20]*wx2+p0_over_gamma[2]*f[17]*wx2+p0_over_gamma[1]*f[16]*wx2+p0_over_gamma[0]*f[6]*wx2+0.2886751345948129*p0_over_gamma[9]*f[76]*dv2+0.2886751345948129*p0_over_gamma[16]*f[69]*dv2+0.2581988897471612*p0_over_gamma[10]*f[63]*dv2+0.2581988897471611*p0_over_gamma[6]*f[60]*dv2+0.2581988897471611*p0_over_gamma[4]*f[59]*dv2+0.2581988897471612*p0_over_gamma[2]*f[53]*dv2+0.2886751345948129*p0_over_gamma[5]*f[31]*dv2+0.2581988897471611*p0_over_gamma[14]*f[28]*dv2+0.2886751345948129*p0_over_gamma[3]*f[28]*dv2+0.2886751345948129*p0_over_gamma[10]*f[27]*dv2+0.2886751345948129*p0_over_gamma[1]*f[26]*dv2+0.2886751345948129*p0_over_gamma[6]*f[20]*dv2+0.2581988897471612*p0_over_gamma[8]*f[17]*dv2+0.2886751345948129*p0_over_gamma[0]*f[17]*dv2+0.2886751345948129*p0_over_gamma[4]*f[16]*dv2+0.2886751345948129*p0_over_gamma[2]*f[6]*dv2); 
  out[24] += volFact*(p0_over_gamma[16]*f[68]*wx3+p0_over_gamma[9]*f[64]*wx3+p0_over_gamma[14]*f[52]*wx3+p0_over_gamma[8]*f[48]*wx3+p0_over_gamma[10]*f[25]*wx3+p0_over_gamma[6]*f[15]*wx3+p0_over_gamma[5]*f[14]*wx3+p0_over_gamma[4]*f[11]*wx3+p0_over_gamma[3]*f[5]*wx3+p0_over_gamma[2]*f[4]*wx3+p0_over_gamma[1]*f[3]*wx3+f[0]*p0_over_gamma[0]*wx3+0.2581988897471612*p0_over_gamma[10]*f[74]*dv3+0.2581988897471611*p0_over_gamma[6]*f[68]*dv3+0.2581988897471611*p0_over_gamma[5]*f[67]*dv3+0.2581988897471612*p0_over_gamma[3]*f[64]*dv3+0.2886751345948129*p0_over_gamma[8]*f[52]*dv3+0.2886751345948129*p0_over_gamma[14]*f[48]*dv3+0.2886751345948129*p0_over_gamma[4]*f[25]*dv3+0.2581988897471611*f[15]*p0_over_gamma[16]*dv3+0.2886751345948129*p0_over_gamma[2]*f[15]*dv3+0.2886751345948129*p0_over_gamma[1]*f[14]*dv3+0.2886751345948129*p0_over_gamma[10]*f[11]*dv3+0.2581988897471612*f[5]*p0_over_gamma[9]*dv3+0.2886751345948129*f[4]*p0_over_gamma[6]*dv3+0.2886751345948129*f[3]*p0_over_gamma[5]*dv3+0.2886751345948129*p0_over_gamma[0]*f[5]*dv3+0.2886751345948129*f[0]*p0_over_gamma[3]*dv3); 
  out[25] += volFact*(1.0*p0_over_gamma[16]*f[72]*wx3+1.0*p0_over_gamma[9]*f[65]*wx3+1.0*p0_over_gamma[14]*f[56]*wx3+1.0*p0_over_gamma[8]*f[49]*wx3+p0_over_gamma[10]*f[29]*wx3+p0_over_gamma[6]*f[23]*wx3+p0_over_gamma[5]*f[21]*wx3+p0_over_gamma[4]*f[18]*wx3+p0_over_gamma[3]*f[12]*wx3+p0_over_gamma[2]*f[9]*wx3+p0_over_gamma[1]*f[7]*wx3+p0_over_gamma[0]*f[1]*wx3+0.2581988897471611*p0_over_gamma[10]*f[77]*dv3+0.2581988897471612*p0_over_gamma[6]*f[72]*dv3+0.2581988897471612*p0_over_gamma[5]*f[70]*dv3+0.2581988897471611*p0_over_gamma[3]*f[65]*dv3+0.2886751345948129*p0_over_gamma[8]*f[56]*dv3+0.2886751345948129*p0_over_gamma[14]*f[49]*dv3+0.2886751345948129*p0_over_gamma[4]*f[29]*dv3+0.2581988897471611*p0_over_gamma[16]*f[23]*dv3+0.2886751345948129*p0_over_gamma[2]*f[23]*dv3+0.2886751345948129*p0_over_gamma[1]*f[21]*dv3+0.2886751345948129*p0_over_gamma[10]*f[18]*dv3+0.2581988897471612*p0_over_gamma[9]*f[12]*dv3+0.2886751345948129*p0_over_gamma[0]*f[12]*dv3+0.2886751345948129*p0_over_gamma[6]*f[9]*dv3+0.2886751345948129*p0_over_gamma[5]*f[7]*dv3+0.2886751345948129*f[1]*p0_over_gamma[3]*dv3); 
  out[26] += volFact*(1.0*p0_over_gamma[16]*f[73]*wx3+1.0*p0_over_gamma[9]*f[66]*wx3+1.0*p0_over_gamma[14]*f[57]*wx3+1.0*p0_over_gamma[8]*f[50]*wx3+p0_over_gamma[10]*f[30]*wx3+p0_over_gamma[6]*f[24]*wx3+p0_over_gamma[5]*f[22]*wx3+p0_over_gamma[4]*f[19]*wx3+p0_over_gamma[3]*f[13]*wx3+p0_over_gamma[2]*f[10]*wx3+p0_over_gamma[1]*f[8]*wx3+p0_over_gamma[0]*f[2]*wx3+0.2581988897471611*p0_over_gamma[10]*f[78]*dv3+0.2581988897471612*p0_over_gamma[6]*f[73]*dv3+0.2581988897471612*p0_over_gamma[5]*f[71]*dv3+0.2581988897471611*p0_over_gamma[3]*f[66]*dv3+0.2886751345948129*p0_over_gamma[8]*f[57]*dv3+0.2886751345948129*p0_over_gamma[14]*f[50]*dv3+0.2886751345948129*p0_over_gamma[4]*f[30]*dv3+0.2581988897471611*p0_over_gamma[16]*f[24]*dv3+0.2886751345948129*p0_over_gamma[2]*f[24]*dv3+0.2886751345948129*p0_over_gamma[1]*f[22]*dv3+0.2886751345948129*p0_over_gamma[10]*f[19]*dv3+0.2581988897471612*p0_over_gamma[9]*f[13]*dv3+0.2886751345948129*p0_over_gamma[0]*f[13]*dv3+0.2886751345948129*p0_over_gamma[6]*f[10]*dv3+0.2886751345948129*p0_over_gamma[5]*f[8]*dv3+0.2886751345948129*f[2]*p0_over_gamma[3]*dv3); 
  out[27] += volFact*(p0_over_gamma[16]*f[76]*wx3+p0_over_gamma[9]*f[69]*wx3+p0_over_gamma[14]*f[60]*wx3+p0_over_gamma[8]*f[53]*wx3+p0_over_gamma[10]*f[31]*wx3+p0_over_gamma[6]*f[28]*wx3+p0_over_gamma[5]*f[27]*wx3+p0_over_gamma[4]*f[26]*wx3+p0_over_gamma[3]*f[20]*wx3+p0_over_gamma[2]*f[17]*wx3+p0_over_gamma[1]*f[16]*wx3+p0_over_gamma[0]*f[6]*wx3+0.2581988897471612*p0_over_gamma[10]*f[79]*dv3+0.2581988897471611*p0_over_gamma[6]*f[76]*dv3+0.2581988897471611*p0_over_gamma[5]*f[75]*dv3+0.2581988897471612*p0_over_gamma[3]*f[69]*dv3+0.2886751345948129*p0_over_gamma[8]*f[60]*dv3+0.2886751345948129*p0_over_gamma[14]*f[53]*dv3+0.2886751345948129*p0_over_gamma[4]*f[31]*dv3+0.2581988897471611*p0_over_gamma[16]*f[28]*dv3+0.2886751345948129*p0_over_gamma[2]*f[28]*dv3+0.2886751345948129*p0_over_gamma[1]*f[27]*dv3+0.2886751345948129*p0_over_gamma[10]*f[26]*dv3+0.2581988897471612*p0_over_gamma[9]*f[20]*dv3+0.2886751345948129*p0_over_gamma[0]*f[20]*dv3+0.2886751345948129*p0_over_gamma[6]*f[17]*dv3+0.2886751345948129*p0_over_gamma[5]*f[16]*dv3+0.2886751345948129*p0_over_gamma[3]*f[6]*dv3); 
  out[28] += volFact*(p1_over_gamma[15]*f[67]*wx2+p1_over_gamma[9]*f[64]*wx2+p1_over_gamma[13]*f[36]*wx2+p1_over_gamma[7]*f[32]*wx2+p1_over_gamma[10]*f[25]*wx2+p1_over_gamma[6]*f[15]*wx2+p1_over_gamma[5]*f[14]*wx2+p1_over_gamma[4]*f[11]*wx2+p1_over_gamma[3]*f[5]*wx2+p1_over_gamma[2]*f[4]*wx2+p1_over_gamma[1]*f[3]*wx2+f[0]*p1_over_gamma[0]*wx2+0.2886751345948129*p1_over_gamma[15]*f[74]*dv2+0.2886751345948129*p1_over_gamma[9]*f[68]*dv2+0.2581988897471612*p1_over_gamma[10]*f[58]*dv2+0.2581988897471611*p1_over_gamma[6]*f[52]*dv2+0.2581988897471611*p1_over_gamma[4]*f[51]*dv2+0.2581988897471612*p1_over_gamma[2]*f[48]*dv2+0.2886751345948129*p1_over_gamma[13]*f[42]*dv2+0.2886751345948129*p1_over_gamma[7]*f[35]*dv2+0.2886751345948129*p1_over_gamma[5]*f[25]*dv2+0.2886751345948129*p1_over_gamma[3]*f[15]*dv2+0.2886751345948129*p1_over_gamma[10]*f[14]*dv2+0.2886751345948129*p1_over_gamma[1]*f[11]*dv2+0.2886751345948129*f[5]*p1_over_gamma[6]*dv2+0.2886751345948129*f[3]*p1_over_gamma[4]*dv2+0.2886751345948129*p1_over_gamma[0]*f[4]*dv2+0.2886751345948129*f[0]*p1_over_gamma[2]*dv2); 
  out[29] += volFact*(1.0*p1_over_gamma[15]*f[70]*wx2+1.0*p1_over_gamma[9]*f[65]*wx2+1.0*p1_over_gamma[13]*f[40]*wx2+1.0*p1_over_gamma[7]*f[33]*wx2+p1_over_gamma[10]*f[29]*wx2+p1_over_gamma[6]*f[23]*wx2+p1_over_gamma[5]*f[21]*wx2+p1_over_gamma[4]*f[18]*wx2+p1_over_gamma[3]*f[12]*wx2+p1_over_gamma[2]*f[9]*wx2+p1_over_gamma[1]*f[7]*wx2+p1_over_gamma[0]*f[1]*wx2+0.2886751345948129*p1_over_gamma[15]*f[77]*dv2+0.2886751345948129*p1_over_gamma[9]*f[72]*dv2+0.2581988897471611*p1_over_gamma[10]*f[61]*dv2+0.2581988897471612*p1_over_gamma[6]*f[56]*dv2+0.2581988897471612*p1_over_gamma[4]*f[54]*dv2+0.2581988897471611*p1_over_gamma[2]*f[49]*dv2+0.2886751345948129*p1_over_gamma[13]*f[45]*dv2+0.2886751345948129*p1_over_gamma[7]*f[38]*dv2+0.2886751345948129*p1_over_gamma[5]*f[29]*dv2+0.2886751345948129*p1_over_gamma[3]*f[23]*dv2+0.2886751345948129*p1_over_gamma[10]*f[21]*dv2+0.2886751345948129*p1_over_gamma[1]*f[18]*dv2+0.2886751345948129*p1_over_gamma[6]*f[12]*dv2+0.2886751345948129*p1_over_gamma[0]*f[9]*dv2+0.2886751345948129*p1_over_gamma[4]*f[7]*dv2+0.2886751345948129*f[1]*p1_over_gamma[2]*dv2); 
  out[30] += volFact*(1.0*p1_over_gamma[15]*f[71]*wx2+1.0*p1_over_gamma[9]*f[66]*wx2+1.0*p1_over_gamma[13]*f[41]*wx2+1.0*p1_over_gamma[7]*f[34]*wx2+p1_over_gamma[10]*f[30]*wx2+p1_over_gamma[6]*f[24]*wx2+p1_over_gamma[5]*f[22]*wx2+p1_over_gamma[4]*f[19]*wx2+p1_over_gamma[3]*f[13]*wx2+p1_over_gamma[2]*f[10]*wx2+p1_over_gamma[1]*f[8]*wx2+p1_over_gamma[0]*f[2]*wx2+0.2886751345948129*p1_over_gamma[15]*f[78]*dv2+0.2886751345948129*p1_over_gamma[9]*f[73]*dv2+0.2581988897471611*p1_over_gamma[10]*f[62]*dv2+0.2581988897471612*p1_over_gamma[6]*f[57]*dv2+0.2581988897471612*p1_over_gamma[4]*f[55]*dv2+0.2581988897471611*p1_over_gamma[2]*f[50]*dv2+0.2886751345948129*p1_over_gamma[13]*f[46]*dv2+0.2886751345948129*p1_over_gamma[7]*f[39]*dv2+0.2886751345948129*p1_over_gamma[5]*f[30]*dv2+0.2886751345948129*p1_over_gamma[3]*f[24]*dv2+0.2886751345948129*p1_over_gamma[10]*f[22]*dv2+0.2886751345948129*p1_over_gamma[1]*f[19]*dv2+0.2886751345948129*p1_over_gamma[6]*f[13]*dv2+0.2886751345948129*p1_over_gamma[0]*f[10]*dv2+0.2886751345948129*p1_over_gamma[4]*f[8]*dv2+0.2886751345948129*f[2]*p1_over_gamma[2]*dv2); 
  out[31] += volFact*(p1_over_gamma[15]*f[75]*wx2+p1_over_gamma[9]*f[69]*wx2+p1_over_gamma[13]*f[44]*wx2+p1_over_gamma[7]*f[37]*wx2+p1_over_gamma[10]*f[31]*wx2+p1_over_gamma[6]*f[28]*wx2+p1_over_gamma[5]*f[27]*wx2+p1_over_gamma[4]*f[26]*wx2+p1_over_gamma[3]*f[20]*wx2+p1_over_gamma[2]*f[17]*wx2+p1_over_gamma[1]*f[16]*wx2+p1_over_gamma[0]*f[6]*wx2+0.2886751345948129*p1_over_gamma[15]*f[79]*dv2+0.2886751345948129*p1_over_gamma[9]*f[76]*dv2+0.2581988897471612*p1_over_gamma[10]*f[63]*dv2+0.2581988897471611*p1_over_gamma[6]*f[60]*dv2+0.2581988897471611*p1_over_gamma[4]*f[59]*dv2+0.2581988897471612*p1_over_gamma[2]*f[53]*dv2+0.2886751345948129*p1_over_gamma[13]*f[47]*dv2+0.2886751345948129*p1_over_gamma[7]*f[43]*dv2+0.2886751345948129*p1_over_gamma[5]*f[31]*dv2+0.2886751345948129*p1_over_gamma[3]*f[28]*dv2+0.2886751345948129*p1_over_gamma[10]*f[27]*dv2+0.2886751345948129*p1_over_gamma[1]*f[26]*dv2+0.2886751345948129*p1_over_gamma[6]*f[20]*dv2+0.2886751345948129*p1_over_gamma[0]*f[17]*dv2+0.2886751345948129*p1_over_gamma[4]*f[16]*dv2+0.2886751345948129*p1_over_gamma[2]*f[6]*dv2); 
  out[32] += volFact*(p1_over_gamma[15]*f[67]*wx3+p1_over_gamma[9]*f[64]*wx3+p1_over_gamma[13]*f[36]*wx3+p1_over_gamma[7]*f[32]*wx3+p1_over_gamma[10]*f[25]*wx3+p1_over_gamma[6]*f[15]*wx3+p1_over_gamma[5]*f[14]*wx3+p1_over_gamma[4]*f[11]*wx3+p1_over_gamma[3]*f[5]*wx3+p1_over_gamma[2]*f[4]*wx3+p1_over_gamma[1]*f[3]*wx3+f[0]*p1_over_gamma[0]*wx3+0.2581988897471612*p1_over_gamma[10]*f[74]*dv3+0.2581988897471611*p1_over_gamma[6]*f[68]*dv3+0.2581988897471611*p1_over_gamma[5]*f[67]*dv3+0.2581988897471612*p1_over_gamma[3]*f[64]*dv3+0.2886751345948129*p1_over_gamma[7]*f[36]*dv3+0.2886751345948129*p1_over_gamma[13]*f[32]*dv3+0.2886751345948129*p1_over_gamma[4]*f[25]*dv3+0.2581988897471611*f[14]*p1_over_gamma[15]*dv3+0.2886751345948129*p1_over_gamma[2]*f[15]*dv3+0.2886751345948129*p1_over_gamma[1]*f[14]*dv3+0.2886751345948129*p1_over_gamma[10]*f[11]*dv3+0.2581988897471612*f[5]*p1_over_gamma[9]*dv3+0.2886751345948129*f[4]*p1_over_gamma[6]*dv3+0.2886751345948129*f[3]*p1_over_gamma[5]*dv3+0.2886751345948129*p1_over_gamma[0]*f[5]*dv3+0.2886751345948129*f[0]*p1_over_gamma[3]*dv3); 
  out[33] += volFact*(1.0*p1_over_gamma[15]*f[70]*wx3+1.0*p1_over_gamma[9]*f[65]*wx3+1.0*p1_over_gamma[13]*f[40]*wx3+1.0*p1_over_gamma[7]*f[33]*wx3+p1_over_gamma[10]*f[29]*wx3+p1_over_gamma[6]*f[23]*wx3+p1_over_gamma[5]*f[21]*wx3+p1_over_gamma[4]*f[18]*wx3+p1_over_gamma[3]*f[12]*wx3+p1_over_gamma[2]*f[9]*wx3+p1_over_gamma[1]*f[7]*wx3+p1_over_gamma[0]*f[1]*wx3+0.2581988897471611*p1_over_gamma[10]*f[77]*dv3+0.2581988897471612*p1_over_gamma[6]*f[72]*dv3+0.2581988897471612*p1_over_gamma[5]*f[70]*dv3+0.2581988897471611*p1_over_gamma[3]*f[65]*dv3+0.2886751345948129*p1_over_gamma[7]*f[40]*dv3+0.2886751345948129*p1_over_gamma[13]*f[33]*dv3+0.2886751345948129*p1_over_gamma[4]*f[29]*dv3+0.2886751345948129*p1_over_gamma[2]*f[23]*dv3+0.2581988897471611*p1_over_gamma[15]*f[21]*dv3+0.2886751345948129*p1_over_gamma[1]*f[21]*dv3+0.2886751345948129*p1_over_gamma[10]*f[18]*dv3+0.2581988897471612*p1_over_gamma[9]*f[12]*dv3+0.2886751345948129*p1_over_gamma[0]*f[12]*dv3+0.2886751345948129*p1_over_gamma[6]*f[9]*dv3+0.2886751345948129*p1_over_gamma[5]*f[7]*dv3+0.2886751345948129*f[1]*p1_over_gamma[3]*dv3); 
  out[34] += volFact*(1.0*p1_over_gamma[15]*f[71]*wx3+1.0*p1_over_gamma[9]*f[66]*wx3+1.0*p1_over_gamma[13]*f[41]*wx3+1.0*p1_over_gamma[7]*f[34]*wx3+p1_over_gamma[10]*f[30]*wx3+p1_over_gamma[6]*f[24]*wx3+p1_over_gamma[5]*f[22]*wx3+p1_over_gamma[4]*f[19]*wx3+p1_over_gamma[3]*f[13]*wx3+p1_over_gamma[2]*f[10]*wx3+p1_over_gamma[1]*f[8]*wx3+p1_over_gamma[0]*f[2]*wx3+0.2581988897471611*p1_over_gamma[10]*f[78]*dv3+0.2581988897471612*p1_over_gamma[6]*f[73]*dv3+0.2581988897471612*p1_over_gamma[5]*f[71]*dv3+0.2581988897471611*p1_over_gamma[3]*f[66]*dv3+0.2886751345948129*p1_over_gamma[7]*f[41]*dv3+0.2886751345948129*p1_over_gamma[13]*f[34]*dv3+0.2886751345948129*p1_over_gamma[4]*f[30]*dv3+0.2886751345948129*p1_over_gamma[2]*f[24]*dv3+0.2581988897471611*p1_over_gamma[15]*f[22]*dv3+0.2886751345948129*p1_over_gamma[1]*f[22]*dv3+0.2886751345948129*p1_over_gamma[10]*f[19]*dv3+0.2581988897471612*p1_over_gamma[9]*f[13]*dv3+0.2886751345948129*p1_over_gamma[0]*f[13]*dv3+0.2886751345948129*p1_over_gamma[6]*f[10]*dv3+0.2886751345948129*p1_over_gamma[5]*f[8]*dv3+0.2886751345948129*f[2]*p1_over_gamma[3]*dv3); 
  out[35] += volFact*(p1_over_gamma[15]*f[75]*wx3+p1_over_gamma[9]*f[69]*wx3+p1_over_gamma[13]*f[44]*wx3+p1_over_gamma[7]*f[37]*wx3+p1_over_gamma[10]*f[31]*wx3+p1_over_gamma[6]*f[28]*wx3+p1_over_gamma[5]*f[27]*wx3+p1_over_gamma[4]*f[26]*wx3+p1_over_gamma[3]*f[20]*wx3+p1_over_gamma[2]*f[17]*wx3+p1_over_gamma[1]*f[16]*wx3+p1_over_gamma[0]*f[6]*wx3+0.2581988897471612*p1_over_gamma[10]*f[79]*dv3+0.2581988897471611*p1_over_gamma[6]*f[76]*dv3+0.2581988897471611*p1_over_gamma[5]*f[75]*dv3+0.2581988897471612*p1_over_gamma[3]*f[69]*dv3+0.2886751345948129*p1_over_gamma[7]*f[44]*dv3+0.2886751345948129*p1_over_gamma[13]*f[37]*dv3+0.2886751345948129*p1_over_gamma[4]*f[31]*dv3+0.2886751345948129*p1_over_gamma[2]*f[28]*dv3+0.2581988897471611*p1_over_gamma[15]*f[27]*dv3+0.2886751345948129*p1_over_gamma[1]*f[27]*dv3+0.2886751345948129*p1_over_gamma[10]*f[26]*dv3+0.2581988897471612*p1_over_gamma[9]*f[20]*dv3+0.2886751345948129*p1_over_gamma[0]*f[20]*dv3+0.2886751345948129*p1_over_gamma[6]*f[17]*dv3+0.2886751345948129*p1_over_gamma[5]*f[16]*dv3+0.2886751345948129*p1_over_gamma[3]*f[6]*dv3); 
  out[36] += volFact*(p2_over_gamma[12]*f[51]*wx3+p2_over_gamma[8]*f[48]*wx3+p2_over_gamma[11]*f[35]*wx3+p2_over_gamma[7]*f[32]*wx3+p2_over_gamma[10]*f[25]*wx3+p2_over_gamma[6]*f[15]*wx3+p2_over_gamma[5]*f[14]*wx3+p2_over_gamma[4]*f[11]*wx3+p2_over_gamma[3]*f[5]*wx3+p2_over_gamma[2]*f[4]*wx3+p2_over_gamma[1]*f[3]*wx3+f[0]*p2_over_gamma[0]*wx3+0.2581988897471612*p2_over_gamma[10]*f[74]*dv3+0.2581988897471611*p2_over_gamma[6]*f[68]*dv3+0.2581988897471611*p2_over_gamma[5]*f[67]*dv3+0.2581988897471612*p2_over_gamma[3]*f[64]*dv3+0.2886751345948129*p2_over_gamma[12]*f[58]*dv3+0.2886751345948129*p2_over_gamma[8]*f[52]*dv3+0.2886751345948129*p2_over_gamma[11]*f[42]*dv3+0.2886751345948129*p2_over_gamma[7]*f[36]*dv3+0.2886751345948129*p2_over_gamma[4]*f[25]*dv3+0.2886751345948129*p2_over_gamma[2]*f[15]*dv3+0.2886751345948129*p2_over_gamma[1]*f[14]*dv3+0.2886751345948129*p2_over_gamma[10]*f[11]*dv3+0.2886751345948129*f[4]*p2_over_gamma[6]*dv3+0.2886751345948129*f[3]*p2_over_gamma[5]*dv3+0.2886751345948129*p2_over_gamma[0]*f[5]*dv3+0.2886751345948129*f[0]*p2_over_gamma[3]*dv3); 
  out[37] += volFact*(1.0*p2_over_gamma[12]*f[54]*wx3+1.0*p2_over_gamma[8]*f[49]*wx3+1.0*p2_over_gamma[11]*f[38]*wx3+1.0*p2_over_gamma[7]*f[33]*wx3+p2_over_gamma[10]*f[29]*wx3+p2_over_gamma[6]*f[23]*wx3+p2_over_gamma[5]*f[21]*wx3+p2_over_gamma[4]*f[18]*wx3+p2_over_gamma[3]*f[12]*wx3+p2_over_gamma[2]*f[9]*wx3+p2_over_gamma[1]*f[7]*wx3+p2_over_gamma[0]*f[1]*wx3+0.2581988897471611*p2_over_gamma[10]*f[77]*dv3+0.2581988897471612*p2_over_gamma[6]*f[72]*dv3+0.2581988897471612*p2_over_gamma[5]*f[70]*dv3+0.2581988897471611*p2_over_gamma[3]*f[65]*dv3+0.2886751345948129*p2_over_gamma[12]*f[61]*dv3+0.2886751345948129*p2_over_gamma[8]*f[56]*dv3+0.2886751345948129*p2_over_gamma[11]*f[45]*dv3+0.2886751345948129*p2_over_gamma[7]*f[40]*dv3+0.2886751345948129*p2_over_gamma[4]*f[29]*dv3+0.2886751345948129*p2_over_gamma[2]*f[23]*dv3+0.2886751345948129*p2_over_gamma[1]*f[21]*dv3+0.2886751345948129*p2_over_gamma[10]*f[18]*dv3+0.2886751345948129*p2_over_gamma[0]*f[12]*dv3+0.2886751345948129*p2_over_gamma[6]*f[9]*dv3+0.2886751345948129*p2_over_gamma[5]*f[7]*dv3+0.2886751345948129*f[1]*p2_over_gamma[3]*dv3); 
  out[38] += volFact*(1.0*p2_over_gamma[12]*f[55]*wx3+1.0*p2_over_gamma[8]*f[50]*wx3+1.0*p2_over_gamma[11]*f[39]*wx3+1.0*p2_over_gamma[7]*f[34]*wx3+p2_over_gamma[10]*f[30]*wx3+p2_over_gamma[6]*f[24]*wx3+p2_over_gamma[5]*f[22]*wx3+p2_over_gamma[4]*f[19]*wx3+p2_over_gamma[3]*f[13]*wx3+p2_over_gamma[2]*f[10]*wx3+p2_over_gamma[1]*f[8]*wx3+p2_over_gamma[0]*f[2]*wx3+0.2581988897471611*p2_over_gamma[10]*f[78]*dv3+0.2581988897471612*p2_over_gamma[6]*f[73]*dv3+0.2581988897471612*p2_over_gamma[5]*f[71]*dv3+0.2581988897471611*p2_over_gamma[3]*f[66]*dv3+0.2886751345948129*p2_over_gamma[12]*f[62]*dv3+0.2886751345948129*p2_over_gamma[8]*f[57]*dv3+0.2886751345948129*p2_over_gamma[11]*f[46]*dv3+0.2886751345948129*p2_over_gamma[7]*f[41]*dv3+0.2886751345948129*p2_over_gamma[4]*f[30]*dv3+0.2886751345948129*p2_over_gamma[2]*f[24]*dv3+0.2886751345948129*p2_over_gamma[1]*f[22]*dv3+0.2886751345948129*p2_over_gamma[10]*f[19]*dv3+0.2886751345948129*p2_over_gamma[0]*f[13]*dv3+0.2886751345948129*p2_over_gamma[6]*f[10]*dv3+0.2886751345948129*p2_over_gamma[5]*f[8]*dv3+0.2886751345948129*f[2]*p2_over_gamma[3]*dv3); 
  out[39] += volFact*(p2_over_gamma[12]*f[59]*wx3+p2_over_gamma[8]*f[53]*wx3+p2_over_gamma[11]*f[43]*wx3+p2_over_gamma[7]*f[37]*wx3+p2_over_gamma[10]*f[31]*wx3+p2_over_gamma[6]*f[28]*wx3+p2_over_gamma[5]*f[27]*wx3+p2_over_gamma[4]*f[26]*wx3+p2_over_gamma[3]*f[20]*wx3+p2_over_gamma[2]*f[17]*wx3+p2_over_gamma[1]*f[16]*wx3+p2_over_gamma[0]*f[6]*wx3+0.2581988897471612*p2_over_gamma[10]*f[79]*dv3+0.2581988897471611*p2_over_gamma[6]*f[76]*dv3+0.2581988897471611*p2_over_gamma[5]*f[75]*dv3+0.2581988897471612*p2_over_gamma[3]*f[69]*dv3+0.2886751345948129*p2_over_gamma[12]*f[63]*dv3+0.2886751345948129*p2_over_gamma[8]*f[60]*dv3+0.2886751345948129*p2_over_gamma[11]*f[47]*dv3+0.2886751345948129*p2_over_gamma[7]*f[44]*dv3+0.2886751345948129*p2_over_gamma[4]*f[31]*dv3+0.2886751345948129*p2_over_gamma[2]*f[28]*dv3+0.2886751345948129*p2_over_gamma[1]*f[27]*dv3+0.2886751345948129*p2_over_gamma[10]*f[26]*dv3+0.2886751345948129*p2_over_gamma[0]*f[20]*dv3+0.2886751345948129*p2_over_gamma[6]*f[17]*dv3+0.2886751345948129*p2_over_gamma[5]*f[16]*dv3+0.2886751345948129*p2_over_gamma[3]*f[6]*dv3); 
} 
GKYL_CU_DH void vlasov_sr_int_mom_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*0.03125; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += (2.0*gamma[19]*f[74]+2.0*gamma[16]*f[68]+2.0*gamma[15]*f[67]+2.0*gamma[9]*f[64]+2.0*gamma[18]*f[58]+2.0*gamma[14]*f[52]+2.0*gamma[12]*f[51]+2.0*gamma[8]*f[48]+2.0*gamma[17]*f[42]+2.0*gamma[13]*f[36]+2.0*gamma[11]*f[35]+2.0*gamma[7]*f[32]+2.0*gamma[10]*f[25]+2.0*gamma[6]*f[15]+2.0*gamma[5]*f[14]+2.0*gamma[4]*f[11]+2.0*gamma[3]*f[5]+2.0*gamma[2]*f[4]+2.0*gamma[1]*f[3]+2.0*f[0]*gamma[0])*volFact; 
  out[2] += volFact*(5.656854249492382*f[0]*wx1+1.632993161855453*f[3]*dv1); 
  out[3] += volFact*(5.656854249492382*f[0]*wx2+1.632993161855453*f[4]*dv2); 
  out[4] += volFact*(5.656854249492382*f[0]*wx3+1.632993161855453*f[5]*dv3); 
} 
