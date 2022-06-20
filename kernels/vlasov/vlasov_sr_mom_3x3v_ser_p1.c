#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[3]*volFact; 
  out[4] += 2.828427124746191*f[7]*volFact; 
  out[5] += 2.828427124746191*f[8]*volFact; 
  out[6] += 2.828427124746191*f[9]*volFact; 
  out[7] += 2.828427124746191*f[22]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 
  out[0] += (p0_over_gamma[7]*f[41]+p0_over_gamma[6]*f[21]+p0_over_gamma[5]*f[20]+p0_over_gamma[4]*f[16]+p0_over_gamma[3]*f[6]+p0_over_gamma[2]*f[5]+p0_over_gamma[1]*f[4]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (p0_over_gamma[7]*f[54]+p0_over_gamma[6]*f[38]+p0_over_gamma[5]*f[35]+p0_over_gamma[4]*f[29]+p0_over_gamma[3]*f[17]+p0_over_gamma[2]*f[13]+p0_over_gamma[1]*f[10]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (p0_over_gamma[7]*f[55]+p0_over_gamma[6]*f[39]+p0_over_gamma[5]*f[36]+p0_over_gamma[4]*f[30]+p0_over_gamma[3]*f[18]+p0_over_gamma[2]*f[14]+p0_over_gamma[1]*f[11]+p0_over_gamma[0]*f[2])*volFact; 
  out[3] += (p0_over_gamma[7]*f[56]+p0_over_gamma[6]*f[40]+p0_over_gamma[5]*f[37]+p0_over_gamma[4]*f[31]+p0_over_gamma[3]*f[19]+p0_over_gamma[2]*f[15]+p0_over_gamma[1]*f[12]+p0_over_gamma[0]*f[3])*volFact; 
  out[4] += (p0_over_gamma[7]*f[60]+p0_over_gamma[6]*f[51]+p0_over_gamma[5]*f[48]+p0_over_gamma[4]*f[44]+p0_over_gamma[3]*f[32]+p0_over_gamma[2]*f[26]+p0_over_gamma[1]*f[23]+p0_over_gamma[0]*f[7])*volFact; 
  out[5] += (p0_over_gamma[7]*f[61]+p0_over_gamma[6]*f[52]+p0_over_gamma[5]*f[49]+p0_over_gamma[4]*f[45]+p0_over_gamma[3]*f[33]+p0_over_gamma[2]*f[27]+p0_over_gamma[1]*f[24]+p0_over_gamma[0]*f[8])*volFact; 
  out[6] += (p0_over_gamma[7]*f[62]+p0_over_gamma[6]*f[53]+p0_over_gamma[5]*f[50]+p0_over_gamma[4]*f[46]+p0_over_gamma[3]*f[34]+p0_over_gamma[2]*f[28]+p0_over_gamma[1]*f[25]+p0_over_gamma[0]*f[9])*volFact; 
  out[7] += (p0_over_gamma[7]*f[63]+p0_over_gamma[6]*f[59]+p0_over_gamma[5]*f[58]+p0_over_gamma[4]*f[57]+p0_over_gamma[3]*f[47]+p0_over_gamma[2]*f[43]+p0_over_gamma[1]*f[42]+p0_over_gamma[0]*f[22])*volFact; 
  out[8] += (p1_over_gamma[7]*f[41]+p1_over_gamma[6]*f[21]+p1_over_gamma[5]*f[20]+p1_over_gamma[4]*f[16]+p1_over_gamma[3]*f[6]+p1_over_gamma[2]*f[5]+p1_over_gamma[1]*f[4]+f[0]*p1_over_gamma[0])*volFact; 
  out[9] += (p1_over_gamma[7]*f[54]+p1_over_gamma[6]*f[38]+p1_over_gamma[5]*f[35]+p1_over_gamma[4]*f[29]+p1_over_gamma[3]*f[17]+p1_over_gamma[2]*f[13]+p1_over_gamma[1]*f[10]+p1_over_gamma[0]*f[1])*volFact; 
  out[10] += (p1_over_gamma[7]*f[55]+p1_over_gamma[6]*f[39]+p1_over_gamma[5]*f[36]+p1_over_gamma[4]*f[30]+p1_over_gamma[3]*f[18]+p1_over_gamma[2]*f[14]+p1_over_gamma[1]*f[11]+p1_over_gamma[0]*f[2])*volFact; 
  out[11] += (p1_over_gamma[7]*f[56]+p1_over_gamma[6]*f[40]+p1_over_gamma[5]*f[37]+p1_over_gamma[4]*f[31]+p1_over_gamma[3]*f[19]+p1_over_gamma[2]*f[15]+p1_over_gamma[1]*f[12]+p1_over_gamma[0]*f[3])*volFact; 
  out[12] += (p1_over_gamma[7]*f[60]+p1_over_gamma[6]*f[51]+p1_over_gamma[5]*f[48]+p1_over_gamma[4]*f[44]+p1_over_gamma[3]*f[32]+p1_over_gamma[2]*f[26]+p1_over_gamma[1]*f[23]+p1_over_gamma[0]*f[7])*volFact; 
  out[13] += (p1_over_gamma[7]*f[61]+p1_over_gamma[6]*f[52]+p1_over_gamma[5]*f[49]+p1_over_gamma[4]*f[45]+p1_over_gamma[3]*f[33]+p1_over_gamma[2]*f[27]+p1_over_gamma[1]*f[24]+p1_over_gamma[0]*f[8])*volFact; 
  out[14] += (p1_over_gamma[7]*f[62]+p1_over_gamma[6]*f[53]+p1_over_gamma[5]*f[50]+p1_over_gamma[4]*f[46]+p1_over_gamma[3]*f[34]+p1_over_gamma[2]*f[28]+p1_over_gamma[1]*f[25]+p1_over_gamma[0]*f[9])*volFact; 
  out[15] += (p1_over_gamma[7]*f[63]+p1_over_gamma[6]*f[59]+p1_over_gamma[5]*f[58]+p1_over_gamma[4]*f[57]+p1_over_gamma[3]*f[47]+p1_over_gamma[2]*f[43]+p1_over_gamma[1]*f[42]+p1_over_gamma[0]*f[22])*volFact; 
  out[16] += (p2_over_gamma[7]*f[41]+p2_over_gamma[6]*f[21]+p2_over_gamma[5]*f[20]+p2_over_gamma[4]*f[16]+p2_over_gamma[3]*f[6]+p2_over_gamma[2]*f[5]+p2_over_gamma[1]*f[4]+f[0]*p2_over_gamma[0])*volFact; 
  out[17] += (p2_over_gamma[7]*f[54]+p2_over_gamma[6]*f[38]+p2_over_gamma[5]*f[35]+p2_over_gamma[4]*f[29]+p2_over_gamma[3]*f[17]+p2_over_gamma[2]*f[13]+p2_over_gamma[1]*f[10]+p2_over_gamma[0]*f[1])*volFact; 
  out[18] += (p2_over_gamma[7]*f[55]+p2_over_gamma[6]*f[39]+p2_over_gamma[5]*f[36]+p2_over_gamma[4]*f[30]+p2_over_gamma[3]*f[18]+p2_over_gamma[2]*f[14]+p2_over_gamma[1]*f[11]+p2_over_gamma[0]*f[2])*volFact; 
  out[19] += (p2_over_gamma[7]*f[56]+p2_over_gamma[6]*f[40]+p2_over_gamma[5]*f[37]+p2_over_gamma[4]*f[31]+p2_over_gamma[3]*f[19]+p2_over_gamma[2]*f[15]+p2_over_gamma[1]*f[12]+p2_over_gamma[0]*f[3])*volFact; 
  out[20] += (p2_over_gamma[7]*f[60]+p2_over_gamma[6]*f[51]+p2_over_gamma[5]*f[48]+p2_over_gamma[4]*f[44]+p2_over_gamma[3]*f[32]+p2_over_gamma[2]*f[26]+p2_over_gamma[1]*f[23]+p2_over_gamma[0]*f[7])*volFact; 
  out[21] += (p2_over_gamma[7]*f[61]+p2_over_gamma[6]*f[52]+p2_over_gamma[5]*f[49]+p2_over_gamma[4]*f[45]+p2_over_gamma[3]*f[33]+p2_over_gamma[2]*f[27]+p2_over_gamma[1]*f[24]+p2_over_gamma[0]*f[8])*volFact; 
  out[22] += (p2_over_gamma[7]*f[62]+p2_over_gamma[6]*f[53]+p2_over_gamma[5]*f[50]+p2_over_gamma[4]*f[46]+p2_over_gamma[3]*f[34]+p2_over_gamma[2]*f[28]+p2_over_gamma[1]*f[25]+p2_over_gamma[0]*f[9])*volFact; 
  out[23] += (p2_over_gamma[7]*f[63]+p2_over_gamma[6]*f[59]+p2_over_gamma[5]*f[58]+p2_over_gamma[4]*f[57]+p2_over_gamma[3]*f[47]+p2_over_gamma[2]*f[43]+p2_over_gamma[1]*f[42]+p2_over_gamma[0]*f[22])*volFact; 
} 
