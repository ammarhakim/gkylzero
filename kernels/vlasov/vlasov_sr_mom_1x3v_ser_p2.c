#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[11]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_1x3v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *p1_over_gamma = &p_over_gamma[20]; 
  const double *p2_over_gamma = &p_over_gamma[40]; 
  out[0] += (p0_over_gamma[19]*f[43]+p0_over_gamma[18]*f[40]+p0_over_gamma[17]*f[38]+p0_over_gamma[16]*f[30]+p0_over_gamma[15]*f[29]+p0_over_gamma[14]*f[27]+p0_over_gamma[13]*f[26]+p0_over_gamma[12]*f[24]+p0_over_gamma[11]*f[22]+p0_over_gamma[10]*f[18]+p0_over_gamma[9]*f[14]+p0_over_gamma[8]*f[13]+p0_over_gamma[7]*f[12]+p0_over_gamma[6]*f[10]+p0_over_gamma[5]*f[9]+p0_over_gamma[4]*f[7]+p0_over_gamma[3]*f[4]+p0_over_gamma[2]*f[3]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (1.0*p0_over_gamma[19]*f[47]+1.0*p0_over_gamma[18]*f[46]+1.0*p0_over_gamma[17]*f[45]+1.0*p0_over_gamma[16]*f[42]+1.0*p0_over_gamma[15]*f[41]+1.0*p0_over_gamma[14]*f[39]+1.0*p0_over_gamma[13]*f[36]+1.0*p0_over_gamma[12]*f[34]+1.0*p0_over_gamma[11]*f[33]+p0_over_gamma[10]*f[31]+1.0*p0_over_gamma[9]*f[28]+1.0*p0_over_gamma[8]*f[23]+1.0*p0_over_gamma[7]*f[20]+p0_over_gamma[6]*f[17]+p0_over_gamma[5]*f[16]+p0_over_gamma[4]*f[15]+p0_over_gamma[3]*f[8]+p0_over_gamma[2]*f[6]+p0_over_gamma[1]*f[5]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (1.0*p0_over_gamma[10]*f[44]+p0_over_gamma[6]*f[37]+p0_over_gamma[5]*f[35]+p0_over_gamma[4]*f[32]+1.0*p0_over_gamma[3]*f[25]+1.0*p0_over_gamma[2]*f[21]+1.0*p0_over_gamma[1]*f[19]+p0_over_gamma[0]*f[11])*volFact; 
  out[3] += (p1_over_gamma[19]*f[43]+p1_over_gamma[18]*f[40]+p1_over_gamma[17]*f[38]+p1_over_gamma[16]*f[30]+p1_over_gamma[15]*f[29]+p1_over_gamma[14]*f[27]+p1_over_gamma[13]*f[26]+p1_over_gamma[12]*f[24]+p1_over_gamma[11]*f[22]+p1_over_gamma[10]*f[18]+p1_over_gamma[9]*f[14]+p1_over_gamma[8]*f[13]+p1_over_gamma[7]*f[12]+p1_over_gamma[6]*f[10]+p1_over_gamma[5]*f[9]+p1_over_gamma[4]*f[7]+p1_over_gamma[3]*f[4]+p1_over_gamma[2]*f[3]+p1_over_gamma[1]*f[2]+f[0]*p1_over_gamma[0])*volFact; 
  out[4] += (1.0*p1_over_gamma[19]*f[47]+1.0*p1_over_gamma[18]*f[46]+1.0*p1_over_gamma[17]*f[45]+1.0*p1_over_gamma[16]*f[42]+1.0*p1_over_gamma[15]*f[41]+1.0*p1_over_gamma[14]*f[39]+1.0*p1_over_gamma[13]*f[36]+1.0*p1_over_gamma[12]*f[34]+1.0*p1_over_gamma[11]*f[33]+p1_over_gamma[10]*f[31]+1.0*p1_over_gamma[9]*f[28]+1.0*p1_over_gamma[8]*f[23]+1.0*p1_over_gamma[7]*f[20]+p1_over_gamma[6]*f[17]+p1_over_gamma[5]*f[16]+p1_over_gamma[4]*f[15]+p1_over_gamma[3]*f[8]+p1_over_gamma[2]*f[6]+p1_over_gamma[1]*f[5]+p1_over_gamma[0]*f[1])*volFact; 
  out[5] += (1.0*p1_over_gamma[10]*f[44]+p1_over_gamma[6]*f[37]+p1_over_gamma[5]*f[35]+p1_over_gamma[4]*f[32]+1.0*p1_over_gamma[3]*f[25]+1.0*p1_over_gamma[2]*f[21]+1.0*p1_over_gamma[1]*f[19]+p1_over_gamma[0]*f[11])*volFact; 
  out[6] += (p2_over_gamma[19]*f[43]+p2_over_gamma[18]*f[40]+p2_over_gamma[17]*f[38]+p2_over_gamma[16]*f[30]+p2_over_gamma[15]*f[29]+p2_over_gamma[14]*f[27]+p2_over_gamma[13]*f[26]+p2_over_gamma[12]*f[24]+p2_over_gamma[11]*f[22]+p2_over_gamma[10]*f[18]+p2_over_gamma[9]*f[14]+p2_over_gamma[8]*f[13]+p2_over_gamma[7]*f[12]+p2_over_gamma[6]*f[10]+p2_over_gamma[5]*f[9]+p2_over_gamma[4]*f[7]+p2_over_gamma[3]*f[4]+p2_over_gamma[2]*f[3]+p2_over_gamma[1]*f[2]+f[0]*p2_over_gamma[0])*volFact; 
  out[7] += (1.0*p2_over_gamma[19]*f[47]+1.0*p2_over_gamma[18]*f[46]+1.0*p2_over_gamma[17]*f[45]+1.0*p2_over_gamma[16]*f[42]+1.0*p2_over_gamma[15]*f[41]+1.0*p2_over_gamma[14]*f[39]+1.0*p2_over_gamma[13]*f[36]+1.0*p2_over_gamma[12]*f[34]+1.0*p2_over_gamma[11]*f[33]+p2_over_gamma[10]*f[31]+1.0*p2_over_gamma[9]*f[28]+1.0*p2_over_gamma[8]*f[23]+1.0*p2_over_gamma[7]*f[20]+p2_over_gamma[6]*f[17]+p2_over_gamma[5]*f[16]+p2_over_gamma[4]*f[15]+p2_over_gamma[3]*f[8]+p2_over_gamma[2]*f[6]+p2_over_gamma[1]*f[5]+p2_over_gamma[0]*f[1])*volFact; 
  out[8] += (1.0*p2_over_gamma[10]*f[44]+p2_over_gamma[6]*f[37]+p2_over_gamma[5]*f[35]+p2_over_gamma[4]*f[32]+1.0*p2_over_gamma[3]*f[25]+1.0*p2_over_gamma[2]*f[21]+1.0*p2_over_gamma[1]*f[19]+p2_over_gamma[0]*f[11])*volFact; 
} 
