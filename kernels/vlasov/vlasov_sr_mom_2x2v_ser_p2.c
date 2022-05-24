#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 2.0*f[11]*volFact; 
  out[5] += 2.0*f[12]*volFact; 
  out[6] += 2.0*f[19]*volFact; 
  out[7] += 2.0*f[20]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  out[0] += (p0_over_gamma[7]*f[30]+p0_over_gamma[6]*f[27]+p0_over_gamma[5]*f[14]+p0_over_gamma[4]*f[13]+p0_over_gamma[3]*f[10]+p0_over_gamma[2]*f[4]+p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (1.0*p0_over_gamma[7]*f[42]+1.0*p0_over_gamma[6]*f[39]+1.0*p0_over_gamma[5]*f[28]+1.0*p0_over_gamma[4]*f[23]+p0_over_gamma[3]*f[17]+p0_over_gamma[2]*f[8]+p0_over_gamma[1]*f[6]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (1.0*p0_over_gamma[7]*f[43]+1.0*p0_over_gamma[6]*f[40]+1.0*p0_over_gamma[5]*f[29]+1.0*p0_over_gamma[4]*f[24]+p0_over_gamma[3]*f[18]+p0_over_gamma[2]*f[9]+p0_over_gamma[1]*f[7]+p0_over_gamma[0]*f[2])*volFact; 
  out[3] += (p0_over_gamma[7]*f[47]+p0_over_gamma[6]*f[46]+p0_over_gamma[5]*f[41]+p0_over_gamma[4]*f[34]+p0_over_gamma[3]*f[31]+p0_over_gamma[2]*f[16]+p0_over_gamma[1]*f[15]+p0_over_gamma[0]*f[5])*volFact; 
  out[4] += (p0_over_gamma[3]*f[37]+1.0*p0_over_gamma[2]*f[25]+1.0*p0_over_gamma[1]*f[21]+p0_over_gamma[0]*f[11])*volFact; 
  out[5] += (p0_over_gamma[3]*f[38]+1.0*p0_over_gamma[2]*f[26]+1.0*p0_over_gamma[1]*f[22]+p0_over_gamma[0]*f[12])*volFact; 
  out[6] += (p0_over_gamma[3]*f[44]+1.0*p0_over_gamma[2]*f[35]+1.0*p0_over_gamma[1]*f[32]+p0_over_gamma[0]*f[19])*volFact; 
  out[7] += (p0_over_gamma[3]*f[45]+1.0*p0_over_gamma[2]*f[36]+1.0*p0_over_gamma[1]*f[33]+p0_over_gamma[0]*f[20])*volFact; 
  out[8] += (p1_over_gamma[7]*f[30]+p1_over_gamma[6]*f[27]+p1_over_gamma[5]*f[14]+p1_over_gamma[4]*f[13]+p1_over_gamma[3]*f[10]+p1_over_gamma[2]*f[4]+p1_over_gamma[1]*f[3]+f[0]*p1_over_gamma[0])*volFact; 
  out[9] += (1.0*p1_over_gamma[7]*f[42]+1.0*p1_over_gamma[6]*f[39]+1.0*p1_over_gamma[5]*f[28]+1.0*p1_over_gamma[4]*f[23]+p1_over_gamma[3]*f[17]+p1_over_gamma[2]*f[8]+p1_over_gamma[1]*f[6]+p1_over_gamma[0]*f[1])*volFact; 
  out[10] += (1.0*p1_over_gamma[7]*f[43]+1.0*p1_over_gamma[6]*f[40]+1.0*p1_over_gamma[5]*f[29]+1.0*p1_over_gamma[4]*f[24]+p1_over_gamma[3]*f[18]+p1_over_gamma[2]*f[9]+p1_over_gamma[1]*f[7]+p1_over_gamma[0]*f[2])*volFact; 
  out[11] += (p1_over_gamma[7]*f[47]+p1_over_gamma[6]*f[46]+p1_over_gamma[5]*f[41]+p1_over_gamma[4]*f[34]+p1_over_gamma[3]*f[31]+p1_over_gamma[2]*f[16]+p1_over_gamma[1]*f[15]+p1_over_gamma[0]*f[5])*volFact; 
  out[12] += (p1_over_gamma[3]*f[37]+1.0*p1_over_gamma[2]*f[25]+1.0*p1_over_gamma[1]*f[21]+p1_over_gamma[0]*f[11])*volFact; 
  out[13] += (p1_over_gamma[3]*f[38]+1.0*p1_over_gamma[2]*f[26]+1.0*p1_over_gamma[1]*f[22]+p1_over_gamma[0]*f[12])*volFact; 
  out[14] += (p1_over_gamma[3]*f[44]+1.0*p1_over_gamma[2]*f[35]+1.0*p1_over_gamma[1]*f[32]+p1_over_gamma[0]*f[19])*volFact; 
  out[15] += (p1_over_gamma[3]*f[45]+1.0*p1_over_gamma[2]*f[36]+1.0*p1_over_gamma[1]*f[33]+p1_over_gamma[0]*f[20])*volFact; 
} 
