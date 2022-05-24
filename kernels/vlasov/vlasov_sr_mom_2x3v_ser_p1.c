#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[6]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 
  out[0] += (p0_over_gamma[7]*f[25]+p0_over_gamma[6]*f[15]+p0_over_gamma[5]*f[14]+p0_over_gamma[4]*f[11]+p0_over_gamma[3]*f[5]+p0_over_gamma[2]*f[4]+p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (p0_over_gamma[7]*f[29]+p0_over_gamma[6]*f[23]+p0_over_gamma[5]*f[21]+p0_over_gamma[4]*f[18]+p0_over_gamma[3]*f[12]+p0_over_gamma[2]*f[9]+p0_over_gamma[1]*f[7]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (p0_over_gamma[7]*f[30]+p0_over_gamma[6]*f[24]+p0_over_gamma[5]*f[22]+p0_over_gamma[4]*f[19]+p0_over_gamma[3]*f[13]+p0_over_gamma[2]*f[10]+p0_over_gamma[1]*f[8]+p0_over_gamma[0]*f[2])*volFact; 
  out[3] += (p0_over_gamma[7]*f[31]+p0_over_gamma[6]*f[28]+p0_over_gamma[5]*f[27]+p0_over_gamma[4]*f[26]+p0_over_gamma[3]*f[20]+p0_over_gamma[2]*f[17]+p0_over_gamma[1]*f[16]+p0_over_gamma[0]*f[6])*volFact; 
  out[4] += (p1_over_gamma[7]*f[25]+p1_over_gamma[6]*f[15]+p1_over_gamma[5]*f[14]+p1_over_gamma[4]*f[11]+p1_over_gamma[3]*f[5]+p1_over_gamma[2]*f[4]+p1_over_gamma[1]*f[3]+f[0]*p1_over_gamma[0])*volFact; 
  out[5] += (p1_over_gamma[7]*f[29]+p1_over_gamma[6]*f[23]+p1_over_gamma[5]*f[21]+p1_over_gamma[4]*f[18]+p1_over_gamma[3]*f[12]+p1_over_gamma[2]*f[9]+p1_over_gamma[1]*f[7]+p1_over_gamma[0]*f[1])*volFact; 
  out[6] += (p1_over_gamma[7]*f[30]+p1_over_gamma[6]*f[24]+p1_over_gamma[5]*f[22]+p1_over_gamma[4]*f[19]+p1_over_gamma[3]*f[13]+p1_over_gamma[2]*f[10]+p1_over_gamma[1]*f[8]+p1_over_gamma[0]*f[2])*volFact; 
  out[7] += (p1_over_gamma[7]*f[31]+p1_over_gamma[6]*f[28]+p1_over_gamma[5]*f[27]+p1_over_gamma[4]*f[26]+p1_over_gamma[3]*f[20]+p1_over_gamma[2]*f[17]+p1_over_gamma[1]*f[16]+p1_over_gamma[0]*f[6])*volFact; 
  out[8] += (p2_over_gamma[7]*f[25]+p2_over_gamma[6]*f[15]+p2_over_gamma[5]*f[14]+p2_over_gamma[4]*f[11]+p2_over_gamma[3]*f[5]+p2_over_gamma[2]*f[4]+p2_over_gamma[1]*f[3]+f[0]*p2_over_gamma[0])*volFact; 
  out[9] += (p2_over_gamma[7]*f[29]+p2_over_gamma[6]*f[23]+p2_over_gamma[5]*f[21]+p2_over_gamma[4]*f[18]+p2_over_gamma[3]*f[12]+p2_over_gamma[2]*f[9]+p2_over_gamma[1]*f[7]+p2_over_gamma[0]*f[1])*volFact; 
  out[10] += (p2_over_gamma[7]*f[30]+p2_over_gamma[6]*f[24]+p2_over_gamma[5]*f[22]+p2_over_gamma[4]*f[19]+p2_over_gamma[3]*f[13]+p2_over_gamma[2]*f[10]+p2_over_gamma[1]*f[8]+p2_over_gamma[0]*f[2])*volFact; 
  out[11] += (p2_over_gamma[7]*f[31]+p2_over_gamma[6]*f[28]+p2_over_gamma[5]*f[27]+p2_over_gamma[4]*f[26]+p2_over_gamma[3]*f[20]+p2_over_gamma[2]*f[17]+p2_over_gamma[1]*f[16]+p2_over_gamma[0]*f[6])*volFact; 
} 
