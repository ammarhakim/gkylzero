#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_1x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 
  out[0] += (p0_over_gamma[7]*f[14]+p0_over_gamma[6]*f[10]+p0_over_gamma[5]*f[9]+p0_over_gamma[4]*f[7]+p0_over_gamma[3]*f[4]+p0_over_gamma[2]*f[3]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (p0_over_gamma[7]*f[15]+p0_over_gamma[6]*f[13]+p0_over_gamma[5]*f[12]+p0_over_gamma[4]*f[11]+p0_over_gamma[3]*f[8]+p0_over_gamma[2]*f[6]+p0_over_gamma[1]*f[5]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (p1_over_gamma[7]*f[14]+p1_over_gamma[6]*f[10]+p1_over_gamma[5]*f[9]+p1_over_gamma[4]*f[7]+p1_over_gamma[3]*f[4]+p1_over_gamma[2]*f[3]+p1_over_gamma[1]*f[2]+f[0]*p1_over_gamma[0])*volFact; 
  out[3] += (p1_over_gamma[7]*f[15]+p1_over_gamma[6]*f[13]+p1_over_gamma[5]*f[12]+p1_over_gamma[4]*f[11]+p1_over_gamma[3]*f[8]+p1_over_gamma[2]*f[6]+p1_over_gamma[1]*f[5]+p1_over_gamma[0]*f[1])*volFact; 
  out[4] += (p2_over_gamma[7]*f[14]+p2_over_gamma[6]*f[10]+p2_over_gamma[5]*f[9]+p2_over_gamma[4]*f[7]+p2_over_gamma[3]*f[4]+p2_over_gamma[2]*f[3]+p2_over_gamma[1]*f[2]+f[0]*p2_over_gamma[0])*volFact; 
  out[5] += (p2_over_gamma[7]*f[15]+p2_over_gamma[6]*f[13]+p2_over_gamma[5]*f[12]+p2_over_gamma[4]*f[11]+p2_over_gamma[3]*f[8]+p2_over_gamma[2]*f[6]+p2_over_gamma[1]*f[5]+p2_over_gamma[0]*f[1])*volFact; 
} 
