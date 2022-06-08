#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_M0_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
} 
GKYL_CU_DH void vlasov_sr_M1i_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *p1_over_gamma = &p_over_gamma[4]; 
  out[0] += (p0_over_gamma[3]*f[6]+p0_over_gamma[2]*f[3]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (p0_over_gamma[3]*f[7]+p0_over_gamma[2]*f[5]+p0_over_gamma[1]*f[4]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (p1_over_gamma[3]*f[6]+p1_over_gamma[2]*f[3]+p1_over_gamma[1]*f[2]+f[0]*p1_over_gamma[0])*volFact; 
  out[3] += (p1_over_gamma[3]*f[7]+p1_over_gamma[2]*f[5]+p1_over_gamma[1]*f[4]+p1_over_gamma[0]*f[1])*volFact; 
} 
