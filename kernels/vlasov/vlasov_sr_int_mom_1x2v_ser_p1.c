#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_int_mom_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double *p1_over_gamma = &p_over_gamma[4]; 
 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += (1.414213562373095*p0_over_gamma[3]*f[6]+1.414213562373095*p0_over_gamma[2]*f[3]+1.414213562373095*p0_over_gamma[1]*f[2]+1.414213562373095*f[0]*p0_over_gamma[0])*volFact; 
  out[2] += (1.414213562373095*p1_over_gamma[3]*f[6]+1.414213562373095*p1_over_gamma[2]*f[3]+1.414213562373095*p1_over_gamma[1]*f[2]+1.414213562373095*f[0]*p1_over_gamma[0])*volFact; 
  out[3] += volFact*(1.414213562373095*p1_over_gamma[3]*f[6]*wx2+1.414213562373095*p1_over_gamma[2]*f[3]*wx2+1.414213562373095*p1_over_gamma[1]*f[2]*wx2+1.414213562373095*f[0]*p1_over_gamma[0]*wx2+1.414213562373095*p0_over_gamma[3]*f[6]*wx1+1.414213562373095*p0_over_gamma[2]*f[3]*wx1+1.414213562373095*p0_over_gamma[1]*f[2]*wx1+1.414213562373095*f[0]*p0_over_gamma[0]*wx1+0.408248290463863*p1_over_gamma[1]*f[6]*dv2+0.408248290463863*f[2]*p1_over_gamma[3]*dv2+0.408248290463863*p1_over_gamma[0]*f[3]*dv2+0.408248290463863*f[0]*p1_over_gamma[2]*dv2+0.408248290463863*p0_over_gamma[2]*f[6]*dv1+0.408248290463863*f[3]*p0_over_gamma[3]*dv1+0.408248290463863*p0_over_gamma[0]*f[2]*dv1+0.408248290463863*f[0]*p0_over_gamma[1]*dv1); 
} 
