#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_int_mom_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double *p1_over_gamma = &p_over_gamma[4]; 
 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (2.0*p0_over_gamma[3]*f[10]+2.0*p0_over_gamma[2]*f[4]+2.0*p0_over_gamma[1]*f[3]+2.0*f[0]*p0_over_gamma[0])*volFact; 
  out[2] += (2.0*p1_over_gamma[3]*f[10]+2.0*p1_over_gamma[2]*f[4]+2.0*p1_over_gamma[1]*f[3]+2.0*f[0]*p1_over_gamma[0])*volFact; 
  out[3] += volFact*(2.0*p1_over_gamma[3]*f[10]*wx2+2.0*p1_over_gamma[2]*f[4]*wx2+2.0*p1_over_gamma[1]*f[3]*wx2+2.0*f[0]*p1_over_gamma[0]*wx2+2.0*p0_over_gamma[3]*f[10]*wx1+2.0*p0_over_gamma[2]*f[4]*wx1+2.0*p0_over_gamma[1]*f[3]*wx1+2.0*f[0]*p0_over_gamma[0]*wx1+0.5773502691896258*p1_over_gamma[1]*f[10]*dv2+0.5773502691896258*p1_over_gamma[0]*f[4]*dv2+0.5773502691896258*f[3]*p1_over_gamma[3]*dv2+0.5773502691896258*f[0]*p1_over_gamma[2]*dv2+0.5773502691896258*p0_over_gamma[2]*f[10]*dv1+0.5773502691896258*p0_over_gamma[3]*f[4]*dv1+0.5773502691896258*p0_over_gamma[0]*f[3]*dv1+0.5773502691896258*f[0]*p0_over_gamma[1]*dv1); 
} 
