#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_int_mom_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*0.25; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += (1.414213562373095*p0_over_gamma[1]*f[2]+1.414213562373095*f[0]*p0_over_gamma[0])*volFact; 
  out[2] += volFact*(1.414213562373095*p0_over_gamma[1]*f[2]*wx1+1.414213562373095*f[0]*p0_over_gamma[0]*wx1+0.408248290463863*p0_over_gamma[0]*f[2]*dv1+0.408248290463863*f[0]*p0_over_gamma[1]*dv1); 
} 
