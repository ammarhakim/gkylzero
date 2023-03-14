#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_int_mom_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *gamma, const double *GammaV_inv, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*0.25; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
 
  out[0] += (1.414213562373095*GammaV_inv[2]*f[4]+1.414213562373095*GammaV_inv[1]*f[1]+1.414213562373095*GammaV_inv[0]*f[0])*volFact; 
  out[1] += (1.0*GammaV_inv[1]*p0_over_gamma[2]*f[7]+1.0*p0_over_gamma[1]*GammaV_inv[2]*f[6]+GammaV_inv[0]*p0_over_gamma[2]*f[5]+p0_over_gamma[0]*GammaV_inv[2]*f[4]+GammaV_inv[1]*p0_over_gamma[1]*f[3]+GammaV_inv[0]*p0_over_gamma[1]*f[2]+p0_over_gamma[0]*GammaV_inv[1]*f[1]+GammaV_inv[0]*f[0]*p0_over_gamma[0])*volFact; 
  out[2] += (1.0*GammaV_inv[1]*gamma[2]*f[7]+1.0*gamma[1]*GammaV_inv[2]*f[6]+GammaV_inv[0]*gamma[2]*f[5]+gamma[0]*GammaV_inv[2]*f[4]+GammaV_inv[1]*gamma[1]*f[3]+GammaV_inv[0]*gamma[1]*f[2]+gamma[0]*GammaV_inv[1]*f[1]+GammaV_inv[0]*f[0]*gamma[0])*volFact; 
} 
