#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_int_mom_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *gamma, const double *GammaV_inv, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
 
  out[0] += (2.0*GammaV_inv[2]*f[7]+2.0*GammaV_inv[1]*f[1]+2.0*GammaV_inv[0]*f[0])*volFact; 
  out[1] += (1.0*GammaV_inv[1]*p0_over_gamma[7]*f[19]+1.0*GammaV_inv[1]*p0_over_gamma[6]*f[18]+GammaV_inv[2]*p0_over_gamma[3]*f[17]+GammaV_inv[0]*p0_over_gamma[7]*f[16]+1.0*GammaV_inv[1]*p0_over_gamma[5]*f[15]+GammaV_inv[0]*p0_over_gamma[6]*f[14]+1.0*GammaV_inv[2]*p0_over_gamma[2]*f[13]+1.0*GammaV_inv[1]*p0_over_gamma[4]*f[12]+1.0*p0_over_gamma[1]*GammaV_inv[2]*f[11]+GammaV_inv[1]*p0_over_gamma[3]*f[10]+GammaV_inv[0]*p0_over_gamma[5]*f[9]+GammaV_inv[0]*p0_over_gamma[4]*f[8]+p0_over_gamma[0]*GammaV_inv[2]*f[7]+GammaV_inv[0]*p0_over_gamma[3]*f[6]+GammaV_inv[1]*p0_over_gamma[2]*f[5]+GammaV_inv[1]*p0_over_gamma[1]*f[4]+GammaV_inv[0]*p0_over_gamma[2]*f[3]+GammaV_inv[0]*p0_over_gamma[1]*f[2]+p0_over_gamma[0]*GammaV_inv[1]*f[1]+GammaV_inv[0]*f[0]*p0_over_gamma[0])*volFact; 
  out[2] += (1.0*GammaV_inv[1]*p1_over_gamma[7]*f[19]+1.0*GammaV_inv[1]*p1_over_gamma[6]*f[18]+GammaV_inv[2]*p1_over_gamma[3]*f[17]+GammaV_inv[0]*p1_over_gamma[7]*f[16]+1.0*GammaV_inv[1]*p1_over_gamma[5]*f[15]+GammaV_inv[0]*p1_over_gamma[6]*f[14]+1.0*GammaV_inv[2]*p1_over_gamma[2]*f[13]+1.0*GammaV_inv[1]*p1_over_gamma[4]*f[12]+1.0*p1_over_gamma[1]*GammaV_inv[2]*f[11]+GammaV_inv[1]*p1_over_gamma[3]*f[10]+GammaV_inv[0]*p1_over_gamma[5]*f[9]+GammaV_inv[0]*p1_over_gamma[4]*f[8]+p1_over_gamma[0]*GammaV_inv[2]*f[7]+GammaV_inv[0]*p1_over_gamma[3]*f[6]+GammaV_inv[1]*p1_over_gamma[2]*f[5]+GammaV_inv[1]*p1_over_gamma[1]*f[4]+GammaV_inv[0]*p1_over_gamma[2]*f[3]+GammaV_inv[0]*p1_over_gamma[1]*f[2]+p1_over_gamma[0]*GammaV_inv[1]*f[1]+GammaV_inv[0]*f[0]*p1_over_gamma[0])*volFact; 
  out[3] += (1.0*GammaV_inv[1]*gamma[7]*f[19]+1.0*GammaV_inv[1]*gamma[6]*f[18]+GammaV_inv[2]*gamma[3]*f[17]+GammaV_inv[0]*gamma[7]*f[16]+1.0*GammaV_inv[1]*gamma[5]*f[15]+GammaV_inv[0]*gamma[6]*f[14]+1.0*GammaV_inv[2]*gamma[2]*f[13]+1.0*GammaV_inv[1]*gamma[4]*f[12]+1.0*gamma[1]*GammaV_inv[2]*f[11]+GammaV_inv[1]*gamma[3]*f[10]+GammaV_inv[0]*gamma[5]*f[9]+GammaV_inv[0]*gamma[4]*f[8]+gamma[0]*GammaV_inv[2]*f[7]+GammaV_inv[0]*gamma[3]*f[6]+GammaV_inv[1]*gamma[2]*f[5]+GammaV_inv[1]*gamma[1]*f[4]+GammaV_inv[0]*gamma[2]*f[3]+GammaV_inv[0]*gamma[1]*f[2]+gamma[0]*GammaV_inv[1]*f[1]+GammaV_inv[0]*f[0]*gamma[0])*volFact; 
} 
