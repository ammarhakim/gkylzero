#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_sr_stream_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // qmem:      q/m*EM fields (unused in streaming-only update).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2/dxv[0]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *p1_over_gamma = &p_over_gamma[4]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[8] = {0.0}; 
  alpha_cdim[0] = 1.414213562373095*p0_over_gamma[0]*dx10; 
  alpha_cdim[2] = 1.414213562373095*p0_over_gamma[1]*dx10; 
  alpha_cdim[3] = 1.414213562373095*p0_over_gamma[2]*dx10; 
  alpha_cdim[6] = 1.414213562373095*p0_over_gamma[3]*dx10; 
  alpha_mid += fabs(0.1767766952966368*alpha_cdim[0]); 

  out[1] += 0.6123724356957944*(alpha_cdim[6]*f[6]+alpha_cdim[3]*f[3]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[4] += 0.6123724356957944*(alpha_cdim[3]*f[6]+f[3]*alpha_cdim[6]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]); 
  out[5] += 0.6123724356957944*(alpha_cdim[2]*f[6]+f[2]*alpha_cdim[6]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]); 
  out[7] += 0.6123724356957944*(alpha_cdim[0]*f[6]+f[0]*alpha_cdim[6]+alpha_cdim[2]*f[3]+f[2]*alpha_cdim[3]); 

  return alpha_mid; 
} 
