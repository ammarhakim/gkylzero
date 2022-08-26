#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_sr_stream_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // qmem:      q/m*EM fields (unused in streaming-only update).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2/dxv[0]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[4] = {0.0}; 
  alpha_cdim[0] = 1.414213562373095*p0_over_gamma[0]*dx10; 
  alpha_cdim[2] = 1.414213562373095*p0_over_gamma[1]*dx10; 
  cflFreq_mid += fabs(0.25*alpha_cdim[0]); 

  out[1] += 0.8660254037844386*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[3] += 0.8660254037844386*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]); 

  return 3.0*cflFreq_mid; 
} 
