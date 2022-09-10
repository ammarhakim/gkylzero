#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_sr_stream_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // qmem:      q/m*EM fields (unused in streaming-only update).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *p1_over_gamma = &p_over_gamma[4]; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[32] = {0.0}; 
  alpha_cdim[0] = 2.0*p0_over_gamma[0]*dx10; 
  alpha_cdim[3] = 2.0*p0_over_gamma[1]*dx10; 
  alpha_cdim[4] = 2.0*p0_over_gamma[2]*dx10; 
  alpha_cdim[10] = 2.0*p0_over_gamma[3]*dx10; 
  cflFreq_mid += fabs(0.125*alpha_cdim[0]); 

  alpha_cdim[16] = 2.0*p1_over_gamma[0]*dx11; 
  alpha_cdim[19] = 2.0*p1_over_gamma[1]*dx11; 
  alpha_cdim[20] = 2.0*p1_over_gamma[2]*dx11; 
  alpha_cdim[26] = 2.0*p1_over_gamma[3]*dx11; 
  cflFreq_mid += fabs(0.125*alpha_cdim[16]); 

  out[1] += 0.4330127018922193*(alpha_cdim[10]*f[10]+alpha_cdim[4]*f[4]+alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(f[10]*alpha_cdim[26]+f[4]*alpha_cdim[20]+f[3]*alpha_cdim[19]+f[0]*alpha_cdim[16]); 
  out[5] += 0.4330127018922193*(f[13]*alpha_cdim[26]+f[8]*alpha_cdim[20]+f[6]*alpha_cdim[19]+f[1]*alpha_cdim[16]+alpha_cdim[10]*f[14]+alpha_cdim[4]*f[9]+alpha_cdim[3]*f[7]+alpha_cdim[0]*f[2]); 
  out[6] += 0.4330127018922193*(alpha_cdim[4]*f[10]+f[4]*alpha_cdim[10]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]); 
  out[7] += 0.4330127018922193*(f[4]*alpha_cdim[26]+f[10]*alpha_cdim[20]+f[0]*alpha_cdim[19]+f[3]*alpha_cdim[16]); 
  out[8] += 0.4330127018922193*(alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]); 
  out[9] += 0.4330127018922193*(f[3]*alpha_cdim[26]+f[0]*alpha_cdim[20]+f[10]*alpha_cdim[19]+f[4]*alpha_cdim[16]); 
  out[11] += 0.4330127018922193*(f[8]*alpha_cdim[26]+f[13]*alpha_cdim[20]+f[1]*alpha_cdim[19]+f[6]*alpha_cdim[16]+alpha_cdim[4]*f[14]+f[9]*alpha_cdim[10]+alpha_cdim[0]*f[7]+f[2]*alpha_cdim[3]); 
  out[12] += 0.4330127018922193*(f[6]*alpha_cdim[26]+f[1]*alpha_cdim[20]+f[13]*alpha_cdim[19]+f[8]*alpha_cdim[16]+alpha_cdim[3]*f[14]+f[7]*alpha_cdim[10]+alpha_cdim[0]*f[9]+f[2]*alpha_cdim[4]); 
  out[13] += 0.4330127018922193*(alpha_cdim[0]*f[10]+f[0]*alpha_cdim[10]+alpha_cdim[3]*f[4]+f[3]*alpha_cdim[4]); 
  out[14] += 0.4330127018922193*(f[0]*alpha_cdim[26]+f[3]*alpha_cdim[20]+f[4]*alpha_cdim[19]+f[10]*alpha_cdim[16]); 
  out[15] += 0.4330127018922193*(f[1]*alpha_cdim[26]+f[6]*alpha_cdim[20]+f[8]*alpha_cdim[19]+f[13]*alpha_cdim[16]+alpha_cdim[0]*f[14]+f[2]*alpha_cdim[10]+alpha_cdim[3]*f[9]+alpha_cdim[4]*f[7]); 

  return 3.0*cflFreq_mid; 
} 
