#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_stream_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // qmem:      q/m*EM fields (unused in streaming-only update).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2/dxv[0]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[16] = {0.0}; 
  alpha_cdim[0] = 1.414213562373095*p0_over_gamma[0]*dx10; 
  alpha_cdim[2] = 1.414213562373095*p0_over_gamma[1]*dx10; 
  alpha_cdim[3] = 1.414213562373095*p0_over_gamma[2]*dx10; 
  alpha_cdim[4] = 1.414213562373095*p0_over_gamma[3]*dx10; 
  alpha_cdim[7] = 1.414213562373095*p0_over_gamma[4]*dx10; 
  alpha_cdim[9] = 1.414213562373095*p0_over_gamma[5]*dx10; 
  alpha_cdim[10] = 1.414213562373095*p0_over_gamma[6]*dx10; 
  alpha_cdim[14] = 1.414213562373095*p0_over_gamma[7]*dx10; 
  cflFreq_mid += fabs(0.125*alpha_cdim[0]); 

  out[1] += 0.4330127018922193*(alpha_cdim[14]*f[14]+alpha_cdim[10]*f[10]+alpha_cdim[9]*f[9]+alpha_cdim[7]*f[7]+alpha_cdim[4]*f[4]+alpha_cdim[3]*f[3]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[5] += 0.4330127018922193*(alpha_cdim[10]*f[14]+f[10]*alpha_cdim[14]+alpha_cdim[4]*f[9]+f[4]*alpha_cdim[9]+alpha_cdim[3]*f[7]+f[3]*alpha_cdim[7]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]); 
  out[6] += 0.4330127018922193*(alpha_cdim[9]*f[14]+f[9]*alpha_cdim[14]+alpha_cdim[4]*f[10]+f[4]*alpha_cdim[10]+alpha_cdim[2]*f[7]+f[2]*alpha_cdim[7]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]); 
  out[8] += 0.4330127018922193*(alpha_cdim[7]*f[14]+f[7]*alpha_cdim[14]+alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10]+alpha_cdim[2]*f[9]+f[2]*alpha_cdim[9]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]); 
  out[11] += 0.4330127018922193*(alpha_cdim[4]*f[14]+f[4]*alpha_cdim[14]+alpha_cdim[9]*f[10]+f[9]*alpha_cdim[10]+alpha_cdim[0]*f[7]+f[0]*alpha_cdim[7]+alpha_cdim[2]*f[3]+f[2]*alpha_cdim[3]); 
  out[12] += 0.4330127018922193*(alpha_cdim[3]*f[14]+f[3]*alpha_cdim[14]+alpha_cdim[7]*f[10]+f[7]*alpha_cdim[10]+alpha_cdim[0]*f[9]+f[0]*alpha_cdim[9]+alpha_cdim[2]*f[4]+f[2]*alpha_cdim[4]); 
  out[13] += 0.4330127018922193*(alpha_cdim[2]*f[14]+f[2]*alpha_cdim[14]+alpha_cdim[0]*f[10]+f[0]*alpha_cdim[10]+alpha_cdim[7]*f[9]+f[7]*alpha_cdim[9]+alpha_cdim[3]*f[4]+f[3]*alpha_cdim[4]); 
  out[15] += 0.4330127018922193*(alpha_cdim[0]*f[14]+f[0]*alpha_cdim[14]+alpha_cdim[2]*f[10]+f[2]*alpha_cdim[10]+alpha_cdim[3]*f[9]+f[3]*alpha_cdim[9]+alpha_cdim[4]*f[7]+f[4]*alpha_cdim[7]); 

  return 3.0*cflFreq_mid; 
} 
