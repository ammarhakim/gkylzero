#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_stream_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
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
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[64] = {0.0}; 
  alpha_cdim[0] = 2.0*p0_over_gamma[0]*dx10; 
  alpha_cdim[3] = 2.0*p0_over_gamma[1]*dx10; 
  alpha_cdim[4] = 2.0*p0_over_gamma[2]*dx10; 
  alpha_cdim[5] = 2.0*p0_over_gamma[3]*dx10; 
  alpha_cdim[11] = 2.0*p0_over_gamma[4]*dx10; 
  alpha_cdim[14] = 2.0*p0_over_gamma[5]*dx10; 
  alpha_cdim[15] = 2.0*p0_over_gamma[6]*dx10; 
  alpha_cdim[25] = 2.0*p0_over_gamma[7]*dx10; 
  cflFreq_mid += fabs(0.0883883476483184*alpha_cdim[0]); 

  alpha_cdim[32] = 2.0*p1_over_gamma[0]*dx11; 
  alpha_cdim[35] = 2.0*p1_over_gamma[1]*dx11; 
  alpha_cdim[36] = 2.0*p1_over_gamma[2]*dx11; 
  alpha_cdim[37] = 2.0*p1_over_gamma[3]*dx11; 
  alpha_cdim[43] = 2.0*p1_over_gamma[4]*dx11; 
  alpha_cdim[46] = 2.0*p1_over_gamma[5]*dx11; 
  alpha_cdim[47] = 2.0*p1_over_gamma[6]*dx11; 
  alpha_cdim[57] = 2.0*p1_over_gamma[7]*dx11; 
  cflFreq_mid += fabs(0.0883883476483184*alpha_cdim[32]); 

  out[1] += 0.3061862178478971*(alpha_cdim[25]*f[25]+alpha_cdim[15]*f[15]+alpha_cdim[14]*f[14]+alpha_cdim[11]*f[11]+alpha_cdim[5]*f[5]+alpha_cdim[4]*f[4]+alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.3061862178478971*(f[25]*alpha_cdim[57]+f[15]*alpha_cdim[47]+f[14]*alpha_cdim[46]+f[11]*alpha_cdim[43]+f[5]*alpha_cdim[37]+f[4]*alpha_cdim[36]+f[3]*alpha_cdim[35]+f[0]*alpha_cdim[32]); 
  out[6] += 0.3061862178478971*(f[29]*alpha_cdim[57]+f[23]*alpha_cdim[47]+f[21]*alpha_cdim[46]+f[18]*alpha_cdim[43]+f[12]*alpha_cdim[37]+f[9]*alpha_cdim[36]+f[7]*alpha_cdim[35]+f[1]*alpha_cdim[32]+alpha_cdim[25]*f[30]+alpha_cdim[15]*f[24]+alpha_cdim[14]*f[22]+alpha_cdim[11]*f[19]+alpha_cdim[5]*f[13]+alpha_cdim[4]*f[10]+alpha_cdim[3]*f[8]+alpha_cdim[0]*f[2]); 
  out[7] += 0.3061862178478971*(alpha_cdim[15]*f[25]+f[15]*alpha_cdim[25]+alpha_cdim[5]*f[14]+f[5]*alpha_cdim[14]+alpha_cdim[4]*f[11]+f[4]*alpha_cdim[11]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]); 
  out[8] += 0.3061862178478971*(f[15]*alpha_cdim[57]+f[25]*alpha_cdim[47]+f[5]*alpha_cdim[46]+f[4]*alpha_cdim[43]+f[14]*alpha_cdim[37]+f[11]*alpha_cdim[36]+f[0]*alpha_cdim[35]+f[3]*alpha_cdim[32]); 
  out[9] += 0.3061862178478971*(alpha_cdim[14]*f[25]+f[14]*alpha_cdim[25]+alpha_cdim[5]*f[15]+f[5]*alpha_cdim[15]+alpha_cdim[3]*f[11]+f[3]*alpha_cdim[11]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]); 
  out[10] += 0.3061862178478971*(f[14]*alpha_cdim[57]+f[5]*alpha_cdim[47]+f[25]*alpha_cdim[46]+f[3]*alpha_cdim[43]+f[15]*alpha_cdim[37]+f[0]*alpha_cdim[36]+f[11]*alpha_cdim[35]+f[4]*alpha_cdim[32]); 
  out[12] += 0.3061862178478971*(alpha_cdim[11]*f[25]+f[11]*alpha_cdim[25]+alpha_cdim[4]*f[15]+f[4]*alpha_cdim[15]+alpha_cdim[3]*f[14]+f[3]*alpha_cdim[14]+alpha_cdim[0]*f[5]+f[0]*alpha_cdim[5]); 
  out[13] += 0.3061862178478971*(f[11]*alpha_cdim[57]+f[4]*alpha_cdim[47]+f[3]*alpha_cdim[46]+f[25]*alpha_cdim[43]+f[0]*alpha_cdim[37]+f[15]*alpha_cdim[36]+f[14]*alpha_cdim[35]+f[5]*alpha_cdim[32]); 
  out[16] += 0.3061862178478971*(f[23]*alpha_cdim[57]+f[29]*alpha_cdim[47]+f[12]*alpha_cdim[46]+f[9]*alpha_cdim[43]+f[21]*alpha_cdim[37]+f[18]*alpha_cdim[36]+f[1]*alpha_cdim[35]+f[7]*alpha_cdim[32]+alpha_cdim[15]*f[30]+f[24]*alpha_cdim[25]+alpha_cdim[5]*f[22]+alpha_cdim[4]*f[19]+f[13]*alpha_cdim[14]+f[10]*alpha_cdim[11]+alpha_cdim[0]*f[8]+f[2]*alpha_cdim[3]); 
  out[17] += 0.3061862178478971*(f[21]*alpha_cdim[57]+f[12]*alpha_cdim[47]+f[29]*alpha_cdim[46]+f[7]*alpha_cdim[43]+f[23]*alpha_cdim[37]+f[1]*alpha_cdim[36]+f[18]*alpha_cdim[35]+f[9]*alpha_cdim[32]+alpha_cdim[14]*f[30]+f[22]*alpha_cdim[25]+alpha_cdim[5]*f[24]+alpha_cdim[3]*f[19]+f[13]*alpha_cdim[15]+f[8]*alpha_cdim[11]+alpha_cdim[0]*f[10]+f[2]*alpha_cdim[4]); 
  out[18] += 0.3061862178478971*(alpha_cdim[5]*f[25]+f[5]*alpha_cdim[25]+alpha_cdim[14]*f[15]+f[14]*alpha_cdim[15]+alpha_cdim[0]*f[11]+f[0]*alpha_cdim[11]+alpha_cdim[3]*f[4]+f[3]*alpha_cdim[4]); 
  out[19] += 0.3061862178478971*(f[5]*alpha_cdim[57]+f[14]*alpha_cdim[47]+f[15]*alpha_cdim[46]+f[0]*alpha_cdim[43]+f[25]*alpha_cdim[37]+f[3]*alpha_cdim[36]+f[4]*alpha_cdim[35]+f[11]*alpha_cdim[32]); 
  out[20] += 0.3061862178478971*(f[18]*alpha_cdim[57]+f[9]*alpha_cdim[47]+f[7]*alpha_cdim[46]+f[29]*alpha_cdim[43]+f[1]*alpha_cdim[37]+f[23]*alpha_cdim[36]+f[21]*alpha_cdim[35]+f[12]*alpha_cdim[32]+alpha_cdim[11]*f[30]+f[19]*alpha_cdim[25]+alpha_cdim[4]*f[24]+alpha_cdim[3]*f[22]+f[10]*alpha_cdim[15]+f[8]*alpha_cdim[14]+alpha_cdim[0]*f[13]+f[2]*alpha_cdim[5]); 
  out[21] += 0.3061862178478971*(alpha_cdim[4]*f[25]+f[4]*alpha_cdim[25]+alpha_cdim[11]*f[15]+f[11]*alpha_cdim[15]+alpha_cdim[0]*f[14]+f[0]*alpha_cdim[14]+alpha_cdim[3]*f[5]+f[3]*alpha_cdim[5]); 
  out[22] += 0.3061862178478971*(f[4]*alpha_cdim[57]+f[11]*alpha_cdim[47]+f[0]*alpha_cdim[46]+f[15]*alpha_cdim[43]+f[3]*alpha_cdim[37]+f[25]*alpha_cdim[36]+f[5]*alpha_cdim[35]+f[14]*alpha_cdim[32]); 
  out[23] += 0.3061862178478971*(alpha_cdim[3]*f[25]+f[3]*alpha_cdim[25]+alpha_cdim[0]*f[15]+f[0]*alpha_cdim[15]+alpha_cdim[11]*f[14]+f[11]*alpha_cdim[14]+alpha_cdim[4]*f[5]+f[4]*alpha_cdim[5]); 
  out[24] += 0.3061862178478971*(f[3]*alpha_cdim[57]+f[0]*alpha_cdim[47]+f[11]*alpha_cdim[46]+f[14]*alpha_cdim[43]+f[4]*alpha_cdim[37]+f[5]*alpha_cdim[36]+f[25]*alpha_cdim[35]+f[15]*alpha_cdim[32]); 
  out[26] += 0.3061862178478971*(f[12]*alpha_cdim[57]+f[21]*alpha_cdim[47]+f[23]*alpha_cdim[46]+f[1]*alpha_cdim[43]+f[29]*alpha_cdim[37]+f[7]*alpha_cdim[36]+f[9]*alpha_cdim[35]+f[18]*alpha_cdim[32]+alpha_cdim[5]*f[30]+f[13]*alpha_cdim[25]+alpha_cdim[14]*f[24]+alpha_cdim[15]*f[22]+alpha_cdim[0]*f[19]+f[2]*alpha_cdim[11]+alpha_cdim[3]*f[10]+alpha_cdim[4]*f[8]); 
  out[27] += 0.3061862178478971*(f[9]*alpha_cdim[57]+f[18]*alpha_cdim[47]+f[1]*alpha_cdim[46]+f[23]*alpha_cdim[43]+f[7]*alpha_cdim[37]+f[29]*alpha_cdim[36]+f[12]*alpha_cdim[35]+f[21]*alpha_cdim[32]+alpha_cdim[4]*f[30]+f[10]*alpha_cdim[25]+alpha_cdim[11]*f[24]+alpha_cdim[0]*f[22]+alpha_cdim[15]*f[19]+f[2]*alpha_cdim[14]+alpha_cdim[3]*f[13]+alpha_cdim[5]*f[8]); 
  out[28] += 0.3061862178478971*(f[7]*alpha_cdim[57]+f[1]*alpha_cdim[47]+f[18]*alpha_cdim[46]+f[21]*alpha_cdim[43]+f[9]*alpha_cdim[37]+f[12]*alpha_cdim[36]+f[29]*alpha_cdim[35]+f[23]*alpha_cdim[32]+alpha_cdim[3]*f[30]+f[8]*alpha_cdim[25]+alpha_cdim[0]*f[24]+alpha_cdim[11]*f[22]+alpha_cdim[14]*f[19]+f[2]*alpha_cdim[15]+alpha_cdim[4]*f[13]+alpha_cdim[5]*f[10]); 
  out[29] += 0.3061862178478971*(alpha_cdim[0]*f[25]+f[0]*alpha_cdim[25]+alpha_cdim[3]*f[15]+f[3]*alpha_cdim[15]+alpha_cdim[4]*f[14]+f[4]*alpha_cdim[14]+alpha_cdim[5]*f[11]+f[5]*alpha_cdim[11]); 
  out[30] += 0.3061862178478971*(f[0]*alpha_cdim[57]+f[3]*alpha_cdim[47]+f[4]*alpha_cdim[46]+f[5]*alpha_cdim[43]+f[11]*alpha_cdim[37]+f[14]*alpha_cdim[36]+f[15]*alpha_cdim[35]+f[25]*alpha_cdim[32]); 
  out[31] += 0.3061862178478971*(f[1]*alpha_cdim[57]+f[7]*alpha_cdim[47]+f[9]*alpha_cdim[46]+f[12]*alpha_cdim[43]+f[18]*alpha_cdim[37]+f[21]*alpha_cdim[36]+f[23]*alpha_cdim[35]+f[29]*alpha_cdim[32]+alpha_cdim[0]*f[30]+f[2]*alpha_cdim[25]+alpha_cdim[3]*f[24]+alpha_cdim[4]*f[22]+alpha_cdim[5]*f[19]+f[8]*alpha_cdim[15]+f[10]*alpha_cdim[14]+alpha_cdim[11]*f[13]); 

  return 3.0*cflFreq_mid; 
} 
