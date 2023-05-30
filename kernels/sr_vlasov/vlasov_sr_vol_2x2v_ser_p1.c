#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // qmem:      q/m*EM fields.
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double dv10 = 2/dxv[2]; 
  const double *E0 = &qmem[0]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double dv11 = 2/dxv[3]; 
  const double *E1 = &qmem[4]; 
  const double *p1_over_gamma = &p_over_gamma[4]; 

  const double *B2 = &qmem[20]; 
  double cflFreq_mid = 0.0; 
  double alpha_cdim[32] = {0.0}; 
  double alpha_vdim[32] = {0.0}; 

  alpha_cdim[0] = 2.0*p0_over_gamma[0]*dx10; 
  alpha_cdim[3] = 2.0*p0_over_gamma[1]*dx10; 
  alpha_cdim[4] = 2.0*p0_over_gamma[2]*dx10; 
  alpha_cdim[10] = 2.0*p0_over_gamma[3]*dx10; 
  cflFreq_mid += 3.0*fabs(0.125*alpha_cdim[0]); 

  alpha_cdim[16] = 2.0*p1_over_gamma[0]*dx11; 
  alpha_cdim[19] = 2.0*p1_over_gamma[1]*dx11; 
  alpha_cdim[20] = 2.0*p1_over_gamma[2]*dx11; 
  alpha_cdim[26] = 2.0*p1_over_gamma[3]*dx11; 
  cflFreq_mid += 3.0*fabs(0.125*alpha_cdim[16]); 

  alpha_vdim[0] = (B2[0]*p1_over_gamma[0]+2.0*E0[0])*dv10; 
  alpha_vdim[1] = (2.0*E0[1]+p1_over_gamma[0]*B2[1])*dv10; 
  alpha_vdim[2] = (2.0*E0[2]+p1_over_gamma[0]*B2[2])*dv10; 
  alpha_vdim[3] = B2[0]*p1_over_gamma[1]*dv10; 
  alpha_vdim[4] = B2[0]*p1_over_gamma[2]*dv10; 
  alpha_vdim[5] = (2.0*E0[3]+p1_over_gamma[0]*B2[3])*dv10; 
  alpha_vdim[6] = B2[1]*p1_over_gamma[1]*dv10; 
  alpha_vdim[7] = p1_over_gamma[1]*B2[2]*dv10; 
  alpha_vdim[8] = B2[1]*p1_over_gamma[2]*dv10; 
  alpha_vdim[9] = B2[2]*p1_over_gamma[2]*dv10; 
  alpha_vdim[10] = B2[0]*p1_over_gamma[3]*dv10; 
  alpha_vdim[11] = p1_over_gamma[1]*B2[3]*dv10; 
  alpha_vdim[12] = p1_over_gamma[2]*B2[3]*dv10; 
  alpha_vdim[13] = B2[1]*p1_over_gamma[3]*dv10; 
  alpha_vdim[14] = B2[2]*p1_over_gamma[3]*dv10; 
  alpha_vdim[15] = B2[3]*p1_over_gamma[3]*dv10; 
  cflFreq_mid += 3.0*fabs(0.125*alpha_vdim[0]); 

  alpha_vdim[16] = (2.0*E1[0]-1.0*B2[0]*p0_over_gamma[0])*dv11; 
  alpha_vdim[17] = (2.0*E1[1]-1.0*p0_over_gamma[0]*B2[1])*dv11; 
  alpha_vdim[18] = (2.0*E1[2]-1.0*p0_over_gamma[0]*B2[2])*dv11; 
  alpha_vdim[19] = -1.0*B2[0]*p0_over_gamma[1]*dv11; 
  alpha_vdim[20] = -1.0*B2[0]*p0_over_gamma[2]*dv11; 
  alpha_vdim[21] = (2.0*E1[3]-1.0*p0_over_gamma[0]*B2[3])*dv11; 
  alpha_vdim[22] = -1.0*B2[1]*p0_over_gamma[1]*dv11; 
  alpha_vdim[23] = -1.0*p0_over_gamma[1]*B2[2]*dv11; 
  alpha_vdim[24] = -1.0*B2[1]*p0_over_gamma[2]*dv11; 
  alpha_vdim[25] = -1.0*B2[2]*p0_over_gamma[2]*dv11; 
  alpha_vdim[26] = -1.0*B2[0]*p0_over_gamma[3]*dv11; 
  alpha_vdim[27] = -1.0*p0_over_gamma[1]*B2[3]*dv11; 
  alpha_vdim[28] = -1.0*p0_over_gamma[2]*B2[3]*dv11; 
  alpha_vdim[29] = -1.0*B2[1]*p0_over_gamma[3]*dv11; 
  alpha_vdim[30] = -1.0*B2[2]*p0_over_gamma[3]*dv11; 
  alpha_vdim[31] = -1.0*B2[3]*p0_over_gamma[3]*dv11; 
  cflFreq_mid += 3.0*fabs(0.125*alpha_vdim[16]); 

  out[1] += 0.4330127018922193*(alpha_cdim[10]*f[10]+alpha_cdim[4]*f[4]+alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(f[10]*alpha_cdim[26]+f[4]*alpha_cdim[20]+f[3]*alpha_cdim[19]+f[0]*alpha_cdim[16]); 
  out[3] += 0.4330127018922193*(alpha_vdim[15]*f[15]+alpha_vdim[14]*f[14]+alpha_vdim[13]*f[13]+alpha_vdim[12]*f[12]+alpha_vdim[11]*f[11]+alpha_vdim[10]*f[10]+alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[15]*alpha_vdim[31]+f[14]*alpha_vdim[30]+f[13]*alpha_vdim[29]+f[12]*alpha_vdim[28]+f[11]*alpha_vdim[27]+f[10]*alpha_vdim[26]+f[9]*alpha_vdim[25]+f[8]*alpha_vdim[24]+f[7]*alpha_vdim[23]+f[6]*alpha_vdim[22]+f[5]*alpha_vdim[21]+f[4]*alpha_vdim[20]+f[3]*alpha_vdim[19]+f[2]*alpha_vdim[18]+f[1]*alpha_vdim[17]+f[0]*alpha_vdim[16]); 
  out[5] += 0.4330127018922193*(f[13]*alpha_cdim[26]+f[8]*alpha_cdim[20]+f[6]*alpha_cdim[19]+f[1]*alpha_cdim[16]+alpha_cdim[10]*f[14]+alpha_cdim[4]*f[9]+alpha_cdim[3]*f[7]+alpha_cdim[0]*f[2]); 
  out[6] += 0.4330127018922193*(alpha_vdim[14]*f[15]+f[14]*alpha_vdim[15]+alpha_vdim[10]*f[13]+f[10]*alpha_vdim[13]+alpha_vdim[9]*f[12]+f[9]*alpha_vdim[12]+alpha_vdim[7]*f[11]+f[7]*alpha_vdim[11]+alpha_cdim[4]*f[10]+f[4]*alpha_cdim[10]+alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[7] += 0.4330127018922193*(f[4]*alpha_cdim[26]+f[10]*alpha_cdim[20]+f[0]*alpha_cdim[19]+f[3]*alpha_cdim[16]+alpha_vdim[13]*f[15]+f[13]*alpha_vdim[15]+alpha_vdim[10]*f[14]+f[10]*alpha_vdim[14]+alpha_vdim[8]*f[12]+f[8]*alpha_vdim[12]+alpha_vdim[6]*f[11]+f[6]*alpha_vdim[11]+alpha_vdim[4]*f[9]+f[4]*alpha_vdim[9]+alpha_vdim[3]*f[7]+f[3]*alpha_vdim[7]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[8] += 0.4330127018922193*(f[14]*alpha_vdim[31]+f[15]*alpha_vdim[30]+f[10]*alpha_vdim[29]+f[9]*alpha_vdim[28]+f[7]*alpha_vdim[27]+f[13]*alpha_vdim[26]+f[12]*alpha_vdim[25]+f[4]*alpha_vdim[24]+f[11]*alpha_vdim[23]+f[3]*alpha_vdim[22]+f[2]*alpha_vdim[21]+f[8]*alpha_vdim[20]+f[6]*alpha_vdim[19]+f[5]*alpha_vdim[18]+f[0]*alpha_vdim[17]+f[1]*alpha_vdim[16]+alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]); 
  out[9] += 0.4330127018922193*(f[13]*alpha_vdim[31]+f[10]*alpha_vdim[30]+f[15]*alpha_vdim[29]+f[8]*alpha_vdim[28]+f[6]*alpha_vdim[27]+f[14]*alpha_vdim[26]+f[3]*alpha_cdim[26]+f[4]*alpha_vdim[25]+f[12]*alpha_vdim[24]+f[3]*alpha_vdim[23]+f[11]*alpha_vdim[22]+f[1]*alpha_vdim[21]+f[9]*alpha_vdim[20]+f[0]*alpha_cdim[20]+f[7]*alpha_vdim[19]+f[10]*alpha_cdim[19]+f[0]*alpha_vdim[18]+f[5]*alpha_vdim[17]+f[2]*alpha_vdim[16]+f[4]*alpha_cdim[16]); 
  out[10] += 0.4330127018922193*(f[12]*alpha_vdim[31]+f[9]*alpha_vdim[30]+f[8]*alpha_vdim[29]+f[15]*alpha_vdim[28]+f[5]*alpha_vdim[27]+f[4]*alpha_vdim[26]+f[14]*alpha_vdim[25]+f[13]*alpha_vdim[24]+f[2]*alpha_vdim[23]+f[1]*alpha_vdim[22]+f[11]*alpha_vdim[21]+f[10]*alpha_vdim[20]+f[0]*alpha_vdim[19]+f[7]*alpha_vdim[18]+f[6]*alpha_vdim[17]+f[3]*alpha_vdim[16]+alpha_vdim[11]*f[15]+f[11]*alpha_vdim[15]+alpha_vdim[7]*f[14]+f[7]*alpha_vdim[14]+alpha_vdim[6]*f[13]+f[6]*alpha_vdim[13]+alpha_vdim[5]*f[12]+f[5]*alpha_vdim[12]+alpha_vdim[3]*f[10]+f[3]*alpha_vdim[10]+alpha_vdim[2]*f[9]+f[2]*alpha_vdim[9]+alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]); 
  out[11] += 0.4330127018922193*(f[8]*alpha_cdim[26]+f[13]*alpha_cdim[20]+f[1]*alpha_cdim[19]+f[6]*alpha_cdim[16]+alpha_vdim[10]*f[15]+f[10]*alpha_vdim[15]+(alpha_vdim[13]+alpha_cdim[4])*f[14]+f[13]*alpha_vdim[14]+alpha_vdim[4]*f[12]+f[4]*alpha_vdim[12]+alpha_vdim[3]*f[11]+f[3]*alpha_vdim[11]+f[9]*(alpha_cdim[10]+alpha_vdim[8])+f[8]*alpha_vdim[9]+(alpha_vdim[6]+alpha_cdim[0])*f[7]+f[6]*alpha_vdim[7]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+f[2]*(alpha_cdim[3]+alpha_vdim[1])+f[1]*alpha_vdim[2]); 
  out[12] += 0.4330127018922193*(f[10]*alpha_vdim[31]+f[13]*alpha_vdim[30]+f[14]*alpha_vdim[29]+f[4]*alpha_vdim[28]+f[3]*alpha_vdim[27]+f[15]*alpha_vdim[26]+f[6]*alpha_cdim[26]+f[8]*alpha_vdim[25]+f[9]*alpha_vdim[24]+f[6]*alpha_vdim[23]+f[7]*alpha_vdim[22]+f[0]*alpha_vdim[21]+f[12]*alpha_vdim[20]+f[1]*alpha_cdim[20]+f[11]*alpha_vdim[19]+f[13]*alpha_cdim[19]+f[1]*alpha_vdim[18]+f[2]*alpha_vdim[17]+f[5]*alpha_vdim[16]+f[8]*alpha_cdim[16]+alpha_cdim[3]*f[14]+f[7]*alpha_cdim[10]+alpha_cdim[0]*f[9]+f[2]*alpha_cdim[4]); 
  out[13] += 0.4330127018922193*(f[9]*alpha_vdim[31]+f[12]*alpha_vdim[30]+f[4]*alpha_vdim[29]+f[14]*alpha_vdim[28]+f[2]*alpha_vdim[27]+f[8]*alpha_vdim[26]+f[15]*alpha_vdim[25]+f[10]*alpha_vdim[24]+f[5]*alpha_vdim[23]+f[0]*alpha_vdim[22]+f[7]*alpha_vdim[21]+f[13]*alpha_vdim[20]+f[1]*alpha_vdim[19]+f[11]*alpha_vdim[18]+f[3]*alpha_vdim[17]+f[6]*alpha_vdim[16]+alpha_vdim[7]*f[15]+f[7]*alpha_vdim[15]+alpha_vdim[11]*f[14]+f[11]*alpha_vdim[14]+alpha_vdim[3]*f[13]+f[3]*alpha_vdim[13]+alpha_vdim[2]*f[12]+f[2]*alpha_vdim[12]+(alpha_vdim[6]+alpha_cdim[0])*f[10]+f[6]*alpha_vdim[10]+f[0]*alpha_cdim[10]+alpha_vdim[5]*f[9]+f[5]*alpha_vdim[9]+alpha_vdim[0]*f[8]+f[0]*alpha_vdim[8]+(alpha_cdim[3]+alpha_vdim[1])*f[4]+f[1]*alpha_vdim[4]+f[3]*alpha_cdim[4]); 
  out[14] += 0.4330127018922193*(f[8]*alpha_vdim[31]+f[4]*alpha_vdim[30]+f[12]*alpha_vdim[29]+f[13]*alpha_vdim[28]+f[1]*alpha_vdim[27]+f[9]*alpha_vdim[26]+f[0]*alpha_cdim[26]+f[10]*alpha_vdim[25]+f[15]*alpha_vdim[24]+f[0]*alpha_vdim[23]+f[5]*alpha_vdim[22]+f[6]*alpha_vdim[21]+f[14]*alpha_vdim[20]+f[3]*alpha_cdim[20]+f[2]*alpha_vdim[19]+f[4]*alpha_cdim[19]+f[3]*alpha_vdim[18]+f[11]*alpha_vdim[17]+f[7]*alpha_vdim[16]+f[10]*alpha_cdim[16]+alpha_vdim[6]*f[15]+f[6]*alpha_vdim[15]+alpha_vdim[3]*f[14]+f[3]*alpha_vdim[14]+alpha_vdim[11]*f[13]+f[11]*alpha_vdim[13]+alpha_vdim[1]*f[12]+f[1]*alpha_vdim[12]+alpha_vdim[7]*f[10]+f[7]*alpha_vdim[10]+alpha_vdim[0]*f[9]+f[0]*alpha_vdim[9]+alpha_vdim[5]*f[8]+f[5]*alpha_vdim[8]+alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]); 
  out[15] += 0.4330127018922193*(f[4]*alpha_vdim[31]+f[8]*alpha_vdim[30]+f[9]*alpha_vdim[29]+f[10]*alpha_vdim[28]+f[0]*alpha_vdim[27]+f[12]*alpha_vdim[26]+f[1]*alpha_cdim[26]+f[13]*alpha_vdim[25]+f[14]*alpha_vdim[24]+f[1]*alpha_vdim[23]+f[2]*alpha_vdim[22]+f[3]*alpha_vdim[21]+f[15]*alpha_vdim[20]+f[6]*alpha_cdim[20]+f[5]*alpha_vdim[19]+f[8]*alpha_cdim[19]+f[6]*alpha_vdim[18]+f[7]*alpha_vdim[17]+f[11]*alpha_vdim[16]+f[13]*alpha_cdim[16]+alpha_vdim[3]*f[15]+f[3]*alpha_vdim[15]+(alpha_vdim[6]+alpha_cdim[0])*f[14]+f[6]*alpha_vdim[14]+alpha_vdim[7]*f[13]+f[7]*alpha_vdim[13]+alpha_vdim[0]*f[12]+f[0]*alpha_vdim[12]+alpha_vdim[10]*f[11]+f[10]*alpha_vdim[11]+f[2]*alpha_cdim[10]+(alpha_cdim[3]+alpha_vdim[1])*f[9]+f[1]*alpha_vdim[9]+alpha_vdim[2]*f[8]+f[2]*alpha_vdim[8]+alpha_cdim[4]*f[7]+alpha_vdim[4]*f[5]+f[4]*alpha_vdim[5]); 

  return cflFreq_mid; 
} 
