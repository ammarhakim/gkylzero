#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // qmem:      q/m*EM fields.
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &qmem[0]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &qmem[2]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double dv12 = 2/dxv[3]; 
  const double *E2 = &qmem[4]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 

  const double *B0 = &qmem[6]; 
  const double *B1 = &qmem[8]; 
  const double *B2 = &qmem[10]; 
  double cflFreq_mid = 0.0; 
  double alpha_cdim[16] = {0.0}; 
  double alpha_vdim[48] = {0.0}; 

  alpha_cdim[0] = 1.414213562373095*p0_over_gamma[0]*dx10; 
  alpha_cdim[2] = 1.414213562373095*p0_over_gamma[1]*dx10; 
  alpha_cdim[3] = 1.414213562373095*p0_over_gamma[2]*dx10; 
  alpha_cdim[4] = 1.414213562373095*p0_over_gamma[3]*dx10; 
  alpha_cdim[7] = 1.414213562373095*p0_over_gamma[4]*dx10; 
  alpha_cdim[9] = 1.414213562373095*p0_over_gamma[5]*dx10; 
  alpha_cdim[10] = 1.414213562373095*p0_over_gamma[6]*dx10; 
  alpha_cdim[14] = 1.414213562373095*p0_over_gamma[7]*dx10; 
  cflFreq_mid += 3.0*fabs(0.125*alpha_cdim[0]); 

  alpha_vdim[0] = ((-1.0*B1[0]*p2_over_gamma[0])+B2[0]*p1_over_gamma[0]+2.828427124746191*E0[0])*dv10; 
  alpha_vdim[1] = (2.828427124746191*E0[1]+p1_over_gamma[0]*B2[1]-1.0*p2_over_gamma[0]*B1[1])*dv10; 
  alpha_vdim[2] = (B2[0]*p1_over_gamma[1]-1.0*B1[0]*p2_over_gamma[1])*dv10; 
  alpha_vdim[3] = (B2[0]*p1_over_gamma[2]-1.0*B1[0]*p2_over_gamma[2])*dv10; 
  alpha_vdim[4] = (B2[0]*p1_over_gamma[3]-1.0*B1[0]*p2_over_gamma[3])*dv10; 
  alpha_vdim[5] = (B2[1]*p1_over_gamma[1]-1.0*B1[1]*p2_over_gamma[1])*dv10; 
  alpha_vdim[6] = (B2[1]*p1_over_gamma[2]-1.0*B1[1]*p2_over_gamma[2])*dv10; 
  alpha_vdim[7] = (B2[0]*p1_over_gamma[4]-1.0*B1[0]*p2_over_gamma[4])*dv10; 
  alpha_vdim[8] = (B2[1]*p1_over_gamma[3]-1.0*B1[1]*p2_over_gamma[3])*dv10; 
  alpha_vdim[9] = (B2[0]*p1_over_gamma[5]-1.0*B1[0]*p2_over_gamma[5])*dv10; 
  alpha_vdim[10] = (B2[0]*p1_over_gamma[6]-1.0*B1[0]*p2_over_gamma[6])*dv10; 
  alpha_vdim[11] = (B2[1]*p1_over_gamma[4]-1.0*B1[1]*p2_over_gamma[4])*dv10; 
  alpha_vdim[12] = (B2[1]*p1_over_gamma[5]-1.0*B1[1]*p2_over_gamma[5])*dv10; 
  alpha_vdim[13] = (B2[1]*p1_over_gamma[6]-1.0*B1[1]*p2_over_gamma[6])*dv10; 
  alpha_vdim[14] = (B2[0]*p1_over_gamma[7]-1.0*B1[0]*p2_over_gamma[7])*dv10; 
  alpha_vdim[15] = (B2[1]*p1_over_gamma[7]-1.0*B1[1]*p2_over_gamma[7])*dv10; 
  cflFreq_mid += 3.0*fabs(0.125*alpha_vdim[0]); 

  alpha_vdim[16] = (B0[0]*p2_over_gamma[0]-1.0*B2[0]*p0_over_gamma[0]+2.828427124746191*E1[0])*dv11; 
  alpha_vdim[17] = (2.828427124746191*E1[1]-1.0*p0_over_gamma[0]*B2[1]+p2_over_gamma[0]*B0[1])*dv11; 
  alpha_vdim[18] = (B0[0]*p2_over_gamma[1]-1.0*B2[0]*p0_over_gamma[1])*dv11; 
  alpha_vdim[19] = (B0[0]*p2_over_gamma[2]-1.0*B2[0]*p0_over_gamma[2])*dv11; 
  alpha_vdim[20] = (B0[0]*p2_over_gamma[3]-1.0*B2[0]*p0_over_gamma[3])*dv11; 
  alpha_vdim[21] = (B0[1]*p2_over_gamma[1]-1.0*B2[1]*p0_over_gamma[1])*dv11; 
  alpha_vdim[22] = (B0[1]*p2_over_gamma[2]-1.0*B2[1]*p0_over_gamma[2])*dv11; 
  alpha_vdim[23] = (B0[0]*p2_over_gamma[4]-1.0*B2[0]*p0_over_gamma[4])*dv11; 
  alpha_vdim[24] = (B0[1]*p2_over_gamma[3]-1.0*B2[1]*p0_over_gamma[3])*dv11; 
  alpha_vdim[25] = (B0[0]*p2_over_gamma[5]-1.0*B2[0]*p0_over_gamma[5])*dv11; 
  alpha_vdim[26] = (B0[0]*p2_over_gamma[6]-1.0*B2[0]*p0_over_gamma[6])*dv11; 
  alpha_vdim[27] = (B0[1]*p2_over_gamma[4]-1.0*B2[1]*p0_over_gamma[4])*dv11; 
  alpha_vdim[28] = (B0[1]*p2_over_gamma[5]-1.0*B2[1]*p0_over_gamma[5])*dv11; 
  alpha_vdim[29] = (B0[1]*p2_over_gamma[6]-1.0*B2[1]*p0_over_gamma[6])*dv11; 
  alpha_vdim[30] = (B0[0]*p2_over_gamma[7]-1.0*B2[0]*p0_over_gamma[7])*dv11; 
  alpha_vdim[31] = (B0[1]*p2_over_gamma[7]-1.0*B2[1]*p0_over_gamma[7])*dv11; 
  cflFreq_mid += 3.0*fabs(0.125*alpha_vdim[16]); 

  alpha_vdim[32] = ((-1.0*B0[0]*p1_over_gamma[0])+B1[0]*p0_over_gamma[0]+2.828427124746191*E2[0])*dv12; 
  alpha_vdim[33] = (2.828427124746191*E2[1]+p0_over_gamma[0]*B1[1]-1.0*p1_over_gamma[0]*B0[1])*dv12; 
  alpha_vdim[34] = (B1[0]*p0_over_gamma[1]-1.0*B0[0]*p1_over_gamma[1])*dv12; 
  alpha_vdim[35] = (B1[0]*p0_over_gamma[2]-1.0*B0[0]*p1_over_gamma[2])*dv12; 
  alpha_vdim[36] = (B1[0]*p0_over_gamma[3]-1.0*B0[0]*p1_over_gamma[3])*dv12; 
  alpha_vdim[37] = (B1[1]*p0_over_gamma[1]-1.0*B0[1]*p1_over_gamma[1])*dv12; 
  alpha_vdim[38] = (B1[1]*p0_over_gamma[2]-1.0*B0[1]*p1_over_gamma[2])*dv12; 
  alpha_vdim[39] = (B1[0]*p0_over_gamma[4]-1.0*B0[0]*p1_over_gamma[4])*dv12; 
  alpha_vdim[40] = (B1[1]*p0_over_gamma[3]-1.0*B0[1]*p1_over_gamma[3])*dv12; 
  alpha_vdim[41] = (B1[0]*p0_over_gamma[5]-1.0*B0[0]*p1_over_gamma[5])*dv12; 
  alpha_vdim[42] = (B1[0]*p0_over_gamma[6]-1.0*B0[0]*p1_over_gamma[6])*dv12; 
  alpha_vdim[43] = (B1[1]*p0_over_gamma[4]-1.0*B0[1]*p1_over_gamma[4])*dv12; 
  alpha_vdim[44] = (B1[1]*p0_over_gamma[5]-1.0*B0[1]*p1_over_gamma[5])*dv12; 
  alpha_vdim[45] = (B1[1]*p0_over_gamma[6]-1.0*B0[1]*p1_over_gamma[6])*dv12; 
  alpha_vdim[46] = (B1[0]*p0_over_gamma[7]-1.0*B0[0]*p1_over_gamma[7])*dv12; 
  alpha_vdim[47] = (B1[1]*p0_over_gamma[7]-1.0*B0[1]*p1_over_gamma[7])*dv12; 
  cflFreq_mid += 3.0*fabs(0.125*alpha_vdim[32]); 

  out[1] += 0.4330127018922193*(alpha_cdim[14]*f[14]+alpha_cdim[10]*f[10]+alpha_cdim[9]*f[9]+alpha_cdim[7]*f[7]+alpha_cdim[4]*f[4]+alpha_cdim[3]*f[3]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha_vdim[15]*f[15]+alpha_vdim[14]*f[14]+alpha_vdim[13]*f[13]+alpha_vdim[12]*f[12]+alpha_vdim[11]*f[11]+alpha_vdim[10]*f[10]+alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[15]*alpha_vdim[31]+f[14]*alpha_vdim[30]+f[13]*alpha_vdim[29]+f[12]*alpha_vdim[28]+f[11]*alpha_vdim[27]+f[10]*alpha_vdim[26]+f[9]*alpha_vdim[25]+f[8]*alpha_vdim[24]+f[7]*alpha_vdim[23]+f[6]*alpha_vdim[22]+f[5]*alpha_vdim[21]+f[4]*alpha_vdim[20]+f[3]*alpha_vdim[19]+f[2]*alpha_vdim[18]+f[1]*alpha_vdim[17]+f[0]*alpha_vdim[16]); 
  out[4] += 0.4330127018922193*(f[15]*alpha_vdim[47]+f[14]*alpha_vdim[46]+f[13]*alpha_vdim[45]+f[12]*alpha_vdim[44]+f[11]*alpha_vdim[43]+f[10]*alpha_vdim[42]+f[9]*alpha_vdim[41]+f[8]*alpha_vdim[40]+f[7]*alpha_vdim[39]+f[6]*alpha_vdim[38]+f[5]*alpha_vdim[37]+f[4]*alpha_vdim[36]+f[3]*alpha_vdim[35]+f[2]*alpha_vdim[34]+f[1]*alpha_vdim[33]+f[0]*alpha_vdim[32]); 
  out[5] += 0.4330127018922193*(alpha_vdim[14]*f[15]+f[14]*(alpha_vdim[15]+alpha_cdim[10])+f[10]*alpha_cdim[14]+alpha_vdim[10]*f[13]+f[10]*alpha_vdim[13]+alpha_vdim[9]*f[12]+f[9]*alpha_vdim[12]+alpha_vdim[7]*f[11]+f[7]*alpha_vdim[11]+alpha_cdim[4]*f[9]+f[4]*alpha_cdim[9]+alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]+alpha_cdim[3]*f[7]+f[3]*alpha_cdim[7]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_vdim[2]*f[5]+f[2]*(alpha_vdim[5]+alpha_cdim[0])+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.4330127018922193*(f[14]*alpha_vdim[31]+f[15]*alpha_vdim[30]+f[10]*alpha_vdim[29]+f[9]*alpha_vdim[28]+f[7]*alpha_vdim[27]+f[13]*alpha_vdim[26]+f[12]*alpha_vdim[25]+f[4]*alpha_vdim[24]+f[11]*alpha_vdim[23]+f[3]*alpha_vdim[22]+f[2]*alpha_vdim[21]+f[8]*alpha_vdim[20]+f[6]*alpha_vdim[19]+f[5]*alpha_vdim[18]+f[0]*alpha_vdim[17]+f[1]*alpha_vdim[16]+alpha_cdim[9]*f[14]+f[9]*alpha_cdim[14]+alpha_cdim[4]*f[10]+f[4]*alpha_cdim[10]+alpha_cdim[2]*f[7]+f[2]*alpha_cdim[7]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]); 
  out[7] += 0.4330127018922193*(f[13]*alpha_vdim[31]+f[10]*alpha_vdim[30]+f[15]*alpha_vdim[29]+f[8]*alpha_vdim[28]+f[6]*alpha_vdim[27]+f[14]*alpha_vdim[26]+f[4]*alpha_vdim[25]+f[12]*alpha_vdim[24]+f[3]*alpha_vdim[23]+f[11]*alpha_vdim[22]+f[1]*alpha_vdim[21]+f[9]*alpha_vdim[20]+f[7]*alpha_vdim[19]+f[0]*alpha_vdim[18]+f[5]*alpha_vdim[17]+f[2]*alpha_vdim[16]+alpha_vdim[12]*f[15]+f[12]*alpha_vdim[15]+alpha_vdim[9]*f[14]+f[9]*alpha_vdim[14]+alpha_vdim[8]*f[13]+f[8]*alpha_vdim[13]+alpha_vdim[5]*f[11]+f[5]*alpha_vdim[11]+alpha_vdim[4]*f[10]+f[4]*alpha_vdim[10]+alpha_vdim[2]*f[7]+f[2]*alpha_vdim[7]+alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[8] += 0.4330127018922193*(f[14]*alpha_vdim[47]+f[15]*alpha_vdim[46]+f[10]*alpha_vdim[45]+f[9]*alpha_vdim[44]+f[7]*alpha_vdim[43]+f[13]*alpha_vdim[42]+f[12]*alpha_vdim[41]+f[4]*alpha_vdim[40]+f[11]*alpha_vdim[39]+f[3]*alpha_vdim[38]+f[2]*alpha_vdim[37]+f[8]*alpha_vdim[36]+f[6]*alpha_vdim[35]+f[5]*alpha_vdim[34]+f[0]*alpha_vdim[33]+f[1]*alpha_vdim[32]+alpha_cdim[7]*f[14]+f[7]*alpha_cdim[14]+alpha_cdim[3]*f[10]+f[3]*alpha_cdim[10]+alpha_cdim[2]*f[9]+f[2]*alpha_cdim[9]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]); 
  out[9] += 0.4330127018922193*(f[13]*alpha_vdim[47]+f[10]*alpha_vdim[46]+f[15]*alpha_vdim[45]+f[8]*alpha_vdim[44]+f[6]*alpha_vdim[43]+f[14]*alpha_vdim[42]+f[4]*alpha_vdim[41]+f[12]*alpha_vdim[40]+f[3]*alpha_vdim[39]+f[11]*alpha_vdim[38]+f[1]*alpha_vdim[37]+f[9]*alpha_vdim[36]+f[7]*alpha_vdim[35]+f[0]*alpha_vdim[34]+f[5]*alpha_vdim[33]+f[2]*alpha_vdim[32]+alpha_vdim[11]*f[15]+f[11]*alpha_vdim[15]+alpha_vdim[7]*f[14]+f[7]*alpha_vdim[14]+alpha_vdim[6]*f[13]+f[6]*alpha_vdim[13]+alpha_vdim[5]*f[12]+f[5]*alpha_vdim[12]+alpha_vdim[3]*f[10]+f[3]*alpha_vdim[10]+alpha_vdim[2]*f[9]+f[2]*alpha_vdim[9]+alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]); 
  out[10] += 0.4330127018922193*(f[12]*alpha_vdim[47]+f[9]*alpha_vdim[46]+f[8]*alpha_vdim[45]+f[15]*alpha_vdim[44]+f[5]*alpha_vdim[43]+f[4]*alpha_vdim[42]+f[14]*alpha_vdim[41]+f[13]*alpha_vdim[40]+f[2]*alpha_vdim[39]+f[1]*alpha_vdim[38]+f[11]*alpha_vdim[37]+f[10]*alpha_vdim[36]+f[0]*alpha_vdim[35]+f[7]*alpha_vdim[34]+f[6]*alpha_vdim[33]+f[3]*alpha_vdim[32]+f[11]*alpha_vdim[31]+f[7]*alpha_vdim[30]+f[6]*alpha_vdim[29]+f[5]*alpha_vdim[28]+f[15]*alpha_vdim[27]+f[3]*alpha_vdim[26]+f[2]*alpha_vdim[25]+f[1]*alpha_vdim[24]+f[14]*alpha_vdim[23]+f[13]*alpha_vdim[22]+f[12]*alpha_vdim[21]+f[0]*alpha_vdim[20]+f[10]*alpha_vdim[19]+f[9]*alpha_vdim[18]+f[8]*alpha_vdim[17]+f[4]*alpha_vdim[16]); 
  out[11] += 0.4330127018922193*(f[10]*alpha_vdim[31]+f[13]*alpha_vdim[30]+f[14]*alpha_vdim[29]+f[4]*alpha_vdim[28]+f[3]*alpha_vdim[27]+f[15]*alpha_vdim[26]+f[8]*alpha_vdim[25]+f[9]*alpha_vdim[24]+f[6]*alpha_vdim[23]+f[7]*alpha_vdim[22]+f[0]*alpha_vdim[21]+f[12]*alpha_vdim[20]+f[11]*alpha_vdim[19]+f[1]*alpha_vdim[18]+f[2]*alpha_vdim[17]+f[5]*alpha_vdim[16]+alpha_vdim[9]*f[15]+f[9]*alpha_vdim[15]+(alpha_vdim[12]+alpha_cdim[4])*f[14]+f[12]*alpha_vdim[14]+f[4]*alpha_cdim[14]+alpha_vdim[4]*f[13]+f[4]*alpha_vdim[13]+alpha_vdim[2]*f[11]+f[2]*alpha_vdim[11]+(alpha_cdim[9]+alpha_vdim[8])*f[10]+f[8]*alpha_vdim[10]+f[9]*alpha_cdim[10]+(alpha_vdim[5]+alpha_cdim[0])*f[7]+f[5]*alpha_vdim[7]+f[0]*alpha_cdim[7]+alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]+f[2]*alpha_cdim[3]); 
  out[12] += 0.4330127018922193*(f[10]*alpha_vdim[47]+f[13]*alpha_vdim[46]+f[14]*alpha_vdim[45]+f[4]*alpha_vdim[44]+f[3]*alpha_vdim[43]+f[15]*alpha_vdim[42]+f[8]*alpha_vdim[41]+f[9]*alpha_vdim[40]+f[6]*alpha_vdim[39]+f[7]*alpha_vdim[38]+f[0]*alpha_vdim[37]+f[12]*alpha_vdim[36]+f[11]*alpha_vdim[35]+f[1]*alpha_vdim[34]+f[2]*alpha_vdim[33]+f[5]*alpha_vdim[32]+alpha_vdim[7]*f[15]+f[7]*alpha_vdim[15]+(alpha_vdim[11]+alpha_cdim[3])*f[14]+f[11]*alpha_vdim[14]+f[3]*alpha_cdim[14]+alpha_vdim[3]*f[13]+f[3]*alpha_vdim[13]+alpha_vdim[2]*f[12]+f[2]*alpha_vdim[12]+(alpha_cdim[7]+alpha_vdim[6])*f[10]+f[6]*alpha_vdim[10]+f[7]*alpha_cdim[10]+(alpha_vdim[5]+alpha_cdim[0])*f[9]+f[5]*alpha_vdim[9]+f[0]*alpha_cdim[9]+alpha_vdim[0]*f[8]+f[0]*alpha_vdim[8]+(alpha_cdim[2]+alpha_vdim[1])*f[4]+f[1]*alpha_vdim[4]+f[2]*alpha_cdim[4]); 
  out[13] += 0.4330127018922193*(f[9]*alpha_vdim[47]+f[12]*alpha_vdim[46]+f[4]*alpha_vdim[45]+f[14]*alpha_vdim[44]+f[2]*alpha_vdim[43]+f[8]*alpha_vdim[42]+f[15]*alpha_vdim[41]+f[10]*alpha_vdim[40]+f[5]*alpha_vdim[39]+f[0]*alpha_vdim[38]+f[7]*alpha_vdim[37]+f[13]*alpha_vdim[36]+f[1]*alpha_vdim[35]+f[11]*alpha_vdim[34]+f[3]*alpha_vdim[33]+f[6]*alpha_vdim[32]+f[7]*alpha_vdim[31]+f[11]*alpha_vdim[30]+f[3]*alpha_vdim[29]+f[2]*alpha_vdim[28]+f[14]*alpha_vdim[27]+f[6]*alpha_vdim[26]+f[5]*alpha_vdim[25]+f[0]*alpha_vdim[24]+f[15]*alpha_vdim[23]+f[10]*alpha_vdim[22]+f[9]*alpha_vdim[21]+f[1]*alpha_vdim[20]+f[13]*alpha_vdim[19]+f[12]*alpha_vdim[18]+f[4]*alpha_vdim[17]+f[8]*alpha_vdim[16]+alpha_cdim[2]*f[14]+f[2]*alpha_cdim[14]+alpha_cdim[0]*f[10]+f[0]*alpha_cdim[10]+alpha_cdim[7]*f[9]+f[7]*alpha_cdim[9]+alpha_cdim[3]*f[4]+f[3]*alpha_cdim[4]); 
  out[14] += 0.4330127018922193*(f[8]*alpha_vdim[47]+f[4]*alpha_vdim[46]+f[12]*alpha_vdim[45]+f[13]*alpha_vdim[44]+f[1]*alpha_vdim[43]+f[9]*alpha_vdim[42]+f[10]*alpha_vdim[41]+f[15]*alpha_vdim[40]+f[0]*alpha_vdim[39]+f[5]*alpha_vdim[38]+f[6]*alpha_vdim[37]+f[14]*alpha_vdim[36]+f[2]*alpha_vdim[35]+f[3]*alpha_vdim[34]+f[11]*alpha_vdim[33]+f[7]*alpha_vdim[32]+f[6]*alpha_vdim[31]+f[3]*alpha_vdim[30]+f[11]*alpha_vdim[29]+f[1]*alpha_vdim[28]+f[13]*alpha_vdim[27]+f[7]*alpha_vdim[26]+f[0]*alpha_vdim[25]+f[5]*alpha_vdim[24]+f[10]*alpha_vdim[23]+f[15]*alpha_vdim[22]+f[8]*alpha_vdim[21]+f[2]*alpha_vdim[20]+f[14]*alpha_vdim[19]+f[4]*alpha_vdim[18]+f[12]*alpha_vdim[17]+f[9]*alpha_vdim[16]+alpha_vdim[5]*f[15]+f[5]*alpha_vdim[15]+alpha_vdim[2]*f[14]+f[2]*alpha_vdim[14]+alpha_vdim[1]*f[13]+f[1]*alpha_vdim[13]+alpha_vdim[11]*f[12]+f[11]*alpha_vdim[12]+alpha_vdim[0]*f[10]+f[0]*alpha_vdim[10]+alpha_vdim[7]*f[9]+f[7]*alpha_vdim[9]+alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8]+alpha_vdim[3]*f[4]+f[3]*alpha_vdim[4]); 
  out[15] += 0.4330127018922193*(f[4]*alpha_vdim[47]+f[8]*alpha_vdim[46]+f[9]*alpha_vdim[45]+f[10]*alpha_vdim[44]+f[0]*alpha_vdim[43]+f[12]*alpha_vdim[42]+f[13]*alpha_vdim[41]+f[14]*alpha_vdim[40]+f[1]*alpha_vdim[39]+f[2]*alpha_vdim[38]+f[3]*alpha_vdim[37]+f[15]*alpha_vdim[36]+f[5]*alpha_vdim[35]+f[6]*alpha_vdim[34]+f[7]*alpha_vdim[33]+f[11]*alpha_vdim[32]+f[3]*alpha_vdim[31]+f[6]*alpha_vdim[30]+f[7]*alpha_vdim[29]+f[0]*alpha_vdim[28]+f[10]*alpha_vdim[27]+f[11]*alpha_vdim[26]+f[1]*alpha_vdim[25]+f[2]*alpha_vdim[24]+f[13]*alpha_vdim[23]+f[14]*alpha_vdim[22]+f[4]*alpha_vdim[21]+f[5]*alpha_vdim[20]+f[15]*alpha_vdim[19]+f[8]*alpha_vdim[18]+f[9]*alpha_vdim[17]+f[12]*alpha_vdim[16]+alpha_vdim[2]*f[15]+f[2]*alpha_vdim[15]+(alpha_vdim[5]+alpha_cdim[0])*f[14]+f[5]*alpha_vdim[14]+f[0]*alpha_cdim[14]+alpha_vdim[0]*f[13]+f[0]*alpha_vdim[13]+alpha_vdim[7]*f[12]+f[7]*alpha_vdim[12]+alpha_vdim[9]*f[11]+f[9]*alpha_vdim[11]+(alpha_cdim[2]+alpha_vdim[1])*f[10]+f[1]*alpha_vdim[10]+f[2]*alpha_cdim[10]+alpha_cdim[3]*f[9]+f[3]*alpha_cdim[9]+alpha_vdim[3]*f[8]+f[3]*alpha_vdim[8]+alpha_cdim[4]*f[7]+f[4]*alpha_cdim[7]+alpha_vdim[4]*f[6]+f[4]*alpha_vdim[6]); 

  return cflFreq_mid; 
} 
