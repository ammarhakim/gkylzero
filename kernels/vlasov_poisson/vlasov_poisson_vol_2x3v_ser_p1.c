#include <gkyl_vlasov_poisson_kernels.h> 

GKYL_CU_DH double vlasov_poisson_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:     potential (scaled by appropriate factors).
  // f:         Input distribution function.
  // out:       Incremented output.

  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv11 = 2/dxv[3]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv12 = 2/dxv[4]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  const double dx11 = 2/dxv[1]; 

  const double *phi = &field[0]; 

  double cflFreq_mid = 0.0; 

  double alpha_cdim[160] = {0.0}; 
  double alpha_vdim[240] = {0.0}; 

  alpha_cdim[0] = 11.31370849898477*w0dx0; 
  alpha_cdim[3] = 3.265986323710906*dv0dx0; 

  cflFreq_mid += 3.0*(fabs(w0dx0)+0.5*dv0dx0); 

  alpha_cdim[80] = 11.31370849898477*w1dx1; 
  alpha_cdim[84] = 3.265986323710906*dv1dx1; 

  cflFreq_mid += 3.0*(fabs(w1dx1)+0.5*dv1dx1); 

  alpha_vdim[0] = -4.898979485566357*phi[1]*dv10*dx10; 
  alpha_vdim[2] = -4.898979485566357*phi[3]*dv10*dx10; 

  cflFreq_mid += 5.0*fabs(0.0883883476483184*alpha_vdim[0]); 

  alpha_vdim[80] = -4.898979485566357*phi[2]*dv11*dx11; 
  alpha_vdim[81] = -4.898979485566357*phi[3]*dv11*dx11; 

  cflFreq_mid += 5.0*fabs(0.0883883476483184*alpha_vdim[0]); 


  out[1] += 0.3061862178478971*(alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.3061862178478971*(alpha_cdim[4]*f[4]+alpha_cdim[0]*f[0]); 
  out[3] += 0.3061862178478971*(alpha_vdim[2]*f[2]+alpha_vdim[0]*f[0]); 
  out[4] += 0.3061862178478971*(alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[6] += 0.3061862178478971*(alpha_cdim[4]*f[9]+alpha_cdim[3]*f[8]+alpha_cdim[0]*(f[2]+f[1])); 
  out[7] += 0.273861278752583*alpha_cdim[3]*f[32]+0.3061862178478971*(alpha_vdim[2]*f[6]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_vdim[0]*f[1]); 
  out[8] += 0.3061862178478971*(alpha_cdim[4]*f[11]+alpha_cdim[0]*f[3]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[9] += 0.3061862178478971*(alpha_cdim[3]*f[11]+alpha_cdim[0]*f[4]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[10] += 0.273861278752583*alpha_cdim[4]*f[48]+0.3061862178478971*(alpha_vdim[1]*f[6]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]+alpha_vdim[0]*f[2]); 
  out[11] += 0.3061862178478971*(alpha_vdim[2]*f[10]+alpha_vdim[1]*f[7]+alpha_vdim[0]*(f[4]+f[3])); 
  out[12] += 0.3061862178478971*(alpha_cdim[3]*f[14]+alpha_cdim[0]*f[5]); 
  out[13] += 0.3061862178478971*(alpha_cdim[4]*f[15]+alpha_cdim[0]*f[5]); 
  out[14] += 0.3061862178478971*(alpha_vdim[2]*f[13]+alpha_vdim[0]*f[5]); 
  out[15] += 0.3061862178478971*(alpha_vdim[1]*f[12]+alpha_vdim[0]*f[5]); 
  out[16] += 0.273861278752583*alpha_cdim[3]*f[34]+0.3061862178478971*(alpha_cdim[4]*f[18]+alpha_cdim[0]*(f[8]+f[7])+alpha_vdim[0]*f[6]+f[2]*alpha_cdim[3]+f[1]*alpha_vdim[2]); 
  out[17] += 0.273861278752583*alpha_cdim[4]*f[49]+0.3061862178478971*(alpha_cdim[3]*f[19]+alpha_cdim[0]*(f[10]+f[9])+alpha_vdim[0]*f[6]+f[1]*alpha_cdim[4]+alpha_vdim[1]*f[2]); 
  out[18] += 0.273861278752583*alpha_cdim[3]*f[35]+0.3061862178478971*(alpha_vdim[2]*f[17]+alpha_cdim[0]*f[11]+alpha_vdim[0]*(f[9]+f[7])+alpha_cdim[3]*f[4]+alpha_vdim[1]*f[3]); 
  out[19] += 0.273861278752583*alpha_cdim[4]*f[51]+0.3061862178478971*(alpha_vdim[1]*f[16]+alpha_cdim[0]*f[11]+alpha_vdim[0]*(f[10]+f[8])+alpha_vdim[2]*f[4]+f[3]*alpha_cdim[4]); 
  out[20] += 0.3061862178478971*(alpha_cdim[4]*f[23]+alpha_cdim[3]*f[22]+alpha_cdim[0]*(f[13]+f[12])); 
  out[21] += 0.273861278752583*alpha_cdim[3]*f[36]+0.3061862178478971*(alpha_vdim[2]*f[20]+alpha_cdim[0]*f[14]+alpha_vdim[0]*f[12]+alpha_cdim[3]*f[5]); 
  out[22] += 0.3061862178478971*(alpha_cdim[4]*f[25]+alpha_cdim[0]*f[14]+alpha_vdim[0]*f[13]+alpha_vdim[2]*f[5]); 
  out[23] += 0.3061862178478971*(alpha_cdim[3]*f[25]+alpha_cdim[0]*f[15]+alpha_vdim[0]*f[12]+alpha_vdim[1]*f[5]); 
  out[24] += 0.273861278752583*alpha_cdim[4]*f[52]+0.3061862178478971*(alpha_vdim[1]*f[20]+alpha_cdim[0]*f[15]+alpha_vdim[0]*f[13]+alpha_cdim[4]*f[5]); 
  out[25] += 0.3061862178478971*(alpha_vdim[2]*f[24]+alpha_vdim[1]*f[21]+alpha_vdim[0]*(f[15]+f[14])); 
  out[26] += 0.273861278752583*(alpha_cdim[4]*f[54]+alpha_cdim[3]*f[39])+0.3061862178478971*(alpha_cdim[0]*(f[19]+f[18])+alpha_vdim[0]*(f[17]+f[16])+alpha_cdim[3]*f[10]+alpha_vdim[2]*f[9]+alpha_vdim[1]*f[8]+alpha_cdim[4]*f[7]); 
  out[27] += 0.273861278752583*alpha_cdim[3]*f[41]+0.3061862178478971*(alpha_cdim[4]*f[29]+alpha_cdim[0]*(f[22]+f[21])+alpha_vdim[0]*f[20]+alpha_cdim[3]*f[13]+alpha_vdim[2]*f[12]); 
  out[28] += 0.273861278752583*alpha_cdim[4]*f[56]+0.3061862178478971*(alpha_cdim[3]*f[30]+alpha_cdim[0]*(f[24]+f[23])+alpha_vdim[0]*f[20]+alpha_vdim[1]*f[13]+alpha_cdim[4]*f[12]); 
  out[29] += 0.273861278752583*alpha_cdim[3]*f[42]+0.3061862178478971*(alpha_vdim[2]*f[28]+alpha_cdim[0]*f[25]+alpha_vdim[0]*(f[23]+f[21])+alpha_cdim[3]*f[15]+alpha_vdim[1]*f[14]); 
  out[30] += 0.273861278752583*alpha_cdim[4]*f[58]+0.3061862178478971*(alpha_vdim[1]*f[27]+alpha_cdim[0]*f[25]+alpha_vdim[0]*(f[24]+f[22])+alpha_vdim[2]*f[15]+alpha_cdim[4]*f[14]); 
  out[31] += 0.273861278752583*(alpha_cdim[4]*f[61]+alpha_cdim[3]*f[46])+0.3061862178478971*(alpha_cdim[0]*(f[30]+f[29])+alpha_vdim[0]*(f[28]+f[27])+alpha_cdim[3]*f[24]+alpha_vdim[2]*f[23]+alpha_vdim[1]*f[22]+alpha_cdim[4]*f[21]); 
  out[32] += 0.6846531968814573*(alpha_vdim[2]*f[8]+alpha_vdim[0]*f[3]); 
  out[33] += 0.3061862178478971*alpha_cdim[0]*f[32]+0.6846531968814573*(alpha_vdim[2]*f[16]+alpha_vdim[0]*f[7])+0.273861278752583*alpha_cdim[3]*f[3]; 
  out[34] += 0.3061862178478971*(alpha_cdim[4]*f[35]+alpha_cdim[0]*f[32])+0.6846531968814573*(alpha_vdim[0]*f[8]+alpha_vdim[2]*f[3]); 
  out[35] += 0.3061862178478971*(alpha_vdim[1]*f[33]+alpha_vdim[0]*f[32])+0.6846531968814573*(alpha_vdim[2]*f[19]+alpha_vdim[0]*f[11]); 
  out[36] += 0.6846531968814573*(alpha_vdim[2]*f[22]+alpha_vdim[0]*f[14]); 
  out[37] += 0.3061862178478971*(alpha_cdim[4]*f[38]+alpha_cdim[0]*(f[34]+f[33]))+0.6846531968814573*alpha_vdim[0]*f[16]+0.273861278752583*alpha_cdim[3]*f[8]+0.6846531968814573*alpha_vdim[2]*f[7]; 
  out[38] += 0.3061862178478971*(alpha_cdim[0]*f[35]+alpha_vdim[0]*f[33]+alpha_vdim[1]*f[32])+0.6846531968814573*(alpha_vdim[2]*f[26]+alpha_vdim[0]*f[18])+0.273861278752583*alpha_cdim[3]*f[11]; 
  out[39] += 0.3061862178478971*(alpha_vdim[1]*f[37]+alpha_cdim[0]*f[35]+alpha_vdim[0]*f[34]+alpha_cdim[4]*f[32])+0.6846531968814573*(alpha_vdim[0]*f[19]+alpha_vdim[2]*f[11]); 
  out[40] += 0.3061862178478971*alpha_cdim[0]*f[36]+0.6846531968814573*(alpha_vdim[2]*f[27]+alpha_vdim[0]*f[21])+0.273861278752583*alpha_cdim[3]*f[14]; 
  out[41] += 0.3061862178478971*(alpha_cdim[4]*f[42]+alpha_cdim[0]*f[36])+0.6846531968814573*(alpha_vdim[0]*f[22]+alpha_vdim[2]*f[14]); 
  out[42] += 0.3061862178478971*(alpha_vdim[1]*f[40]+alpha_vdim[0]*f[36])+0.6846531968814573*(alpha_vdim[2]*f[30]+alpha_vdim[0]*f[25]); 
  out[43] += 0.3061862178478971*(alpha_cdim[0]*(f[39]+f[38])+alpha_vdim[0]*f[37]+alpha_vdim[1]*f[34]+alpha_cdim[4]*f[33])+0.6846531968814573*alpha_vdim[0]*f[26]+0.273861278752583*alpha_cdim[3]*f[19]+0.6846531968814573*alpha_vdim[2]*f[18]; 
  out[44] += 0.3061862178478971*(alpha_cdim[4]*f[45]+alpha_cdim[0]*(f[41]+f[40]))+0.6846531968814573*alpha_vdim[0]*f[27]+0.273861278752583*alpha_cdim[3]*f[22]+0.6846531968814573*alpha_vdim[2]*f[21]; 
  out[45] += 0.3061862178478971*(alpha_cdim[0]*f[42]+alpha_vdim[0]*f[40]+alpha_vdim[1]*f[36])+0.6846531968814573*(alpha_vdim[2]*f[31]+alpha_vdim[0]*f[29])+0.273861278752583*alpha_cdim[3]*f[25]; 
  out[46] += 0.3061862178478971*(alpha_vdim[1]*f[44]+alpha_cdim[0]*f[42]+alpha_vdim[0]*f[41]+alpha_cdim[4]*f[36])+0.6846531968814573*(alpha_vdim[0]*f[30]+alpha_vdim[2]*f[25]); 
  out[47] += 0.3061862178478971*(alpha_cdim[0]*(f[46]+f[45])+alpha_vdim[0]*f[44]+alpha_vdim[1]*f[41]+alpha_cdim[4]*f[40])+0.6846531968814573*alpha_vdim[0]*f[31]+0.273861278752583*alpha_cdim[3]*f[30]+0.6846531968814573*alpha_vdim[2]*f[29]; 
  out[48] += 0.6846531968814573*(alpha_vdim[1]*f[9]+alpha_vdim[0]*f[4]); 
  out[49] += 0.3061862178478971*(alpha_cdim[3]*f[51]+alpha_cdim[0]*f[48])+0.6846531968814573*(alpha_vdim[0]*f[9]+alpha_vdim[1]*f[4]); 
  out[50] += 0.3061862178478971*alpha_cdim[0]*f[48]+0.6846531968814573*(alpha_vdim[1]*f[17]+alpha_vdim[0]*f[10])+0.273861278752583*alpha_cdim[4]*f[4]; 
  out[51] += 0.3061862178478971*(alpha_vdim[2]*f[50]+alpha_vdim[0]*f[48])+0.6846531968814573*(alpha_vdim[1]*f[18]+alpha_vdim[0]*f[11]); 
  out[52] += 0.6846531968814573*(alpha_vdim[1]*f[23]+alpha_vdim[0]*f[15]); 
  out[53] += 0.3061862178478971*(alpha_cdim[3]*f[55]+alpha_cdim[0]*(f[50]+f[49]))+0.6846531968814573*(alpha_vdim[0]*f[17]+alpha_vdim[1]*f[10])+0.273861278752583*alpha_cdim[4]*f[9]; 
  out[54] += 0.3061862178478971*(alpha_vdim[2]*f[53]+alpha_cdim[0]*f[51]+alpha_vdim[0]*f[49]+alpha_cdim[3]*f[48])+0.6846531968814573*(alpha_vdim[0]*f[18]+alpha_vdim[1]*f[11]); 
  out[55] += 0.3061862178478971*(alpha_cdim[0]*f[51]+alpha_vdim[0]*f[50]+alpha_vdim[2]*f[48])+0.6846531968814573*(alpha_vdim[1]*f[26]+alpha_vdim[0]*f[19])+0.273861278752583*alpha_cdim[4]*f[11]; 
  out[56] += 0.3061862178478971*(alpha_cdim[3]*f[58]+alpha_cdim[0]*f[52])+0.6846531968814573*(alpha_vdim[0]*f[23]+alpha_vdim[1]*f[15]); 
  out[57] += 0.3061862178478971*alpha_cdim[0]*f[52]+0.6846531968814573*(alpha_vdim[1]*f[28]+alpha_vdim[0]*f[24])+0.273861278752583*alpha_cdim[4]*f[15]; 
  out[58] += 0.3061862178478971*(alpha_vdim[2]*f[57]+alpha_vdim[0]*f[52])+0.6846531968814573*(alpha_vdim[1]*f[29]+alpha_vdim[0]*f[25]); 
  out[59] += 0.3061862178478971*(alpha_cdim[0]*(f[55]+f[54])+alpha_vdim[0]*f[53]+alpha_cdim[3]*f[50]+alpha_vdim[2]*f[49])+0.6846531968814573*(alpha_vdim[0]*f[26]+alpha_vdim[1]*f[19])+0.273861278752583*alpha_cdim[4]*f[18]; 
  out[60] += 0.3061862178478971*(alpha_cdim[3]*f[62]+alpha_cdim[0]*(f[57]+f[56]))+0.6846531968814573*(alpha_vdim[0]*f[28]+alpha_vdim[1]*f[24])+0.273861278752583*alpha_cdim[4]*f[23]; 
  out[61] += 0.3061862178478971*(alpha_vdim[2]*f[60]+alpha_cdim[0]*f[58]+alpha_vdim[0]*f[56]+alpha_cdim[3]*f[52])+0.6846531968814573*(alpha_vdim[0]*f[29]+alpha_vdim[1]*f[25]); 
  out[62] += 0.3061862178478971*(alpha_cdim[0]*f[58]+alpha_vdim[0]*f[57]+alpha_vdim[2]*f[52])+0.6846531968814573*(alpha_vdim[1]*f[31]+alpha_vdim[0]*f[30])+0.273861278752583*alpha_cdim[4]*f[25]; 
  out[63] += 0.3061862178478971*(alpha_cdim[0]*(f[62]+f[61])+alpha_vdim[0]*f[60]+alpha_cdim[3]*f[57]+alpha_vdim[2]*f[56])+0.6846531968814573*(alpha_vdim[0]*f[31]+alpha_vdim[1]*f[30])+0.273861278752583*alpha_cdim[4]*f[29]; 
  out[65] += 0.3061862178478971*(alpha_cdim[3]*f[67]+alpha_cdim[0]*f[64]); 
  out[66] += 0.3061862178478971*(alpha_cdim[4]*f[68]+alpha_cdim[0]*f[64]); 
  out[67] += 0.3061862178478971*(alpha_vdim[2]*f[66]+alpha_vdim[0]*f[64]); 
  out[68] += 0.3061862178478971*(alpha_vdim[1]*f[65]+alpha_vdim[0]*f[64]); 
  out[69] += 0.3061862178478971*(alpha_cdim[4]*f[72]+alpha_cdim[3]*f[71]+alpha_cdim[0]*(f[66]+f[65])); 
  out[70] += 0.3061862178478971*(alpha_vdim[2]*f[69]+alpha_cdim[0]*f[67]+alpha_vdim[0]*f[65]+alpha_cdim[3]*f[64]); 
  out[71] += 0.3061862178478971*(alpha_cdim[4]*f[74]+alpha_cdim[0]*f[67]+alpha_vdim[0]*f[66]+alpha_vdim[2]*f[64]); 
  out[72] += 0.3061862178478971*(alpha_cdim[3]*f[74]+alpha_cdim[0]*f[68]+alpha_vdim[0]*f[65]+alpha_vdim[1]*f[64]); 
  out[73] += 0.3061862178478971*(alpha_vdim[1]*f[69]+alpha_cdim[0]*f[68]+alpha_vdim[0]*f[66]+alpha_cdim[4]*f[64]); 
  out[74] += 0.3061862178478971*(alpha_vdim[2]*f[73]+alpha_vdim[1]*f[70]+alpha_vdim[0]*(f[68]+f[67])); 
  out[75] += 0.3061862178478971*(alpha_cdim[4]*f[77]+alpha_cdim[0]*(f[71]+f[70])+alpha_vdim[0]*f[69]+alpha_cdim[3]*f[66]+alpha_vdim[2]*f[65]); 
  out[76] += 0.3061862178478971*(alpha_cdim[3]*f[78]+alpha_cdim[0]*(f[73]+f[72])+alpha_vdim[0]*f[69]+alpha_vdim[1]*f[66]+alpha_cdim[4]*f[65]); 
  out[77] += 0.3061862178478971*(alpha_vdim[2]*f[76]+alpha_cdim[0]*f[74]+alpha_vdim[0]*(f[72]+f[70])+alpha_cdim[3]*f[68]+alpha_vdim[1]*f[67]); 
  out[78] += 0.3061862178478971*(alpha_vdim[1]*f[75]+alpha_cdim[0]*f[74]+alpha_vdim[0]*(f[73]+f[71])+alpha_vdim[2]*f[68]+alpha_cdim[4]*f[67]); 
  out[79] += 0.3061862178478971*(alpha_cdim[0]*(f[78]+f[77])+alpha_vdim[0]*(f[76]+f[75])+alpha_cdim[3]*f[73]+alpha_vdim[2]*f[72]+alpha_vdim[1]*f[71]+alpha_cdim[4]*f[70]); 

  return cflFreq_mid; 
} 
