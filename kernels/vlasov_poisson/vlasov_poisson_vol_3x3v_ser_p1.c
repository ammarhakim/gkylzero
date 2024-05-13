#include <gkyl_vlasov_poisson_kernels.h> 

GKYL_CU_DH double vlasov_poisson_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:     potential (scaled by appropriate factors).
  // f:         Input distribution function.
  // out:       Incremented output.

  const double dv10 = 2/dxv[3]; 
  const double dv1 = dxv[3], wv1 = w[3]; 
  const double dv11 = 2/dxv[4]; 
  const double dv2 = dxv[4], wv2 = w[4]; 
  const double dv12 = 2/dxv[5]; 
  const double dv3 = dxv[5], wv3 = w[5]; 
  double dv0dx0 = dxv[3]/dxv[0]; 
  double w0dx0 = w[3]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 
  double dv1dx1 = dxv[4]/dxv[1]; 
  double w1dx1 = w[4]/dxv[1]; 
  const double dx11 = 2/dxv[1]; 
  double dv2dx2 = dxv[5]/dxv[2]; 
  double w2dx2 = w[5]/dxv[2]; 
  const double dx12 = 2/dxv[2]; 

  const double *phi = &field[0]; 

  double cflFreq_mid = 0.0; 

  double alpha_cdim[480] = {0.0}; 
  double alpha_vdim[480] = {0.0}; 

  alpha_cdim[0] = 16.0*w0dx0; 
  alpha_cdim[4] = 4.618802153517007*dv0dx0; 

  cflFreq_mid += 3.0*(fabs(w0dx0)+0.5*dv0dx0); 

  alpha_cdim[160] = 16.0*w1dx1; 
  alpha_cdim[165] = 4.618802153517007*dv1dx1; 

  cflFreq_mid += 3.0*(fabs(w1dx1)+0.5*dv1dx1); 

  alpha_cdim[320] = 16.0*w2dx2; 
  alpha_cdim[326] = 4.618802153517007*dv2dx2; 

  cflFreq_mid += 3.0*(fabs(w2dx2)+0.5*dv2dx2); 

  alpha_vdim[0] = -4.898979485566357*phi[1]*dv10*dx10; 
  alpha_vdim[2] = -4.898979485566357*phi[4]*dv10*dx10; 
  alpha_vdim[3] = -4.898979485566357*phi[5]*dv10*dx10; 
  alpha_vdim[9] = -4.898979485566357*phi[7]*dv10*dx10; 

  cflFreq_mid += 5.0*fabs(0.0625*alpha_vdim[0]); 

  alpha_vdim[160] = -4.898979485566357*phi[2]*dv11*dx11; 
  alpha_vdim[161] = -4.898979485566357*phi[4]*dv11*dx11; 
  alpha_vdim[163] = -4.898979485566357*phi[6]*dv11*dx11; 
  alpha_vdim[168] = -4.898979485566357*phi[7]*dv11*dx11; 

  cflFreq_mid += 5.0*fabs(0.0625*alpha_vdim[0]); 

  alpha_vdim[320] = -4.898979485566357*phi[3]*dv12*dx12; 
  alpha_vdim[321] = -4.898979485566357*phi[5]*dv12*dx12; 
  alpha_vdim[322] = -4.898979485566357*phi[6]*dv12*dx12; 
  alpha_vdim[327] = -4.898979485566357*phi[7]*dv12*dx12; 

  cflFreq_mid += 5.0*fabs(0.0625*alpha_vdim[0]); 

  out[1] += 0.2165063509461096*(alpha_cdim[4]*f[4]+alpha_cdim[0]*f[0]); 
  out[2] += 0.2165063509461096*(alpha_cdim[5]*f[5]+alpha_cdim[0]*f[0]); 
  out[3] += 0.2165063509461096*(alpha_cdim[6]*f[6]+alpha_cdim[0]*f[0]); 
  out[4] += 0.2165063509461096*(alpha_vdim[9]*f[9]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[0]*f[0]); 
  out[5] += 0.2165063509461096*(alpha_vdim[8]*f[8]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[6] += 0.2165063509461096*(alpha_vdim[7]*f[7]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[7] += 0.2165063509461096*(alpha_cdim[5]*f[13]+alpha_cdim[4]*f[11]+alpha_cdim[0]*(f[2]+f[1])); 
  out[8] += 0.2165063509461096*(alpha_cdim[6]*f[17]+alpha_cdim[4]*f[12]+alpha_cdim[0]*(f[3]+f[1])); 
  out[9] += 0.2165063509461096*(alpha_cdim[6]*f[18]+alpha_cdim[5]*f[15]+alpha_cdim[0]*(f[3]+f[2])); 
  out[10] += 0.1936491673103708*alpha_cdim[4]*f[64]+0.2165063509461096*(alpha_vdim[9]*f[22]+alpha_vdim[3]*f[8]+alpha_vdim[2]*f[7]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]+alpha_vdim[0]*f[1]); 
  out[11] += 0.2165063509461096*(alpha_cdim[5]*f[16]+alpha_vdim[3]*f[9]+f[3]*alpha_vdim[9]+alpha_cdim[0]*f[4]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[12] += 0.2165063509461096*(alpha_cdim[6]*f[20]+alpha_vdim[2]*f[9]+f[2]*alpha_vdim[9]+alpha_cdim[0]*f[4]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[13] += 0.2165063509461096*(alpha_cdim[4]*f[16]+alpha_vdim[3]*f[8]+f[3]*alpha_vdim[8]+alpha_cdim[0]*f[5]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[14] += 0.1936491673103708*alpha_cdim[5]*f[96]+0.2165063509461096*(alpha_vdim[8]*f[22]+alpha_vdim[3]*f[9]+alpha_vdim[1]*f[7]+alpha_cdim[0]*f[5]+f[0]*alpha_cdim[5]+alpha_vdim[0]*f[2]); 
  out[15] += 0.2165063509461096*(alpha_cdim[6]*f[21]+alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8]+alpha_cdim[0]*f[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[16] += 0.2165063509461096*(alpha_vdim[9]*f[28]+alpha_vdim[8]*f[24]+alpha_vdim[3]*f[15]+alpha_vdim[2]*f[14]+alpha_vdim[3]*f[12]+alpha_vdim[1]*f[10]+alpha_vdim[0]*(f[5]+f[4])); 
  out[17] += 0.2165063509461096*(alpha_cdim[4]*f[20]+alpha_vdim[2]*f[7]+f[2]*alpha_vdim[7]+alpha_cdim[0]*f[6]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[18] += 0.2165063509461096*(alpha_cdim[5]*f[21]+alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7]+alpha_cdim[0]*f[6]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[19] += 0.1936491673103708*alpha_cdim[6]*f[128]+0.2165063509461096*(alpha_vdim[7]*f[22]+alpha_vdim[2]*f[9]+alpha_vdim[1]*f[8]+alpha_cdim[0]*f[6]+f[0]*alpha_cdim[6]+alpha_vdim[0]*f[3]); 
  out[20] += 0.2165063509461096*(alpha_vdim[9]*f[34]+alpha_vdim[7]*f[23]+alpha_vdim[3]*f[19]+alpha_vdim[2]*(f[18]+f[11])+alpha_vdim[1]*f[10]+alpha_vdim[0]*(f[6]+f[4])); 
  out[21] += 0.2165063509461096*(alpha_vdim[8]*f[33]+alpha_vdim[7]*f[26]+alpha_vdim[3]*f[19]+alpha_vdim[1]*f[17]+alpha_vdim[2]*f[14]+alpha_vdim[1]*f[13]+alpha_vdim[0]*(f[6]+f[5])); 
  out[22] += 0.2165063509461096*(alpha_cdim[6]*f[32]+alpha_cdim[5]*f[27]+alpha_cdim[4]*f[25]+alpha_cdim[0]*(f[9]+f[8]+f[7])); 
  out[23] += 0.1936491673103708*alpha_cdim[4]*f[66]+0.2165063509461096*(alpha_cdim[5]*f[29]+alpha_vdim[3]*f[22]+alpha_cdim[0]*(f[11]+f[10])+f[8]*alpha_vdim[9]+alpha_vdim[0]*f[7]+f[2]*alpha_cdim[4]+f[1]*alpha_vdim[2]); 
  out[24] += 0.1936491673103708*alpha_cdim[4]*f[67]+0.2165063509461096*(alpha_cdim[6]*f[35]+alpha_vdim[2]*f[22]+alpha_cdim[0]*(f[12]+f[10])+f[7]*alpha_vdim[9]+alpha_vdim[0]*f[8]+f[3]*alpha_cdim[4]+f[1]*alpha_vdim[3]); 
  out[25] += 0.2165063509461096*(alpha_cdim[6]*f[36]+alpha_cdim[5]*f[31]+alpha_cdim[0]*(f[12]+f[11])+alpha_vdim[0]*f[9]+f[0]*alpha_vdim[9]+alpha_vdim[2]*f[3]+f[2]*alpha_vdim[3]); 
  out[26] += 0.1936491673103708*alpha_cdim[5]*f[97]+0.2165063509461096*(alpha_cdim[4]*f[30]+alpha_vdim[3]*f[22]+alpha_cdim[0]*(f[14]+f[13])+alpha_vdim[8]*f[9]+alpha_vdim[0]*f[7]+f[1]*alpha_cdim[5]+alpha_vdim[1]*f[2]); 
  out[27] += 0.2165063509461096*(alpha_cdim[6]*f[38]+alpha_cdim[4]*f[31]+alpha_cdim[0]*(f[15]+f[13])+alpha_vdim[0]*f[8]+f[0]*alpha_vdim[8]+alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]); 
  out[28] += 0.1936491673103708*alpha_cdim[5]*f[99]+0.2165063509461096*(alpha_cdim[6]*f[39]+alpha_vdim[1]*f[22]+alpha_cdim[0]*(f[15]+f[14])+alpha_vdim[0]*f[9]+f[7]*alpha_vdim[8]+f[3]*alpha_cdim[5]+f[2]*alpha_vdim[3]); 
  out[29] += 0.1936491673103708*alpha_cdim[4]*f[68]+0.2165063509461096*(alpha_vdim[9]*f[43]+alpha_vdim[3]*f[27]+alpha_vdim[2]*f[26]+alpha_vdim[3]*f[24]+alpha_cdim[0]*f[16]+alpha_vdim[0]*f[13]+alpha_vdim[8]*f[12]+alpha_vdim[0]*f[10]+alpha_cdim[4]*f[5]+alpha_vdim[1]*f[4]); 
  out[30] += 0.1936491673103708*alpha_cdim[5]*f[100]+0.2165063509461096*(alpha_vdim[8]*f[42]+alpha_vdim[3]*(f[28]+f[25])+alpha_vdim[1]*f[23]+alpha_cdim[0]*f[16]+alpha_vdim[9]*f[15]+alpha_vdim[0]*(f[14]+f[11])+alpha_vdim[2]*f[5]+f[4]*alpha_cdim[5]); 
  out[31] += 0.2165063509461096*(alpha_cdim[6]*f[41]+alpha_vdim[2]*f[28]+alpha_vdim[1]*f[24]+alpha_cdim[0]*f[16]+alpha_vdim[0]*f[15]+alpha_vdim[9]*f[14]+alpha_vdim[0]*f[12]+alpha_vdim[8]*f[10]+alpha_vdim[3]*(f[5]+f[4])); 
  out[32] += 0.2165063509461096*(alpha_cdim[5]*f[38]+alpha_cdim[4]*f[36]+alpha_cdim[0]*(f[18]+f[17])+alpha_vdim[0]*f[7]+f[0]*alpha_vdim[7]+alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[33] += 0.1936491673103708*alpha_cdim[6]*f[129]+0.2165063509461096*(alpha_cdim[4]*f[37]+alpha_vdim[2]*f[22]+alpha_cdim[0]*(f[19]+f[17])+alpha_vdim[7]*f[9]+alpha_vdim[0]*f[8]+f[1]*alpha_cdim[6]+alpha_vdim[1]*f[3]); 
  out[34] += 0.1936491673103708*alpha_cdim[6]*f[130]+0.2165063509461096*(alpha_cdim[5]*f[40]+alpha_vdim[1]*f[22]+alpha_cdim[0]*(f[19]+f[18])+alpha_vdim[0]*f[9]+alpha_vdim[7]*f[8]+f[2]*alpha_cdim[6]+alpha_vdim[2]*f[3]); 
  out[35] += 0.1936491673103708*alpha_cdim[4]*f[69]+0.2165063509461096*(alpha_vdim[9]*f[47]+alpha_vdim[3]*f[33]+alpha_vdim[2]*(f[32]+f[23])+alpha_cdim[0]*f[20]+alpha_vdim[0]*f[17]+alpha_vdim[7]*f[11]+alpha_vdim[0]*f[10]+alpha_cdim[4]*f[6]+alpha_vdim[1]*f[4]); 
  out[36] += 0.2165063509461096*(alpha_cdim[5]*f[41]+alpha_vdim[3]*f[34]+alpha_vdim[1]*f[23]+alpha_cdim[0]*f[20]+alpha_vdim[9]*f[19]+alpha_vdim[0]*(f[18]+f[11])+alpha_vdim[7]*f[10]+alpha_vdim[2]*(f[6]+f[4])); 
  out[37] += 0.1936491673103708*alpha_cdim[6]*f[132]+0.2165063509461096*(alpha_vdim[7]*f[42]+alpha_vdim[2]*(f[34]+f[25])+alpha_vdim[1]*f[24]+alpha_cdim[0]*f[20]+alpha_vdim[0]*f[19]+alpha_vdim[9]*f[18]+alpha_vdim[0]*f[12]+alpha_vdim[3]*f[6]+f[4]*alpha_cdim[6]); 
  out[38] += 0.2165063509461096*(alpha_cdim[4]*f[41]+alpha_vdim[3]*f[33]+alpha_vdim[2]*f[26]+alpha_cdim[0]*f[21]+alpha_vdim[8]*f[19]+alpha_vdim[0]*f[17]+alpha_vdim[7]*f[14]+alpha_vdim[0]*f[13]+alpha_vdim[1]*(f[6]+f[5])); 
  out[39] += 0.1936491673103708*alpha_cdim[5]*f[101]+0.2165063509461096*(alpha_vdim[8]*f[47]+alpha_vdim[3]*f[34]+alpha_vdim[1]*(f[32]+f[26])+alpha_cdim[0]*f[21]+alpha_vdim[0]*(f[18]+f[14])+alpha_vdim[7]*f[13]+alpha_cdim[5]*f[6]+alpha_vdim[2]*f[5]); 
  out[40] += 0.1936491673103708*alpha_cdim[6]*f[133]+0.2165063509461096*(alpha_vdim[7]*f[43]+alpha_vdim[1]*f[33]+alpha_vdim[2]*f[28]+alpha_vdim[1]*f[27]+alpha_cdim[0]*f[21]+alpha_vdim[0]*f[19]+alpha_vdim[8]*f[17]+alpha_vdim[0]*f[15]+alpha_vdim[3]*f[6]+f[5]*alpha_cdim[6]); 
  out[41] += 0.2165063509461096*(alpha_vdim[9]*f[53]+alpha_vdim[8]*f[49]+alpha_vdim[7]*f[44]+alpha_vdim[3]*f[40]+alpha_vdim[2]*f[39]+alpha_vdim[3]*f[37]+alpha_vdim[1]*f[35]+alpha_vdim[2]*f[30]+alpha_vdim[1]*f[29]+alpha_vdim[0]*(f[21]+f[20]+f[16])); 
  out[42] += 0.1936491673103708*alpha_cdim[4]*f[72]+0.2165063509461096*(alpha_cdim[6]*f[48]+alpha_cdim[5]*f[45]+alpha_cdim[0]*(f[25]+f[24]+f[23])+alpha_vdim[0]*f[22]+alpha_cdim[4]*f[9]+f[1]*alpha_vdim[9]+alpha_vdim[2]*f[8]+alpha_vdim[3]*f[7]); 
  out[43] += 0.1936491673103708*alpha_cdim[5]*f[103]+0.2165063509461096*(alpha_cdim[6]*f[51]+alpha_cdim[4]*f[46]+alpha_cdim[0]*(f[28]+f[27]+f[26])+alpha_vdim[0]*f[22]+alpha_vdim[1]*f[9]+alpha_cdim[5]*f[8]+f[2]*alpha_vdim[8]+alpha_vdim[3]*f[7]); 
  out[44] += 0.1936491673103708*(alpha_cdim[5]*f[105]+alpha_cdim[4]*f[74])+0.2165063509461096*(alpha_vdim[3]*(f[43]+f[42])+alpha_cdim[0]*(f[30]+f[29])+alpha_vdim[9]*f[27]+alpha_vdim[0]*f[26]+alpha_vdim[8]*f[25]+alpha_vdim[0]*f[23]+alpha_cdim[4]*f[14]+alpha_vdim[2]*f[13]+alpha_vdim[1]*f[11]+alpha_cdim[5]*f[10]); 
  out[45] += 0.1936491673103708*alpha_cdim[4]*f[75]+0.2165063509461096*(alpha_cdim[6]*f[54]+alpha_vdim[2]*f[43]+alpha_cdim[0]*(f[31]+f[29])+alpha_vdim[0]*f[27]+alpha_vdim[9]*f[26]+alpha_vdim[0]*f[24]+alpha_cdim[4]*f[15]+alpha_vdim[3]*f[13]+alpha_vdim[1]*f[12]+alpha_vdim[3]*f[10]+f[4]*alpha_vdim[8]); 
  out[46] += 0.1936491673103708*alpha_cdim[5]*f[107]+0.2165063509461096*(alpha_cdim[6]*f[55]+alpha_vdim[1]*f[42]+alpha_cdim[0]*(f[31]+f[30])+alpha_vdim[0]*(f[28]+f[25])+alpha_vdim[8]*f[23]+alpha_vdim[2]*f[15]+alpha_vdim[3]*f[14]+alpha_cdim[5]*f[12]+alpha_vdim[3]*f[11]+f[5]*alpha_vdim[9]); 
  out[47] += 0.1936491673103708*alpha_cdim[6]*f[134]+0.2165063509461096*(alpha_cdim[5]*f[52]+alpha_cdim[4]*f[50]+alpha_cdim[0]*(f[34]+f[33]+f[32])+alpha_vdim[0]*f[22]+alpha_vdim[1]*f[9]+alpha_vdim[2]*f[8]+alpha_cdim[6]*f[7]+f[3]*alpha_vdim[7]); 
  out[48] += 0.1936491673103708*alpha_cdim[4]*f[77]+0.2165063509461096*(alpha_cdim[5]*f[54]+alpha_vdim[3]*f[47]+alpha_cdim[0]*(f[36]+f[35])+alpha_vdim[9]*f[33]+alpha_vdim[0]*(f[32]+f[23])+alpha_cdim[4]*f[18]+alpha_vdim[2]*f[17]+alpha_vdim[1]*f[11]+alpha_vdim[2]*f[10]+f[4]*alpha_vdim[7]); 
  out[49] += 0.1936491673103708*(alpha_cdim[6]*f[137]+alpha_cdim[4]*f[78])+0.2165063509461096*(alpha_vdim[2]*(f[47]+f[42])+alpha_cdim[0]*(f[37]+f[35])+alpha_vdim[0]*f[33]+alpha_vdim[9]*f[32]+alpha_vdim[7]*f[25]+alpha_vdim[0]*f[24]+alpha_cdim[4]*f[19]+alpha_vdim[3]*f[17]+alpha_vdim[1]*f[12]+alpha_cdim[6]*f[10]); 
  out[50] += 0.1936491673103708*alpha_cdim[6]*f[138]+0.2165063509461096*(alpha_cdim[5]*f[56]+alpha_vdim[1]*f[42]+alpha_cdim[0]*(f[37]+f[36])+alpha_vdim[0]*(f[34]+f[25])+alpha_vdim[7]*f[24]+alpha_vdim[2]*f[19]+alpha_vdim[3]*f[18]+alpha_vdim[2]*f[12]+alpha_cdim[6]*f[11]+f[6]*alpha_vdim[9]); 
  out[51] += 0.1936491673103708*alpha_cdim[5]*f[108]+0.2165063509461096*(alpha_cdim[4]*f[55]+alpha_vdim[3]*f[47]+alpha_cdim[0]*(f[39]+f[38])+alpha_vdim[8]*f[34]+alpha_vdim[0]*(f[32]+f[26])+alpha_vdim[1]*f[18]+alpha_cdim[5]*f[17]+alpha_vdim[1]*f[14]+alpha_vdim[2]*f[13]+f[5]*alpha_vdim[7]); 
  out[52] += 0.1936491673103708*alpha_cdim[6]*f[140]+0.2165063509461096*(alpha_cdim[4]*f[56]+alpha_vdim[2]*f[43]+alpha_cdim[0]*(f[40]+f[38])+alpha_vdim[0]*f[33]+alpha_vdim[7]*f[28]+alpha_vdim[0]*f[27]+alpha_vdim[1]*f[19]+alpha_vdim[3]*f[17]+alpha_vdim[1]*f[15]+alpha_cdim[6]*f[13]+f[6]*alpha_vdim[8]); 
  out[53] += 0.1936491673103708*(alpha_cdim[6]*f[141]+alpha_cdim[5]*f[110])+0.2165063509461096*(alpha_vdim[1]*(f[47]+f[43])+alpha_cdim[0]*(f[40]+f[39])+alpha_vdim[0]*f[34]+alpha_vdim[8]*f[32]+alpha_vdim[0]*f[28]+alpha_vdim[7]*f[27]+alpha_cdim[5]*f[19]+alpha_vdim[3]*f[18]+alpha_vdim[2]*f[15]+alpha_cdim[6]*f[14]); 
  out[54] += 0.1936491673103708*alpha_cdim[4]*f[79]+0.2165063509461096*(alpha_vdim[9]*f[59]+alpha_vdim[3]*f[52]+alpha_vdim[2]*f[51]+alpha_vdim[3]*f[49]+alpha_vdim[2]*f[44]+alpha_cdim[0]*f[41]+alpha_vdim[0]*f[38]+alpha_vdim[8]*f[37]+alpha_vdim[0]*f[35]+alpha_vdim[7]*f[30]+alpha_vdim[0]*f[29]+alpha_cdim[4]*f[21]+alpha_vdim[1]*(f[20]+f[16])); 
  out[55] += 0.1936491673103708*alpha_cdim[5]*f[111]+0.2165063509461096*(alpha_vdim[8]*f[58]+alpha_vdim[3]*(f[53]+f[50])+alpha_vdim[1]*(f[48]+f[44])+alpha_cdim[0]*f[41]+alpha_vdim[9]*f[40]+alpha_vdim[0]*(f[39]+f[36]+f[30])+alpha_vdim[7]*f[29]+alpha_vdim[2]*f[21]+alpha_cdim[5]*f[20]+alpha_vdim[2]*f[16]); 
  out[56] += 0.1936491673103708*alpha_cdim[6]*f[143]+0.2165063509461096*(alpha_vdim[7]*f[57]+alpha_vdim[2]*f[53]+alpha_vdim[1]*f[49]+alpha_vdim[2]*f[46]+alpha_vdim[1]*f[45]+alpha_cdim[0]*f[41]+alpha_vdim[0]*f[40]+alpha_vdim[9]*f[39]+alpha_vdim[0]*f[37]+alpha_vdim[8]*f[35]+alpha_vdim[0]*f[31]+alpha_vdim[3]*(f[21]+f[20])+alpha_cdim[6]*f[16]); 
  out[57] += 0.1936491673103708*(alpha_cdim[5]*f[114]+alpha_cdim[4]*f[83])+0.2165063509461096*(alpha_cdim[6]*f[60]+alpha_cdim[0]*(f[46]+f[45]+f[44])+alpha_vdim[0]*(f[43]+f[42])+alpha_cdim[4]*f[28]+alpha_vdim[2]*f[27]+alpha_vdim[3]*f[26]+alpha_vdim[1]*f[25]+alpha_cdim[5]*f[24]+alpha_vdim[3]*f[23]+alpha_vdim[9]*f[13]+alpha_vdim[8]*f[11]); 
  out[58] += 0.1936491673103708*(alpha_cdim[6]*f[145]+alpha_cdim[4]*f[86])+0.2165063509461096*(alpha_cdim[5]*f[61]+alpha_cdim[0]*(f[50]+f[49]+f[48])+alpha_vdim[0]*(f[47]+f[42])+alpha_cdim[4]*f[34]+alpha_vdim[2]*f[33]+alpha_vdim[3]*f[32]+alpha_vdim[1]*f[25]+alpha_vdim[2]*f[24]+alpha_cdim[6]*f[23]+alpha_vdim[9]*f[17]+alpha_vdim[7]*f[12]); 
  out[59] += 0.1936491673103708*(alpha_cdim[6]*f[148]+alpha_cdim[5]*f[117])+0.2165063509461096*(alpha_cdim[4]*f[62]+alpha_cdim[0]*(f[53]+f[52]+f[51])+alpha_vdim[0]*(f[47]+f[43])+alpha_vdim[1]*f[34]+alpha_cdim[5]*f[33]+alpha_vdim[3]*f[32]+alpha_vdim[1]*f[28]+alpha_vdim[2]*f[27]+alpha_cdim[6]*f[26]+alpha_vdim[8]*f[18]+alpha_vdim[7]*f[15]); 
  out[60] += 0.1936491673103708*(alpha_cdim[5]*f[119]+alpha_cdim[4]*f[88])+0.2165063509461096*(alpha_vdim[3]*(f[59]+f[58])+alpha_cdim[0]*(f[55]+f[54])+alpha_vdim[9]*f[52]+alpha_vdim[0]*f[51]+alpha_vdim[8]*f[50]+alpha_vdim[0]*(f[48]+f[44])+alpha_cdim[4]*f[39]+alpha_vdim[2]*f[38]+alpha_vdim[1]*f[36]+alpha_cdim[5]*f[35]+alpha_vdim[1]*f[30]+alpha_vdim[2]*f[29]+alpha_vdim[7]*f[16]); 
  out[61] += 0.1936491673103708*(alpha_cdim[6]*f[151]+alpha_cdim[4]*f[89])+0.2165063509461096*(alpha_vdim[2]*(f[59]+f[57])+alpha_cdim[0]*(f[56]+f[54])+alpha_vdim[0]*f[52]+alpha_vdim[9]*f[51]+alpha_vdim[0]*f[49]+alpha_vdim[7]*f[46]+alpha_vdim[0]*f[45]+alpha_cdim[4]*f[40]+alpha_vdim[3]*f[38]+alpha_vdim[1]*f[37]+alpha_vdim[3]*f[35]+alpha_vdim[1]*f[31]+alpha_cdim[6]*f[29]+alpha_vdim[8]*f[20]); 
  out[62] += 0.1936491673103708*(alpha_cdim[6]*f[152]+alpha_cdim[5]*f[121])+0.2165063509461096*(alpha_vdim[1]*(f[58]+f[57])+alpha_cdim[0]*(f[56]+f[55])+alpha_vdim[0]*(f[53]+f[50])+alpha_vdim[8]*f[48]+alpha_vdim[0]*f[46]+alpha_vdim[7]*f[45]+alpha_vdim[2]*f[40]+alpha_vdim[3]*f[39]+alpha_cdim[5]*f[37]+alpha_vdim[3]*f[36]+alpha_vdim[2]*f[31]+alpha_cdim[6]*f[30]+alpha_vdim[9]*f[21]); 
  out[63] += 0.1936491673103708*(alpha_cdim[6]*f[156]+alpha_cdim[5]*f[125]+alpha_cdim[4]*f[94])+0.2165063509461096*(alpha_cdim[0]*(f[62]+f[61]+f[60])+alpha_vdim[0]*(f[59]+f[58]+f[57])+alpha_cdim[4]*f[53]+alpha_vdim[2]*f[52]+alpha_vdim[3]*f[51]+alpha_vdim[1]*f[50]+alpha_cdim[5]*f[49]+alpha_vdim[3]*f[48]+alpha_vdim[1]*f[46]+alpha_vdim[2]*f[45]+alpha_cdim[6]*f[44]+alpha_vdim[9]*f[38]+alpha_vdim[8]*f[36]+alpha_vdim[7]*f[31]); 
  out[64] += 0.4841229182759271*(alpha_vdim[9]*f[25]+alpha_vdim[3]*f[12]+alpha_vdim[2]*f[11]+alpha_vdim[0]*f[4]); 
  out[65] += 0.2165063509461096*alpha_cdim[0]*f[64]+0.4841229182759271*(alpha_vdim[9]*f[42]+alpha_vdim[3]*f[24]+alpha_vdim[2]*f[23]+alpha_vdim[0]*f[10])+0.1936491673103708*alpha_cdim[4]*f[4]; 
  out[66] += 0.2165063509461096*(alpha_cdim[5]*f[68]+alpha_cdim[0]*f[64])+0.4841229182759271*(alpha_vdim[3]*f[25]+alpha_vdim[9]*f[12]+alpha_vdim[0]*f[11]+alpha_vdim[2]*f[4]); 
  out[67] += 0.2165063509461096*(alpha_cdim[6]*f[69]+alpha_cdim[0]*f[64])+0.4841229182759271*(alpha_vdim[2]*f[25]+alpha_vdim[0]*f[12]+alpha_vdim[9]*f[11]+alpha_vdim[3]*f[4]); 
  out[68] += 0.2165063509461096*(alpha_vdim[8]*f[71]+alpha_vdim[3]*f[67]+alpha_vdim[1]*f[65]+alpha_vdim[0]*f[64])+0.4841229182759271*(alpha_vdim[9]*f[46]+alpha_vdim[3]*f[31]+alpha_vdim[2]*f[30]+alpha_vdim[0]*f[16]); 
  out[69] += 0.2165063509461096*(alpha_vdim[7]*f[70]+alpha_vdim[2]*f[66]+alpha_vdim[1]*f[65]+alpha_vdim[0]*f[64])+0.4841229182759271*(alpha_vdim[9]*f[50]+alpha_vdim[3]*f[37]+alpha_vdim[2]*f[36]+alpha_vdim[0]*f[20]); 
  out[70] += 0.2165063509461096*(alpha_cdim[5]*f[73]+alpha_cdim[0]*(f[66]+f[65]))+0.4841229182759271*(alpha_vdim[3]*f[42]+alpha_vdim[9]*f[24]+alpha_vdim[0]*f[23])+0.1936491673103708*alpha_cdim[4]*f[11]+0.4841229182759271*alpha_vdim[2]*f[10]; 
  out[71] += 0.2165063509461096*(alpha_cdim[6]*f[76]+alpha_cdim[0]*(f[67]+f[65]))+0.4841229182759271*(alpha_vdim[2]*f[42]+alpha_vdim[0]*f[24]+alpha_vdim[9]*f[23])+0.1936491673103708*alpha_cdim[4]*f[12]+0.4841229182759271*alpha_vdim[3]*f[10]; 
  out[72] += 0.2165063509461096*(alpha_cdim[6]*f[77]+alpha_cdim[5]*f[75]+alpha_cdim[0]*(f[67]+f[66]))+0.4841229182759271*(alpha_vdim[0]*f[25]+alpha_vdim[2]*f[12]+alpha_vdim[3]*f[11]+f[4]*alpha_vdim[9]); 
  out[73] += 0.2165063509461096*(alpha_vdim[3]*f[71]+alpha_cdim[0]*f[68]+alpha_vdim[8]*f[67]+alpha_vdim[0]*f[65]+alpha_vdim[1]*f[64])+0.4841229182759271*(alpha_vdim[9]*f[57]+alpha_vdim[3]*f[45]+alpha_vdim[2]*f[44]+alpha_vdim[0]*f[29])+0.1936491673103708*alpha_cdim[4]*f[16]; 
  out[74] += 0.2165063509461096*(alpha_vdim[8]*f[80]+alpha_vdim[3]*f[72]+alpha_vdim[1]*f[70]+alpha_cdim[0]*f[68]+alpha_vdim[0]*f[66]+alpha_cdim[5]*f[64])+0.4841229182759271*(alpha_vdim[3]*f[46]+alpha_vdim[9]*f[31]+alpha_vdim[0]*f[30]+alpha_vdim[2]*f[16]); 
  out[75] += 0.2165063509461096*(alpha_cdim[6]*f[79]+alpha_vdim[1]*f[71]+alpha_cdim[0]*f[68]+alpha_vdim[0]*f[67]+alpha_vdim[8]*f[65]+alpha_vdim[3]*f[64])+0.4841229182759271*(alpha_vdim[2]*f[46]+alpha_vdim[0]*f[31]+alpha_vdim[9]*f[30]+alpha_vdim[3]*f[16]); 
  out[76] += 0.2165063509461096*(alpha_vdim[2]*f[70]+alpha_cdim[0]*f[69]+alpha_vdim[7]*f[66]+alpha_vdim[0]*f[65]+alpha_vdim[1]*f[64])+0.4841229182759271*(alpha_vdim[9]*f[58]+alpha_vdim[3]*f[49]+alpha_vdim[2]*f[48]+alpha_vdim[0]*f[35])+0.1936491673103708*alpha_cdim[4]*f[20]; 
  out[77] += 0.2165063509461096*(alpha_cdim[5]*f[79]+alpha_vdim[1]*f[70]+alpha_cdim[0]*f[69]+alpha_vdim[0]*f[66]+alpha_vdim[7]*f[65]+alpha_vdim[2]*f[64])+0.4841229182759271*(alpha_vdim[3]*f[50]+alpha_vdim[9]*f[37]+alpha_vdim[0]*f[36]+alpha_vdim[2]*f[20]); 
  out[78] += 0.2165063509461096*(alpha_vdim[7]*f[80]+alpha_vdim[2]*f[72]+alpha_vdim[1]*f[71]+alpha_cdim[0]*f[69]+alpha_vdim[0]*f[67]+alpha_cdim[6]*f[64])+0.4841229182759271*(alpha_vdim[2]*f[50]+alpha_vdim[0]*f[37]+alpha_vdim[9]*f[36]+alpha_vdim[3]*f[20]); 
  out[79] += 0.2165063509461096*(alpha_vdim[8]*f[85]+alpha_vdim[7]*f[81]+alpha_vdim[3]*f[78]+alpha_vdim[1]*f[76]+alpha_vdim[2]*f[74]+alpha_vdim[1]*f[73]+alpha_vdim[0]*(f[69]+f[68]))+0.4841229182759271*(alpha_vdim[9]*f[62]+alpha_vdim[3]*f[56]+alpha_vdim[2]*f[55]+alpha_vdim[0]*f[41]); 
  out[80] += 0.2165063509461096*(alpha_cdim[6]*f[84]+alpha_cdim[5]*f[82]+alpha_cdim[0]*(f[72]+f[71]+f[70]))+0.4841229182759271*alpha_vdim[0]*f[42]+0.1936491673103708*alpha_cdim[4]*f[25]+0.4841229182759271*(alpha_vdim[2]*f[24]+alpha_vdim[3]*f[23]+alpha_vdim[9]*f[10]); 
  out[81] += 0.2165063509461096*(alpha_vdim[3]*f[80]+alpha_cdim[0]*(f[74]+f[73])+alpha_vdim[8]*f[72]+alpha_vdim[0]*f[70]+alpha_vdim[1]*f[66]+alpha_cdim[5]*f[65])+0.4841229182759271*(alpha_vdim[3]*f[57]+alpha_vdim[9]*f[45]+alpha_vdim[0]*f[44])+0.1936491673103708*alpha_cdim[4]*f[30]+0.4841229182759271*alpha_vdim[2]*f[29]; 
  out[82] += 0.2165063509461096*(alpha_cdim[6]*f[87]+alpha_cdim[0]*(f[75]+f[73])+alpha_vdim[0]*f[71]+alpha_vdim[1]*f[67]+alpha_vdim[3]*f[65]+alpha_vdim[8]*f[64])+0.4841229182759271*(alpha_vdim[2]*f[57]+alpha_vdim[0]*f[45]+alpha_vdim[9]*f[44])+0.1936491673103708*alpha_cdim[4]*f[31]+0.4841229182759271*alpha_vdim[3]*f[29]; 
  out[83] += 0.2165063509461096*(alpha_cdim[6]*f[88]+alpha_vdim[1]*f[80]+alpha_cdim[0]*(f[75]+f[74])+alpha_vdim[0]*f[72]+alpha_vdim[8]*f[70]+alpha_cdim[5]*f[67]+alpha_vdim[3]*f[66])+0.4841229182759271*(alpha_vdim[0]*f[46]+alpha_vdim[2]*f[31]+alpha_vdim[3]*f[30]+alpha_vdim[9]*f[16]); 
  out[84] += 0.2165063509461096*(alpha_cdim[5]*f[87]+alpha_cdim[0]*(f[77]+f[76])+alpha_vdim[0]*f[70]+alpha_vdim[1]*f[66]+alpha_vdim[2]*f[65]+alpha_vdim[7]*f[64])+0.4841229182759271*(alpha_vdim[3]*f[58]+alpha_vdim[9]*f[49]+alpha_vdim[0]*f[48])+0.1936491673103708*alpha_cdim[4]*f[36]+0.4841229182759271*alpha_vdim[2]*f[35]; 
  out[85] += 0.2165063509461096*(alpha_vdim[2]*f[80]+alpha_cdim[0]*(f[78]+f[76])+alpha_vdim[7]*f[72]+alpha_vdim[0]*f[71]+alpha_vdim[1]*f[67]+alpha_cdim[6]*f[65])+0.4841229182759271*(alpha_vdim[2]*f[58]+alpha_vdim[0]*f[49]+alpha_vdim[9]*f[48])+0.1936491673103708*alpha_cdim[4]*f[37]+0.4841229182759271*alpha_vdim[3]*f[35]; 
  out[86] += 0.2165063509461096*(alpha_cdim[5]*f[89]+alpha_vdim[1]*f[80]+alpha_cdim[0]*(f[78]+f[77])+alpha_vdim[0]*f[72]+alpha_vdim[7]*f[71]+alpha_vdim[2]*f[67]+alpha_cdim[6]*f[66])+0.4841229182759271*(alpha_vdim[0]*f[50]+alpha_vdim[2]*f[37]+alpha_vdim[3]*f[36]+alpha_vdim[9]*f[20]); 
  out[87] += 0.2165063509461096*(alpha_vdim[3]*f[85]+alpha_vdim[2]*f[81]+alpha_cdim[0]*f[79]+alpha_vdim[8]*f[78]+alpha_vdim[0]*f[76]+alpha_vdim[7]*f[74]+alpha_vdim[0]*f[73]+alpha_vdim[1]*(f[69]+f[68]))+0.4841229182759271*(alpha_vdim[9]*f[63]+alpha_vdim[3]*f[61]+alpha_vdim[2]*f[60]+alpha_vdim[0]*f[54])+0.1936491673103708*alpha_cdim[4]*f[41]; 
  out[88] += 0.2165063509461096*(alpha_vdim[8]*f[91]+alpha_vdim[3]*f[86]+alpha_vdim[1]*(f[84]+f[81])+alpha_cdim[0]*f[79]+alpha_vdim[0]*(f[77]+f[74])+alpha_vdim[7]*f[73]+alpha_cdim[5]*f[69]+alpha_vdim[2]*f[68])+0.4841229182759271*(alpha_vdim[3]*f[62]+alpha_vdim[9]*f[56]+alpha_vdim[0]*f[55]+alpha_vdim[2]*f[41]); 
  out[89] += 0.2165063509461096*(alpha_vdim[7]*f[90]+alpha_vdim[1]*f[85]+alpha_vdim[2]*f[83]+alpha_vdim[1]*f[82]+alpha_cdim[0]*f[79]+alpha_vdim[0]*f[78]+alpha_vdim[8]*f[76]+alpha_vdim[0]*f[75]+alpha_vdim[3]*f[69]+alpha_cdim[6]*f[68])+0.4841229182759271*(alpha_vdim[2]*f[62]+alpha_vdim[0]*f[56]+alpha_vdim[9]*f[55]+alpha_vdim[3]*f[41]); 
  out[90] += 0.2165063509461096*(alpha_cdim[6]*f[92]+alpha_cdim[0]*(f[83]+f[82]+f[81])+alpha_vdim[0]*f[80]+alpha_vdim[1]*f[72]+alpha_cdim[5]*f[71]+alpha_vdim[3]*f[70]+alpha_vdim[8]*f[66])+0.4841229182759271*alpha_vdim[0]*f[57]+0.1936491673103708*alpha_cdim[4]*f[46]+0.4841229182759271*(alpha_vdim[2]*f[45]+alpha_vdim[3]*f[44]+alpha_vdim[9]*f[29]); 
  out[91] += 0.2165063509461096*(alpha_cdim[5]*f[93]+alpha_cdim[0]*(f[86]+f[85]+f[84])+alpha_vdim[0]*f[80]+alpha_vdim[1]*f[72]+alpha_vdim[2]*f[71]+alpha_cdim[6]*f[70]+alpha_vdim[7]*f[67])+0.4841229182759271*alpha_vdim[0]*f[58]+0.1936491673103708*alpha_cdim[4]*f[50]+0.4841229182759271*(alpha_vdim[2]*f[49]+alpha_vdim[3]*f[48]+alpha_vdim[9]*f[35]); 
  out[92] += 0.2165063509461096*(alpha_vdim[3]*f[91]+alpha_cdim[0]*(f[88]+f[87])+alpha_vdim[8]*f[86]+alpha_vdim[0]*(f[84]+f[81])+alpha_vdim[1]*f[77]+alpha_cdim[5]*f[76]+alpha_vdim[1]*f[74]+alpha_vdim[2]*f[73]+alpha_vdim[7]*f[68])+0.4841229182759271*(alpha_vdim[3]*f[63]+alpha_vdim[9]*f[61]+alpha_vdim[0]*f[60])+0.1936491673103708*alpha_cdim[4]*f[55]+0.4841229182759271*alpha_vdim[2]*f[54]; 
  out[93] += 0.2165063509461096*(alpha_vdim[2]*f[90]+alpha_cdim[0]*(f[89]+f[87])+alpha_vdim[0]*f[85]+alpha_vdim[7]*f[83]+alpha_vdim[0]*f[82]+alpha_vdim[1]*f[78]+alpha_vdim[3]*f[76]+alpha_vdim[1]*f[75]+alpha_cdim[6]*f[73]+alpha_vdim[8]*f[69])+0.4841229182759271*(alpha_vdim[2]*f[63]+alpha_vdim[0]*f[61]+alpha_vdim[9]*f[60])+0.1936491673103708*alpha_cdim[4]*f[56]+0.4841229182759271*alpha_vdim[3]*f[54]; 
  out[94] += 0.2165063509461096*(alpha_vdim[1]*(f[91]+f[90])+alpha_cdim[0]*(f[89]+f[88])+alpha_vdim[0]*f[86]+alpha_vdim[8]*f[84]+alpha_vdim[0]*f[83]+alpha_vdim[7]*f[82]+alpha_cdim[5]*f[78]+alpha_vdim[3]*f[77]+alpha_vdim[2]*f[75]+alpha_cdim[6]*f[74])+0.4841229182759271*(alpha_vdim[0]*f[62]+alpha_vdim[2]*f[56]+alpha_vdim[3]*f[55]+alpha_vdim[9]*f[41]); 
  out[95] += 0.2165063509461096*(alpha_cdim[0]*(f[94]+f[93]+f[92])+alpha_vdim[0]*(f[91]+f[90])+alpha_vdim[1]*f[86]+alpha_cdim[5]*f[85]+alpha_vdim[3]*f[84]+alpha_vdim[1]*f[83]+alpha_vdim[2]*f[82]+alpha_cdim[6]*f[81]+alpha_vdim[8]*f[77]+alpha_vdim[7]*f[75])+0.4841229182759271*alpha_vdim[0]*f[63]+0.1936491673103708*alpha_cdim[4]*f[62]+0.4841229182759271*(alpha_vdim[2]*f[61]+alpha_vdim[3]*f[60]+alpha_vdim[9]*f[54]); 
  out[96] += 0.4841229182759271*(alpha_vdim[8]*f[27]+alpha_vdim[3]*f[15]+alpha_vdim[1]*f[13]+alpha_vdim[0]*f[5]); 
  out[97] += 0.2165063509461096*(alpha_cdim[4]*f[100]+alpha_cdim[0]*f[96])+0.4841229182759271*(alpha_vdim[3]*f[27]+alpha_vdim[8]*f[15]+alpha_vdim[0]*f[13]+alpha_vdim[1]*f[5]); 
  out[98] += 0.2165063509461096*alpha_cdim[0]*f[96]+0.4841229182759271*(alpha_vdim[8]*f[43]+alpha_vdim[3]*f[28]+alpha_vdim[1]*f[26]+alpha_vdim[0]*f[14])+0.1936491673103708*alpha_cdim[5]*f[5]; 
  out[99] += 0.2165063509461096*(alpha_cdim[6]*f[101]+alpha_cdim[0]*f[96])+0.4841229182759271*(alpha_vdim[1]*f[27]+alpha_vdim[0]*f[15]+alpha_vdim[8]*f[13]+alpha_vdim[3]*f[5]); 
  out[100] += 0.2165063509461096*(alpha_vdim[9]*f[104]+alpha_vdim[3]*f[99]+alpha_vdim[2]*f[98]+alpha_vdim[0]*f[96])+0.4841229182759271*(alpha_vdim[8]*f[45]+alpha_vdim[3]*f[31]+alpha_vdim[1]*f[29]+alpha_vdim[0]*f[16]); 
  out[101] += 0.2165063509461096*(alpha_vdim[7]*f[102]+alpha_vdim[2]*f[98]+alpha_vdim[1]*f[97]+alpha_vdim[0]*f[96])+0.4841229182759271*(alpha_vdim[8]*f[52]+alpha_vdim[3]*f[40]+alpha_vdim[1]*f[38]+alpha_vdim[0]*f[21]); 
  out[102] += 0.2165063509461096*(alpha_cdim[4]*f[106]+alpha_cdim[0]*(f[98]+f[97]))+0.4841229182759271*(alpha_vdim[3]*f[43]+alpha_vdim[8]*f[28]+alpha_vdim[0]*f[26]+alpha_vdim[1]*f[14])+0.1936491673103708*alpha_cdim[5]*f[13]; 
  out[103] += 0.2165063509461096*(alpha_cdim[6]*f[108]+alpha_cdim[4]*f[107]+alpha_cdim[0]*(f[99]+f[97]))+0.4841229182759271*(alpha_vdim[0]*f[27]+alpha_vdim[1]*f[15]+alpha_vdim[3]*f[13]+f[5]*alpha_vdim[8]); 
  out[104] += 0.2165063509461096*(alpha_cdim[6]*f[109]+alpha_cdim[0]*(f[99]+f[98]))+0.4841229182759271*(alpha_vdim[1]*f[43]+alpha_vdim[0]*f[28]+alpha_vdim[8]*f[26])+0.1936491673103708*alpha_cdim[5]*f[15]+0.4841229182759271*alpha_vdim[3]*f[14]; 
  out[105] += 0.2165063509461096*(alpha_vdim[9]*f[112]+alpha_vdim[3]*f[103]+alpha_vdim[2]*f[102]+alpha_cdim[0]*f[100]+alpha_vdim[0]*f[97]+alpha_cdim[4]*f[96])+0.4841229182759271*(alpha_vdim[3]*f[45]+alpha_vdim[8]*f[31]+alpha_vdim[0]*f[29]+alpha_vdim[1]*f[16]); 
  out[106] += 0.2165063509461096*(alpha_vdim[3]*f[104]+alpha_cdim[0]*f[100]+alpha_vdim[9]*f[99]+alpha_vdim[0]*f[98]+alpha_vdim[2]*f[96])+0.4841229182759271*(alpha_vdim[8]*f[57]+alpha_vdim[3]*f[46]+alpha_vdim[1]*f[44]+alpha_vdim[0]*f[30])+0.1936491673103708*alpha_cdim[5]*f[16]; 
  out[107] += 0.2165063509461096*(alpha_cdim[6]*f[111]+alpha_vdim[2]*f[104]+alpha_cdim[0]*f[100]+alpha_vdim[0]*f[99]+alpha_vdim[9]*f[98]+alpha_vdim[3]*f[96])+0.4841229182759271*(alpha_vdim[1]*f[45]+alpha_vdim[0]*f[31]+alpha_vdim[8]*f[29]+alpha_vdim[3]*f[16]); 
  out[108] += 0.2165063509461096*(alpha_cdim[4]*f[111]+alpha_vdim[2]*f[102]+alpha_cdim[0]*f[101]+alpha_vdim[7]*f[98]+alpha_vdim[0]*f[97]+alpha_vdim[1]*f[96])+0.4841229182759271*(alpha_vdim[3]*f[52]+alpha_vdim[8]*f[40]+alpha_vdim[0]*f[38]+alpha_vdim[1]*f[21]); 
  out[109] += 0.2165063509461096*(alpha_vdim[1]*f[102]+alpha_cdim[0]*f[101]+alpha_vdim[0]*f[98]+alpha_vdim[7]*f[97]+alpha_vdim[2]*f[96])+0.4841229182759271*(alpha_vdim[8]*f[59]+alpha_vdim[3]*f[53]+alpha_vdim[1]*f[51]+alpha_vdim[0]*f[39])+0.1936491673103708*alpha_cdim[5]*f[21]; 
  out[110] += 0.2165063509461096*(alpha_vdim[7]*f[112]+alpha_vdim[2]*f[104]+alpha_vdim[1]*f[103]+alpha_cdim[0]*f[101]+alpha_vdim[0]*f[99]+alpha_cdim[6]*f[96])+0.4841229182759271*(alpha_vdim[1]*f[52]+alpha_vdim[0]*f[40]+alpha_vdim[8]*f[38]+alpha_vdim[3]*f[21]); 
  out[111] += 0.2165063509461096*(alpha_vdim[9]*f[118]+alpha_vdim[7]*f[113]+alpha_vdim[3]*f[110]+alpha_vdim[2]*(f[109]+f[106])+alpha_vdim[1]*f[105]+alpha_vdim[0]*(f[101]+f[100]))+0.4841229182759271*(alpha_vdim[8]*f[61]+alpha_vdim[3]*f[56]+alpha_vdim[1]*f[54]+alpha_vdim[0]*f[41]); 
  out[112] += 0.2165063509461096*(alpha_cdim[6]*f[116]+alpha_cdim[4]*f[115]+alpha_cdim[0]*(f[104]+f[103]+f[102]))+0.4841229182759271*(alpha_vdim[0]*f[43]+alpha_vdim[1]*f[28])+0.1936491673103708*alpha_cdim[5]*f[27]+0.4841229182759271*(alpha_vdim[3]*f[26]+alpha_vdim[8]*f[14]); 
  out[113] += 0.2165063509461096*(alpha_vdim[3]*f[112]+alpha_cdim[0]*(f[106]+f[105])+alpha_vdim[9]*f[103]+alpha_vdim[0]*f[102]+alpha_cdim[4]*f[98]+alpha_vdim[2]*f[97])+0.4841229182759271*(alpha_vdim[3]*f[57]+alpha_vdim[8]*f[46]+alpha_vdim[0]*f[44]+alpha_vdim[1]*f[30])+0.1936491673103708*alpha_cdim[5]*f[29]; 
  out[114] += 0.2165063509461096*(alpha_cdim[6]*f[119]+alpha_vdim[2]*f[112]+alpha_cdim[0]*(f[107]+f[105])+alpha_vdim[0]*f[103]+alpha_vdim[9]*f[102]+alpha_cdim[4]*f[99]+alpha_vdim[3]*f[97])+0.4841229182759271*(alpha_vdim[0]*f[45]+alpha_vdim[1]*f[31]+alpha_vdim[3]*f[29]+alpha_vdim[8]*f[16]); 
  out[115] += 0.2165063509461096*(alpha_cdim[6]*f[120]+alpha_cdim[0]*(f[107]+f[106])+alpha_vdim[0]*f[104]+alpha_vdim[2]*f[99]+alpha_vdim[3]*f[98]+alpha_vdim[9]*f[96])+0.4841229182759271*(alpha_vdim[1]*f[57]+alpha_vdim[0]*f[46]+alpha_vdim[8]*f[44])+0.1936491673103708*alpha_cdim[5]*f[31]+0.4841229182759271*alpha_vdim[3]*f[30]; 
  out[116] += 0.2165063509461096*(alpha_cdim[4]*f[120]+alpha_cdim[0]*(f[109]+f[108])+alpha_vdim[0]*f[102]+alpha_vdim[1]*f[98]+alpha_vdim[2]*f[97]+alpha_vdim[7]*f[96])+0.4841229182759271*(alpha_vdim[3]*f[59]+alpha_vdim[8]*f[53]+alpha_vdim[0]*f[51]+alpha_vdim[1]*f[39])+0.1936491673103708*alpha_cdim[5]*f[38]; 
  out[117] += 0.2165063509461096*(alpha_cdim[4]*f[121]+alpha_vdim[2]*f[112]+alpha_cdim[0]*(f[110]+f[108])+alpha_vdim[7]*f[104]+alpha_vdim[0]*f[103]+alpha_vdim[1]*f[99]+alpha_cdim[6]*f[97])+0.4841229182759271*(alpha_vdim[0]*f[52]+alpha_vdim[1]*f[40]+alpha_vdim[3]*f[38]+alpha_vdim[8]*f[21]); 
  out[118] += 0.2165063509461096*(alpha_vdim[1]*f[112]+alpha_cdim[0]*(f[110]+f[109])+alpha_vdim[0]*f[104]+alpha_vdim[7]*f[103]+alpha_vdim[2]*f[99]+alpha_cdim[6]*f[98])+0.4841229182759271*(alpha_vdim[1]*f[59]+alpha_vdim[0]*f[53]+alpha_vdim[8]*f[51])+0.1936491673103708*alpha_cdim[5]*f[40]+0.4841229182759271*alpha_vdim[3]*f[39]; 
  out[119] += 0.2165063509461096*(alpha_vdim[9]*f[123]+alpha_vdim[3]*f[117]+alpha_vdim[2]*(f[116]+f[113])+alpha_cdim[0]*f[111]+alpha_vdim[0]*f[108]+alpha_vdim[7]*f[106]+alpha_vdim[0]*f[105]+alpha_cdim[4]*f[101]+alpha_vdim[1]*f[100])+0.4841229182759271*(alpha_vdim[3]*f[61]+alpha_vdim[8]*f[56]+alpha_vdim[0]*f[54]+alpha_vdim[1]*f[41]); 
  out[120] += 0.2165063509461096*(alpha_vdim[3]*f[118]+alpha_vdim[1]*f[113]+alpha_cdim[0]*f[111]+alpha_vdim[9]*f[110]+alpha_vdim[0]*(f[109]+f[106])+alpha_vdim[7]*f[105]+alpha_vdim[2]*(f[101]+f[100]))+0.4841229182759271*(alpha_vdim[8]*f[63]+alpha_vdim[3]*f[62]+alpha_vdim[1]*f[60]+alpha_vdim[0]*f[55])+0.1936491673103708*alpha_cdim[5]*f[41]; 
  out[121] += 0.2165063509461096*(alpha_vdim[7]*f[122]+alpha_vdim[2]*(f[118]+f[115])+alpha_vdim[1]*f[114]+alpha_cdim[0]*f[111]+alpha_vdim[0]*f[110]+alpha_vdim[9]*f[109]+alpha_vdim[0]*f[107]+alpha_vdim[3]*f[101]+alpha_cdim[6]*f[100])+0.4841229182759271*(alpha_vdim[1]*f[61]+alpha_vdim[0]*f[56]+alpha_vdim[8]*f[54]+alpha_vdim[3]*f[41]); 
  out[122] += 0.2165063509461096*(alpha_cdim[6]*f[124]+alpha_cdim[0]*(f[115]+f[114]+f[113])+alpha_vdim[0]*f[112]+alpha_cdim[4]*f[104]+alpha_vdim[2]*f[103]+alpha_vdim[3]*f[102]+alpha_vdim[9]*f[97])+0.4841229182759271*(alpha_vdim[0]*f[57]+alpha_vdim[1]*f[46])+0.1936491673103708*alpha_cdim[5]*f[45]+0.4841229182759271*(alpha_vdim[3]*f[44]+alpha_vdim[8]*f[30]); 
  out[123] += 0.2165063509461096*(alpha_cdim[4]*f[126]+alpha_cdim[0]*(f[118]+f[117]+f[116])+alpha_vdim[0]*f[112]+alpha_vdim[1]*f[104]+alpha_vdim[2]*f[103]+alpha_cdim[6]*f[102]+alpha_vdim[7]*f[99])+0.4841229182759271*(alpha_vdim[0]*f[59]+alpha_vdim[1]*f[53])+0.1936491673103708*alpha_cdim[5]*f[52]+0.4841229182759271*(alpha_vdim[3]*f[51]+alpha_vdim[8]*f[39]); 
  out[124] += 0.2165063509461096*(alpha_vdim[3]*f[123]+alpha_cdim[0]*(f[120]+f[119])+alpha_vdim[9]*f[117]+alpha_vdim[0]*(f[116]+f[113])+alpha_cdim[4]*f[109]+alpha_vdim[2]*f[108]+alpha_vdim[1]*f[106]+alpha_vdim[2]*f[105]+alpha_vdim[7]*f[100])+0.4841229182759271*(alpha_vdim[3]*f[63]+alpha_vdim[8]*f[62]+alpha_vdim[0]*f[60]+alpha_vdim[1]*f[55])+0.1936491673103708*alpha_cdim[5]*f[54]; 
  out[125] += 0.2165063509461096*(alpha_vdim[2]*(f[123]+f[122])+alpha_cdim[0]*(f[121]+f[119])+alpha_vdim[0]*f[117]+alpha_vdim[9]*f[116]+alpha_vdim[7]*f[115]+alpha_vdim[0]*f[114]+alpha_cdim[4]*f[110]+alpha_vdim[3]*f[108]+alpha_vdim[1]*f[107]+alpha_cdim[6]*f[105])+0.4841229182759271*(alpha_vdim[0]*f[61]+alpha_vdim[1]*f[56]+alpha_vdim[3]*f[54]+alpha_vdim[8]*f[41]); 
  out[126] += 0.2165063509461096*(alpha_vdim[1]*f[122]+alpha_cdim[0]*(f[121]+f[120])+alpha_vdim[0]*(f[118]+f[115])+alpha_vdim[7]*f[114]+alpha_vdim[2]*f[110]+alpha_vdim[3]*f[109]+alpha_vdim[2]*f[107]+alpha_cdim[6]*f[106]+alpha_vdim[9]*f[101])+0.4841229182759271*(alpha_vdim[1]*f[63]+alpha_vdim[0]*f[62]+alpha_vdim[8]*f[60])+0.1936491673103708*alpha_cdim[5]*f[56]+0.4841229182759271*alpha_vdim[3]*f[55]; 
  out[127] += 0.2165063509461096*(alpha_cdim[0]*(f[126]+f[125]+f[124])+alpha_vdim[0]*(f[123]+f[122])+alpha_cdim[4]*f[118]+alpha_vdim[2]*f[117]+alpha_vdim[3]*f[116]+alpha_vdim[1]*f[115]+alpha_vdim[2]*f[114]+alpha_cdim[6]*f[113]+alpha_vdim[9]*f[108]+alpha_vdim[7]*f[107])+0.4841229182759271*(alpha_vdim[0]*f[63]+alpha_vdim[1]*f[62])+0.1936491673103708*alpha_cdim[5]*f[61]+0.4841229182759271*(alpha_vdim[3]*f[60]+alpha_vdim[8]*f[55]); 
  out[128] += 0.4841229182759271*(alpha_vdim[7]*f[32]+alpha_vdim[2]*f[18]+alpha_vdim[1]*f[17]+alpha_vdim[0]*f[6]); 
  out[129] += 0.2165063509461096*(alpha_cdim[4]*f[132]+alpha_cdim[0]*f[128])+0.4841229182759271*(alpha_vdim[2]*f[32]+alpha_vdim[7]*f[18]+alpha_vdim[0]*f[17]+alpha_vdim[1]*f[6]); 
  out[130] += 0.2165063509461096*(alpha_cdim[5]*f[133]+alpha_cdim[0]*f[128])+0.4841229182759271*(alpha_vdim[1]*f[32]+alpha_vdim[0]*f[18]+alpha_vdim[7]*f[17]+alpha_vdim[2]*f[6]); 
  out[131] += 0.2165063509461096*alpha_cdim[0]*f[128]+0.4841229182759271*(alpha_vdim[7]*f[47]+alpha_vdim[2]*f[34]+alpha_vdim[1]*f[33]+alpha_vdim[0]*f[19])+0.1936491673103708*alpha_cdim[6]*f[6]; 
  out[132] += 0.2165063509461096*(alpha_vdim[9]*f[136]+alpha_vdim[3]*f[131]+alpha_vdim[2]*f[130]+alpha_vdim[0]*f[128])+0.4841229182759271*(alpha_vdim[7]*f[48]+alpha_vdim[2]*f[36]+alpha_vdim[1]*f[35]+alpha_vdim[0]*f[20]); 
  out[133] += 0.2165063509461096*(alpha_vdim[8]*f[135]+alpha_vdim[3]*f[131]+alpha_vdim[1]*f[129]+alpha_vdim[0]*f[128])+0.4841229182759271*(alpha_vdim[7]*f[51]+alpha_vdim[2]*f[39]+alpha_vdim[1]*f[38]+alpha_vdim[0]*f[21]); 
  out[134] += 0.2165063509461096*(alpha_cdim[5]*f[140]+alpha_cdim[4]*f[138]+alpha_cdim[0]*(f[130]+f[129]))+0.4841229182759271*(alpha_vdim[0]*f[32]+alpha_vdim[1]*f[18]+alpha_vdim[2]*f[17]+f[6]*alpha_vdim[7]); 
  out[135] += 0.2165063509461096*(alpha_cdim[4]*f[139]+alpha_cdim[0]*(f[131]+f[129]))+0.4841229182759271*(alpha_vdim[2]*f[47]+alpha_vdim[7]*f[34]+alpha_vdim[0]*f[33]+alpha_vdim[1]*f[19])+0.1936491673103708*alpha_cdim[6]*f[17]; 
  out[136] += 0.2165063509461096*(alpha_cdim[5]*f[142]+alpha_cdim[0]*(f[131]+f[130]))+0.4841229182759271*(alpha_vdim[1]*f[47]+alpha_vdim[0]*f[34]+alpha_vdim[7]*f[33]+alpha_vdim[2]*f[19])+0.1936491673103708*alpha_cdim[6]*f[18]; 
  out[137] += 0.2165063509461096*(alpha_vdim[9]*f[144]+alpha_vdim[3]*f[135]+alpha_vdim[2]*f[134]+alpha_cdim[0]*f[132]+alpha_vdim[0]*f[129]+alpha_cdim[4]*f[128])+0.4841229182759271*(alpha_vdim[2]*f[48]+alpha_vdim[7]*f[36]+alpha_vdim[0]*f[35]+alpha_vdim[1]*f[20]); 
  out[138] += 0.2165063509461096*(alpha_cdim[5]*f[143]+alpha_vdim[3]*f[136]+alpha_cdim[0]*f[132]+alpha_vdim[9]*f[131]+alpha_vdim[0]*f[130]+alpha_vdim[2]*f[128])+0.4841229182759271*(alpha_vdim[1]*f[48]+alpha_vdim[0]*f[36]+alpha_vdim[7]*f[35]+alpha_vdim[2]*f[20]); 
  out[139] += 0.2165063509461096*(alpha_vdim[2]*f[136]+alpha_cdim[0]*f[132]+alpha_vdim[0]*f[131]+alpha_vdim[9]*f[130]+alpha_vdim[3]*f[128])+0.4841229182759271*(alpha_vdim[7]*f[58]+alpha_vdim[2]*f[50]+alpha_vdim[1]*f[49]+alpha_vdim[0]*f[37])+0.1936491673103708*alpha_cdim[6]*f[20]; 
  out[140] += 0.2165063509461096*(alpha_cdim[4]*f[143]+alpha_vdim[3]*f[135]+alpha_cdim[0]*f[133]+alpha_vdim[8]*f[131]+alpha_vdim[0]*f[129]+alpha_vdim[1]*f[128])+0.4841229182759271*(alpha_vdim[2]*f[51]+alpha_vdim[7]*f[39]+alpha_vdim[0]*f[38]+alpha_vdim[1]*f[21]); 
  out[141] += 0.2165063509461096*(alpha_vdim[8]*f[144]+alpha_vdim[3]*f[136]+alpha_vdim[1]*f[134]+alpha_cdim[0]*f[133]+alpha_vdim[0]*f[130]+alpha_cdim[5]*f[128])+0.4841229182759271*(alpha_vdim[1]*f[51]+alpha_vdim[0]*f[39]+alpha_vdim[7]*f[38]+alpha_vdim[2]*f[21]); 
  out[142] += 0.2165063509461096*(alpha_vdim[1]*f[135]+alpha_cdim[0]*f[133]+alpha_vdim[0]*f[131]+alpha_vdim[8]*f[129]+alpha_vdim[3]*f[128])+0.4841229182759271*(alpha_vdim[7]*f[59]+alpha_vdim[2]*f[53]+alpha_vdim[1]*f[52]+alpha_vdim[0]*f[40])+0.1936491673103708*alpha_cdim[6]*f[21]; 
  out[143] += 0.2165063509461096*(alpha_vdim[9]*f[150]+alpha_vdim[8]*f[146]+alpha_vdim[3]*f[142]+alpha_vdim[2]*f[141]+alpha_vdim[3]*f[139]+alpha_vdim[1]*f[137]+alpha_vdim[0]*(f[133]+f[132]))+0.4841229182759271*(alpha_vdim[7]*f[60]+alpha_vdim[2]*f[55]+alpha_vdim[1]*f[54]+alpha_vdim[0]*f[41]); 
  out[144] += 0.2165063509461096*(alpha_cdim[5]*f[149]+alpha_cdim[4]*f[147]+alpha_cdim[0]*(f[136]+f[135]+f[134]))+0.4841229182759271*(alpha_vdim[0]*f[47]+alpha_vdim[1]*f[34]+alpha_vdim[2]*f[33])+0.1936491673103708*alpha_cdim[6]*f[32]+0.4841229182759271*alpha_vdim[7]*f[19]; 
  out[145] += 0.2165063509461096*(alpha_cdim[5]*f[151]+alpha_vdim[3]*f[144]+alpha_cdim[0]*(f[138]+f[137])+alpha_vdim[9]*f[135]+alpha_vdim[0]*f[134]+alpha_cdim[4]*f[130]+alpha_vdim[2]*f[129])+0.4841229182759271*(alpha_vdim[0]*f[48]+alpha_vdim[1]*f[36]+alpha_vdim[2]*f[35]+alpha_vdim[7]*f[20]); 
  out[146] += 0.2165063509461096*(alpha_vdim[2]*f[144]+alpha_cdim[0]*(f[139]+f[137])+alpha_vdim[0]*f[135]+alpha_vdim[9]*f[134]+alpha_cdim[4]*f[131]+alpha_vdim[3]*f[129])+0.4841229182759271*(alpha_vdim[2]*f[58]+alpha_vdim[7]*f[50]+alpha_vdim[0]*f[49]+alpha_vdim[1]*f[37])+0.1936491673103708*alpha_cdim[6]*f[35]; 
  out[147] += 0.2165063509461096*(alpha_cdim[5]*f[153]+alpha_cdim[0]*(f[139]+f[138])+alpha_vdim[0]*f[136]+alpha_vdim[2]*f[131]+alpha_vdim[3]*f[130]+alpha_vdim[9]*f[128])+0.4841229182759271*(alpha_vdim[1]*f[58]+alpha_vdim[0]*f[50]+alpha_vdim[7]*f[49]+alpha_vdim[2]*f[37])+0.1936491673103708*alpha_cdim[6]*f[36]; 
  out[148] += 0.2165063509461096*(alpha_cdim[4]*f[152]+alpha_vdim[3]*f[144]+alpha_cdim[0]*(f[141]+f[140])+alpha_vdim[8]*f[136]+alpha_vdim[0]*f[134]+alpha_vdim[1]*f[130]+alpha_cdim[5]*f[129])+0.4841229182759271*(alpha_vdim[0]*f[51]+alpha_vdim[1]*f[39]+alpha_vdim[2]*f[38]+alpha_vdim[7]*f[21]); 
  out[149] += 0.2165063509461096*(alpha_cdim[4]*f[153]+alpha_cdim[0]*(f[142]+f[140])+alpha_vdim[0]*f[135]+alpha_vdim[1]*f[131]+alpha_vdim[3]*f[129]+alpha_vdim[8]*f[128])+0.4841229182759271*(alpha_vdim[2]*f[59]+alpha_vdim[7]*f[53]+alpha_vdim[0]*f[52]+alpha_vdim[1]*f[40])+0.1936491673103708*alpha_cdim[6]*f[38]; 
  out[150] += 0.2165063509461096*(alpha_vdim[1]*f[144]+alpha_cdim[0]*(f[142]+f[141])+alpha_vdim[0]*f[136]+alpha_vdim[8]*f[134]+alpha_cdim[5]*f[131]+alpha_vdim[3]*f[130])+0.4841229182759271*(alpha_vdim[1]*f[59]+alpha_vdim[0]*f[53]+alpha_vdim[7]*f[52]+alpha_vdim[2]*f[40])+0.1936491673103708*alpha_cdim[6]*f[39]; 
  out[151] += 0.2165063509461096*(alpha_vdim[9]*f[155]+alpha_vdim[3]*f[149]+alpha_vdim[2]*f[148]+alpha_vdim[3]*f[146]+alpha_cdim[0]*f[143]+alpha_vdim[0]*f[140]+alpha_vdim[8]*f[139]+alpha_vdim[0]*f[137]+alpha_cdim[4]*f[133]+alpha_vdim[1]*f[132])+0.4841229182759271*(alpha_vdim[2]*f[60]+alpha_vdim[7]*f[55]+alpha_vdim[0]*f[54]+alpha_vdim[1]*f[41]); 
  out[152] += 0.2165063509461096*(alpha_vdim[8]*f[154]+alpha_vdim[3]*(f[150]+f[147])+alpha_vdim[1]*f[145]+alpha_cdim[0]*f[143]+alpha_vdim[9]*f[142]+alpha_vdim[0]*(f[141]+f[138])+alpha_vdim[2]*f[133]+alpha_cdim[5]*f[132])+0.4841229182759271*(alpha_vdim[1]*f[60]+alpha_vdim[0]*f[55]+alpha_vdim[7]*f[54]+alpha_vdim[2]*f[41]); 
  out[153] += 0.2165063509461096*(alpha_vdim[2]*f[150]+alpha_vdim[1]*f[146]+alpha_cdim[0]*f[143]+alpha_vdim[0]*f[142]+alpha_vdim[9]*f[141]+alpha_vdim[0]*f[139]+alpha_vdim[8]*f[137]+alpha_vdim[3]*(f[133]+f[132]))+0.4841229182759271*(alpha_vdim[7]*f[63]+alpha_vdim[2]*f[62]+alpha_vdim[1]*f[61]+alpha_vdim[0]*f[56])+0.1936491673103708*alpha_cdim[6]*f[41]; 
  out[154] += 0.2165063509461096*(alpha_cdim[5]*f[157]+alpha_cdim[0]*(f[147]+f[146]+f[145])+alpha_vdim[0]*f[144]+alpha_cdim[4]*f[136]+alpha_vdim[2]*f[135]+alpha_vdim[3]*f[134]+alpha_vdim[9]*f[129])+0.4841229182759271*(alpha_vdim[0]*f[58]+alpha_vdim[1]*f[50]+alpha_vdim[2]*f[49])+0.1936491673103708*alpha_cdim[6]*f[48]+0.4841229182759271*alpha_vdim[7]*f[37]; 
  out[155] += 0.2165063509461096*(alpha_cdim[4]*f[158]+alpha_cdim[0]*(f[150]+f[149]+f[148])+alpha_vdim[0]*f[144]+alpha_vdim[1]*f[136]+alpha_cdim[5]*f[135]+alpha_vdim[3]*f[134]+alpha_vdim[8]*f[130])+0.4841229182759271*(alpha_vdim[0]*f[59]+alpha_vdim[1]*f[53]+alpha_vdim[2]*f[52])+0.1936491673103708*alpha_cdim[6]*f[51]+0.4841229182759271*alpha_vdim[7]*f[40]; 
  out[156] += 0.2165063509461096*(alpha_vdim[3]*(f[155]+f[154])+alpha_cdim[0]*(f[152]+f[151])+alpha_vdim[9]*f[149]+alpha_vdim[0]*f[148]+alpha_vdim[8]*f[147]+alpha_vdim[0]*f[145]+alpha_cdim[4]*f[141]+alpha_vdim[2]*f[140]+alpha_vdim[1]*f[138]+alpha_cdim[5]*f[137])+0.4841229182759271*(alpha_vdim[0]*f[60]+alpha_vdim[1]*f[55]+alpha_vdim[2]*f[54]+alpha_vdim[7]*f[41]); 
  out[157] += 0.2165063509461096*(alpha_vdim[2]*f[155]+alpha_cdim[0]*(f[153]+f[151])+alpha_vdim[0]*f[149]+alpha_vdim[9]*f[148]+alpha_vdim[0]*f[146]+alpha_cdim[4]*f[142]+alpha_vdim[3]*f[140]+alpha_vdim[1]*f[139]+alpha_vdim[3]*f[137]+alpha_vdim[8]*f[132])+0.4841229182759271*(alpha_vdim[2]*f[63]+alpha_vdim[7]*f[62]+alpha_vdim[0]*f[61]+alpha_vdim[1]*f[56])+0.1936491673103708*alpha_cdim[6]*f[54]; 
  out[158] += 0.2165063509461096*(alpha_vdim[1]*f[154]+alpha_cdim[0]*(f[153]+f[152])+alpha_vdim[0]*(f[150]+f[147])+alpha_vdim[8]*f[145]+alpha_vdim[2]*f[142]+alpha_vdim[3]*f[141]+alpha_cdim[5]*f[139]+alpha_vdim[3]*f[138]+alpha_vdim[9]*f[133])+0.4841229182759271*(alpha_vdim[1]*f[63]+alpha_vdim[0]*f[62]+alpha_vdim[7]*f[61]+alpha_vdim[2]*f[56])+0.1936491673103708*alpha_cdim[6]*f[55]; 
  out[159] += 0.2165063509461096*(alpha_cdim[0]*(f[158]+f[157]+f[156])+alpha_vdim[0]*(f[155]+f[154])+alpha_cdim[4]*f[150]+alpha_vdim[2]*f[149]+alpha_vdim[3]*f[148]+alpha_vdim[1]*f[147]+alpha_cdim[5]*f[146]+alpha_vdim[3]*f[145]+alpha_vdim[9]*f[140]+alpha_vdim[8]*f[138])+0.4841229182759271*(alpha_vdim[0]*f[63]+alpha_vdim[1]*f[62]+alpha_vdim[2]*f[61])+0.1936491673103708*alpha_cdim[6]*f[60]+0.4841229182759271*alpha_vdim[7]*f[56]; 

  return cflFreq_mid; 
} 

