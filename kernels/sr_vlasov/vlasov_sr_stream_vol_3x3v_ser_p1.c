#include <gkyl_vlasov_sr_kernels.h> 
GKYL_CU_DH double vlasov_sr_stream_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // qmem:      q/m*EM fields (unused in streaming-only update).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double dx12 = 2/dxv[2]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[192] = {0.0}; 
  alpha_cdim[0] = 2.828427124746191*p0_over_gamma[0]*dx10; 
  alpha_cdim[4] = 2.828427124746191*p0_over_gamma[1]*dx10; 
  alpha_cdim[5] = 2.828427124746191*p0_over_gamma[2]*dx10; 
  alpha_cdim[6] = 2.828427124746191*p0_over_gamma[3]*dx10; 
  alpha_cdim[16] = 2.828427124746191*p0_over_gamma[4]*dx10; 
  alpha_cdim[20] = 2.828427124746191*p0_over_gamma[5]*dx10; 
  alpha_cdim[21] = 2.828427124746191*p0_over_gamma[6]*dx10; 
  alpha_cdim[41] = 2.828427124746191*p0_over_gamma[7]*dx10; 
  cflFreq_mid += fabs(0.0625*alpha_cdim[0]); 

  alpha_cdim[64] = 2.828427124746191*p1_over_gamma[0]*dx11; 
  alpha_cdim[68] = 2.828427124746191*p1_over_gamma[1]*dx11; 
  alpha_cdim[69] = 2.828427124746191*p1_over_gamma[2]*dx11; 
  alpha_cdim[70] = 2.828427124746191*p1_over_gamma[3]*dx11; 
  alpha_cdim[80] = 2.828427124746191*p1_over_gamma[4]*dx11; 
  alpha_cdim[84] = 2.828427124746191*p1_over_gamma[5]*dx11; 
  alpha_cdim[85] = 2.828427124746191*p1_over_gamma[6]*dx11; 
  alpha_cdim[105] = 2.828427124746191*p1_over_gamma[7]*dx11; 
  cflFreq_mid += fabs(0.0625*alpha_cdim[64]); 

  alpha_cdim[128] = 2.828427124746191*p2_over_gamma[0]*dx12; 
  alpha_cdim[132] = 2.828427124746191*p2_over_gamma[1]*dx12; 
  alpha_cdim[133] = 2.828427124746191*p2_over_gamma[2]*dx12; 
  alpha_cdim[134] = 2.828427124746191*p2_over_gamma[3]*dx12; 
  alpha_cdim[144] = 2.828427124746191*p2_over_gamma[4]*dx12; 
  alpha_cdim[148] = 2.828427124746191*p2_over_gamma[5]*dx12; 
  alpha_cdim[149] = 2.828427124746191*p2_over_gamma[6]*dx12; 
  alpha_cdim[169] = 2.828427124746191*p2_over_gamma[7]*dx12; 
  cflFreq_mid += fabs(0.0625*alpha_cdim[128]); 

  out[1] += 0.2165063509461096*(alpha_cdim[41]*f[41]+alpha_cdim[21]*f[21]+alpha_cdim[20]*f[20]+alpha_cdim[16]*f[16]+alpha_cdim[6]*f[6]+alpha_cdim[5]*f[5]+alpha_cdim[4]*f[4]+alpha_cdim[0]*f[0]); 
  out[2] += 0.2165063509461096*(f[41]*alpha_cdim[105]+f[21]*alpha_cdim[85]+f[20]*alpha_cdim[84]+f[16]*alpha_cdim[80]+f[6]*alpha_cdim[70]+f[5]*alpha_cdim[69]+f[4]*alpha_cdim[68]+f[0]*alpha_cdim[64]); 
  out[3] += 0.2165063509461096*(f[41]*alpha_cdim[169]+f[21]*alpha_cdim[149]+f[20]*alpha_cdim[148]+f[16]*alpha_cdim[144]+f[6]*alpha_cdim[134]+f[5]*alpha_cdim[133]+f[4]*alpha_cdim[132]+f[0]*alpha_cdim[128]); 
  out[7] += 0.2165063509461096*(f[54]*alpha_cdim[105]+f[38]*alpha_cdim[85]+f[35]*alpha_cdim[84]+f[29]*alpha_cdim[80]+f[17]*alpha_cdim[70]+f[13]*alpha_cdim[69]+f[10]*alpha_cdim[68]+f[1]*alpha_cdim[64]+alpha_cdim[41]*f[55]+alpha_cdim[21]*f[39]+alpha_cdim[20]*f[36]+alpha_cdim[16]*f[30]+alpha_cdim[6]*f[18]+alpha_cdim[5]*f[14]+alpha_cdim[4]*f[11]+alpha_cdim[0]*f[2]); 
  out[8] += 0.2165063509461096*(f[54]*alpha_cdim[169]+f[38]*alpha_cdim[149]+f[35]*alpha_cdim[148]+f[29]*alpha_cdim[144]+f[17]*alpha_cdim[134]+f[13]*alpha_cdim[133]+f[10]*alpha_cdim[132]+f[1]*alpha_cdim[128]+alpha_cdim[41]*f[56]+alpha_cdim[21]*f[40]+alpha_cdim[20]*f[37]+alpha_cdim[16]*f[31]+alpha_cdim[6]*f[19]+alpha_cdim[5]*f[15]+alpha_cdim[4]*f[12]+alpha_cdim[0]*f[3]); 
  out[9] += 0.2165063509461096*(f[55]*alpha_cdim[169]+f[39]*alpha_cdim[149]+f[36]*alpha_cdim[148]+f[30]*alpha_cdim[144]+f[18]*alpha_cdim[134]+f[14]*alpha_cdim[133]+f[11]*alpha_cdim[132]+f[2]*alpha_cdim[128]+f[56]*alpha_cdim[105]+f[40]*alpha_cdim[85]+f[37]*alpha_cdim[84]+f[31]*alpha_cdim[80]+f[19]*alpha_cdim[70]+f[15]*alpha_cdim[69]+f[12]*alpha_cdim[68]+f[3]*alpha_cdim[64]); 
  out[10] += 0.2165063509461096*(alpha_cdim[21]*f[41]+f[21]*alpha_cdim[41]+alpha_cdim[6]*f[20]+f[6]*alpha_cdim[20]+alpha_cdim[5]*f[16]+f[5]*alpha_cdim[16]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]); 
  out[11] += 0.2165063509461096*(f[21]*alpha_cdim[105]+f[41]*alpha_cdim[85]+f[6]*alpha_cdim[84]+f[5]*alpha_cdim[80]+f[20]*alpha_cdim[70]+f[16]*alpha_cdim[69]+f[0]*alpha_cdim[68]+f[4]*alpha_cdim[64]); 
  out[12] += 0.2165063509461096*(f[21]*alpha_cdim[169]+f[41]*alpha_cdim[149]+f[6]*alpha_cdim[148]+f[5]*alpha_cdim[144]+f[20]*alpha_cdim[134]+f[16]*alpha_cdim[133]+f[0]*alpha_cdim[132]+f[4]*alpha_cdim[128]); 
  out[13] += 0.2165063509461096*(alpha_cdim[20]*f[41]+f[20]*alpha_cdim[41]+alpha_cdim[6]*f[21]+f[6]*alpha_cdim[21]+alpha_cdim[4]*f[16]+f[4]*alpha_cdim[16]+alpha_cdim[0]*f[5]+f[0]*alpha_cdim[5]); 
  out[14] += 0.2165063509461096*(f[20]*alpha_cdim[105]+f[6]*alpha_cdim[85]+f[41]*alpha_cdim[84]+f[4]*alpha_cdim[80]+f[21]*alpha_cdim[70]+f[0]*alpha_cdim[69]+f[16]*alpha_cdim[68]+f[5]*alpha_cdim[64]); 
  out[15] += 0.2165063509461096*(f[20]*alpha_cdim[169]+f[6]*alpha_cdim[149]+f[41]*alpha_cdim[148]+f[4]*alpha_cdim[144]+f[21]*alpha_cdim[134]+f[0]*alpha_cdim[133]+f[16]*alpha_cdim[132]+f[5]*alpha_cdim[128]); 
  out[17] += 0.2165063509461096*(alpha_cdim[16]*f[41]+f[16]*alpha_cdim[41]+alpha_cdim[5]*f[21]+f[5]*alpha_cdim[21]+alpha_cdim[4]*f[20]+f[4]*alpha_cdim[20]+alpha_cdim[0]*f[6]+f[0]*alpha_cdim[6]); 
  out[18] += 0.2165063509461096*(f[16]*alpha_cdim[105]+f[5]*alpha_cdim[85]+f[4]*alpha_cdim[84]+f[41]*alpha_cdim[80]+f[0]*alpha_cdim[70]+f[21]*alpha_cdim[69]+f[20]*alpha_cdim[68]+f[6]*alpha_cdim[64]); 
  out[19] += 0.2165063509461096*(f[16]*alpha_cdim[169]+f[5]*alpha_cdim[149]+f[4]*alpha_cdim[148]+f[41]*alpha_cdim[144]+f[0]*alpha_cdim[134]+f[21]*alpha_cdim[133]+f[20]*alpha_cdim[132]+f[6]*alpha_cdim[128]); 
  out[22] += 0.2165063509461096*(f[60]*alpha_cdim[169]+f[51]*alpha_cdim[149]+f[48]*alpha_cdim[148]+f[44]*alpha_cdim[144]+f[32]*alpha_cdim[134]+f[26]*alpha_cdim[133]+f[23]*alpha_cdim[132]+f[7]*alpha_cdim[128]+f[61]*alpha_cdim[105]+f[52]*alpha_cdim[85]+f[49]*alpha_cdim[84]+f[45]*alpha_cdim[80]+f[33]*alpha_cdim[70]+f[27]*alpha_cdim[69]+f[24]*alpha_cdim[68]+f[8]*alpha_cdim[64]+alpha_cdim[41]*f[62]+alpha_cdim[21]*f[53]+alpha_cdim[20]*f[50]+alpha_cdim[16]*f[46]+alpha_cdim[6]*f[34]+alpha_cdim[5]*f[28]+alpha_cdim[4]*f[25]+alpha_cdim[0]*f[9]); 
  out[23] += 0.2165063509461096*(f[38]*alpha_cdim[105]+f[54]*alpha_cdim[85]+f[17]*alpha_cdim[84]+f[13]*alpha_cdim[80]+f[35]*alpha_cdim[70]+f[29]*alpha_cdim[69]+f[1]*alpha_cdim[68]+f[10]*alpha_cdim[64]+alpha_cdim[21]*f[55]+f[39]*alpha_cdim[41]+alpha_cdim[6]*f[36]+alpha_cdim[5]*f[30]+f[18]*alpha_cdim[20]+f[14]*alpha_cdim[16]+alpha_cdim[0]*f[11]+f[2]*alpha_cdim[4]); 
  out[24] += 0.2165063509461096*(f[38]*alpha_cdim[169]+f[54]*alpha_cdim[149]+f[17]*alpha_cdim[148]+f[13]*alpha_cdim[144]+f[35]*alpha_cdim[134]+f[29]*alpha_cdim[133]+f[1]*alpha_cdim[132]+f[10]*alpha_cdim[128]+alpha_cdim[21]*f[56]+f[40]*alpha_cdim[41]+alpha_cdim[6]*f[37]+alpha_cdim[5]*f[31]+f[19]*alpha_cdim[20]+f[15]*alpha_cdim[16]+alpha_cdim[0]*f[12]+f[3]*alpha_cdim[4]); 
  out[25] += 0.2165063509461096*(f[39]*alpha_cdim[169]+f[55]*alpha_cdim[149]+f[18]*alpha_cdim[148]+f[14]*alpha_cdim[144]+f[36]*alpha_cdim[134]+f[30]*alpha_cdim[133]+f[2]*alpha_cdim[132]+f[11]*alpha_cdim[128]+f[40]*alpha_cdim[105]+f[56]*alpha_cdim[85]+f[19]*alpha_cdim[84]+f[15]*alpha_cdim[80]+f[37]*alpha_cdim[70]+f[31]*alpha_cdim[69]+f[3]*alpha_cdim[68]+f[12]*alpha_cdim[64]); 
  out[26] += 0.2165063509461096*(f[35]*alpha_cdim[105]+f[17]*alpha_cdim[85]+f[54]*alpha_cdim[84]+f[10]*alpha_cdim[80]+f[38]*alpha_cdim[70]+f[1]*alpha_cdim[69]+f[29]*alpha_cdim[68]+f[13]*alpha_cdim[64]+alpha_cdim[20]*f[55]+f[36]*alpha_cdim[41]+alpha_cdim[6]*f[39]+alpha_cdim[4]*f[30]+f[18]*alpha_cdim[21]+f[11]*alpha_cdim[16]+alpha_cdim[0]*f[14]+f[2]*alpha_cdim[5]); 
  out[27] += 0.2165063509461096*(f[35]*alpha_cdim[169]+f[17]*alpha_cdim[149]+f[54]*alpha_cdim[148]+f[10]*alpha_cdim[144]+f[38]*alpha_cdim[134]+f[1]*alpha_cdim[133]+f[29]*alpha_cdim[132]+f[13]*alpha_cdim[128]+alpha_cdim[20]*f[56]+f[37]*alpha_cdim[41]+alpha_cdim[6]*f[40]+alpha_cdim[4]*f[31]+f[19]*alpha_cdim[21]+f[12]*alpha_cdim[16]+alpha_cdim[0]*f[15]+f[3]*alpha_cdim[5]); 
  out[28] += 0.2165063509461096*(f[36]*alpha_cdim[169]+f[18]*alpha_cdim[149]+f[55]*alpha_cdim[148]+f[11]*alpha_cdim[144]+f[39]*alpha_cdim[134]+f[2]*alpha_cdim[133]+f[30]*alpha_cdim[132]+f[14]*alpha_cdim[128]+f[37]*alpha_cdim[105]+f[19]*alpha_cdim[85]+f[56]*alpha_cdim[84]+f[12]*alpha_cdim[80]+f[40]*alpha_cdim[70]+f[3]*alpha_cdim[69]+f[31]*alpha_cdim[68]+f[15]*alpha_cdim[64]); 
  out[29] += 0.2165063509461096*(alpha_cdim[6]*f[41]+f[6]*alpha_cdim[41]+alpha_cdim[20]*f[21]+f[20]*alpha_cdim[21]+alpha_cdim[0]*f[16]+f[0]*alpha_cdim[16]+alpha_cdim[4]*f[5]+f[4]*alpha_cdim[5]); 
  out[30] += 0.2165063509461096*(f[6]*alpha_cdim[105]+f[20]*alpha_cdim[85]+f[21]*alpha_cdim[84]+f[0]*alpha_cdim[80]+f[41]*alpha_cdim[70]+f[4]*alpha_cdim[69]+f[5]*alpha_cdim[68]+f[16]*alpha_cdim[64]); 
  out[31] += 0.2165063509461096*(f[6]*alpha_cdim[169]+f[20]*alpha_cdim[149]+f[21]*alpha_cdim[148]+f[0]*alpha_cdim[144]+f[41]*alpha_cdim[134]+f[4]*alpha_cdim[133]+f[5]*alpha_cdim[132]+f[16]*alpha_cdim[128]); 
  out[32] += 0.2165063509461096*(f[29]*alpha_cdim[105]+f[13]*alpha_cdim[85]+f[10]*alpha_cdim[84]+f[54]*alpha_cdim[80]+f[1]*alpha_cdim[70]+f[38]*alpha_cdim[69]+f[35]*alpha_cdim[68]+f[17]*alpha_cdim[64]+alpha_cdim[16]*f[55]+f[30]*alpha_cdim[41]+alpha_cdim[5]*f[39]+alpha_cdim[4]*f[36]+f[14]*alpha_cdim[21]+f[11]*alpha_cdim[20]+alpha_cdim[0]*f[18]+f[2]*alpha_cdim[6]); 
  out[33] += 0.2165063509461096*(f[29]*alpha_cdim[169]+f[13]*alpha_cdim[149]+f[10]*alpha_cdim[148]+f[54]*alpha_cdim[144]+f[1]*alpha_cdim[134]+f[38]*alpha_cdim[133]+f[35]*alpha_cdim[132]+f[17]*alpha_cdim[128]+alpha_cdim[16]*f[56]+f[31]*alpha_cdim[41]+alpha_cdim[5]*f[40]+alpha_cdim[4]*f[37]+f[15]*alpha_cdim[21]+f[12]*alpha_cdim[20]+alpha_cdim[0]*f[19]+f[3]*alpha_cdim[6]); 
  out[34] += 0.2165063509461096*(f[30]*alpha_cdim[169]+f[14]*alpha_cdim[149]+f[11]*alpha_cdim[148]+f[55]*alpha_cdim[144]+f[2]*alpha_cdim[134]+f[39]*alpha_cdim[133]+f[36]*alpha_cdim[132]+f[18]*alpha_cdim[128]+f[31]*alpha_cdim[105]+f[15]*alpha_cdim[85]+f[12]*alpha_cdim[84]+f[56]*alpha_cdim[80]+f[3]*alpha_cdim[70]+f[40]*alpha_cdim[69]+f[37]*alpha_cdim[68]+f[19]*alpha_cdim[64]); 
  out[35] += 0.2165063509461096*(alpha_cdim[5]*f[41]+f[5]*alpha_cdim[41]+alpha_cdim[16]*f[21]+f[16]*alpha_cdim[21]+alpha_cdim[0]*f[20]+f[0]*alpha_cdim[20]+alpha_cdim[4]*f[6]+f[4]*alpha_cdim[6]); 
  out[36] += 0.2165063509461096*(f[5]*alpha_cdim[105]+f[16]*alpha_cdim[85]+f[0]*alpha_cdim[84]+f[21]*alpha_cdim[80]+f[4]*alpha_cdim[70]+f[41]*alpha_cdim[69]+f[6]*alpha_cdim[68]+f[20]*alpha_cdim[64]); 
  out[37] += 0.2165063509461096*(f[5]*alpha_cdim[169]+f[16]*alpha_cdim[149]+f[0]*alpha_cdim[148]+f[21]*alpha_cdim[144]+f[4]*alpha_cdim[134]+f[41]*alpha_cdim[133]+f[6]*alpha_cdim[132]+f[20]*alpha_cdim[128]); 
  out[38] += 0.2165063509461096*(alpha_cdim[4]*f[41]+f[4]*alpha_cdim[41]+alpha_cdim[0]*f[21]+f[0]*alpha_cdim[21]+alpha_cdim[16]*f[20]+f[16]*alpha_cdim[20]+alpha_cdim[5]*f[6]+f[5]*alpha_cdim[6]); 
  out[39] += 0.2165063509461096*(f[4]*alpha_cdim[105]+f[0]*alpha_cdim[85]+f[16]*alpha_cdim[84]+f[20]*alpha_cdim[80]+f[5]*alpha_cdim[70]+f[6]*alpha_cdim[69]+f[41]*alpha_cdim[68]+f[21]*alpha_cdim[64]); 
  out[40] += 0.2165063509461096*(f[4]*alpha_cdim[169]+f[0]*alpha_cdim[149]+f[16]*alpha_cdim[148]+f[20]*alpha_cdim[144]+f[5]*alpha_cdim[134]+f[6]*alpha_cdim[133]+f[41]*alpha_cdim[132]+f[21]*alpha_cdim[128]); 
  out[42] += 0.2165063509461096*(f[51]*alpha_cdim[169]+f[60]*alpha_cdim[149]+f[32]*alpha_cdim[148]+f[26]*alpha_cdim[144]+f[48]*alpha_cdim[134]+f[44]*alpha_cdim[133]+f[7]*alpha_cdim[132]+f[23]*alpha_cdim[128]+f[52]*alpha_cdim[105]+f[61]*alpha_cdim[85]+f[33]*alpha_cdim[84]+f[27]*alpha_cdim[80]+f[49]*alpha_cdim[70]+f[45]*alpha_cdim[69]+f[8]*alpha_cdim[68]+f[24]*alpha_cdim[64]+alpha_cdim[21]*f[62]+alpha_cdim[41]*f[53]+alpha_cdim[6]*f[50]+alpha_cdim[5]*f[46]+alpha_cdim[20]*f[34]+alpha_cdim[16]*f[28]+alpha_cdim[0]*f[25]+alpha_cdim[4]*f[9]); 
  out[43] += 0.2165063509461096*(f[48]*alpha_cdim[169]+f[32]*alpha_cdim[149]+f[60]*alpha_cdim[148]+f[23]*alpha_cdim[144]+f[51]*alpha_cdim[134]+f[7]*alpha_cdim[133]+f[44]*alpha_cdim[132]+f[26]*alpha_cdim[128]+f[49]*alpha_cdim[105]+f[33]*alpha_cdim[85]+f[61]*alpha_cdim[84]+f[24]*alpha_cdim[80]+f[52]*alpha_cdim[70]+f[8]*alpha_cdim[69]+f[45]*alpha_cdim[68]+f[27]*alpha_cdim[64]+alpha_cdim[20]*f[62]+alpha_cdim[6]*f[53]+alpha_cdim[41]*f[50]+alpha_cdim[4]*f[46]+alpha_cdim[21]*f[34]+alpha_cdim[0]*f[28]+alpha_cdim[16]*f[25]+alpha_cdim[5]*f[9]); 
  out[44] += 0.2165063509461096*(f[17]*alpha_cdim[105]+f[35]*alpha_cdim[85]+f[38]*alpha_cdim[84]+f[1]*alpha_cdim[80]+f[54]*alpha_cdim[70]+f[10]*alpha_cdim[69]+f[13]*alpha_cdim[68]+f[29]*alpha_cdim[64]+alpha_cdim[6]*f[55]+f[18]*alpha_cdim[41]+alpha_cdim[20]*f[39]+alpha_cdim[21]*f[36]+alpha_cdim[0]*f[30]+f[2]*alpha_cdim[16]+alpha_cdim[4]*f[14]+alpha_cdim[5]*f[11]); 
  out[45] += 0.2165063509461096*(f[17]*alpha_cdim[169]+f[35]*alpha_cdim[149]+f[38]*alpha_cdim[148]+f[1]*alpha_cdim[144]+f[54]*alpha_cdim[134]+f[10]*alpha_cdim[133]+f[13]*alpha_cdim[132]+f[29]*alpha_cdim[128]+alpha_cdim[6]*f[56]+f[19]*alpha_cdim[41]+alpha_cdim[20]*f[40]+alpha_cdim[21]*f[37]+alpha_cdim[0]*f[31]+f[3]*alpha_cdim[16]+alpha_cdim[4]*f[15]+alpha_cdim[5]*f[12]); 
  out[46] += 0.2165063509461096*(f[18]*alpha_cdim[169]+f[36]*alpha_cdim[149]+f[39]*alpha_cdim[148]+f[2]*alpha_cdim[144]+f[55]*alpha_cdim[134]+f[11]*alpha_cdim[133]+f[14]*alpha_cdim[132]+f[30]*alpha_cdim[128]+f[19]*alpha_cdim[105]+f[37]*alpha_cdim[85]+f[40]*alpha_cdim[84]+f[3]*alpha_cdim[80]+f[56]*alpha_cdim[70]+f[12]*alpha_cdim[69]+f[15]*alpha_cdim[68]+f[31]*alpha_cdim[64]); 
  out[47] += 0.2165063509461096*(f[44]*alpha_cdim[169]+f[26]*alpha_cdim[149]+f[23]*alpha_cdim[148]+f[60]*alpha_cdim[144]+f[7]*alpha_cdim[134]+f[51]*alpha_cdim[133]+f[48]*alpha_cdim[132]+f[32]*alpha_cdim[128]+f[45]*alpha_cdim[105]+f[27]*alpha_cdim[85]+f[24]*alpha_cdim[84]+f[61]*alpha_cdim[80]+f[8]*alpha_cdim[70]+f[52]*alpha_cdim[69]+f[49]*alpha_cdim[68]+f[33]*alpha_cdim[64]+alpha_cdim[16]*f[62]+alpha_cdim[5]*f[53]+alpha_cdim[4]*f[50]+alpha_cdim[41]*f[46]+alpha_cdim[0]*f[34]+alpha_cdim[21]*f[28]+alpha_cdim[20]*f[25]+alpha_cdim[6]*f[9]); 
  out[48] += 0.2165063509461096*(f[13]*alpha_cdim[105]+f[29]*alpha_cdim[85]+f[1]*alpha_cdim[84]+f[38]*alpha_cdim[80]+f[10]*alpha_cdim[70]+f[54]*alpha_cdim[69]+f[17]*alpha_cdim[68]+f[35]*alpha_cdim[64]+alpha_cdim[5]*f[55]+f[14]*alpha_cdim[41]+alpha_cdim[16]*f[39]+alpha_cdim[0]*f[36]+alpha_cdim[21]*f[30]+f[2]*alpha_cdim[20]+alpha_cdim[4]*f[18]+alpha_cdim[6]*f[11]); 
  out[49] += 0.2165063509461096*(f[13]*alpha_cdim[169]+f[29]*alpha_cdim[149]+f[1]*alpha_cdim[148]+f[38]*alpha_cdim[144]+f[10]*alpha_cdim[134]+f[54]*alpha_cdim[133]+f[17]*alpha_cdim[132]+f[35]*alpha_cdim[128]+alpha_cdim[5]*f[56]+f[15]*alpha_cdim[41]+alpha_cdim[16]*f[40]+alpha_cdim[0]*f[37]+alpha_cdim[21]*f[31]+f[3]*alpha_cdim[20]+alpha_cdim[4]*f[19]+alpha_cdim[6]*f[12]); 
  out[50] += 0.2165063509461096*(f[14]*alpha_cdim[169]+f[30]*alpha_cdim[149]+f[2]*alpha_cdim[148]+f[39]*alpha_cdim[144]+f[11]*alpha_cdim[134]+f[55]*alpha_cdim[133]+f[18]*alpha_cdim[132]+f[36]*alpha_cdim[128]+f[15]*alpha_cdim[105]+f[31]*alpha_cdim[85]+f[3]*alpha_cdim[84]+f[40]*alpha_cdim[80]+f[12]*alpha_cdim[70]+f[56]*alpha_cdim[69]+f[19]*alpha_cdim[68]+f[37]*alpha_cdim[64]); 
  out[51] += 0.2165063509461096*(f[10]*alpha_cdim[105]+f[1]*alpha_cdim[85]+f[29]*alpha_cdim[84]+f[35]*alpha_cdim[80]+f[13]*alpha_cdim[70]+f[17]*alpha_cdim[69]+f[54]*alpha_cdim[68]+f[38]*alpha_cdim[64]+alpha_cdim[4]*f[55]+f[11]*alpha_cdim[41]+alpha_cdim[0]*f[39]+alpha_cdim[16]*f[36]+alpha_cdim[20]*f[30]+f[2]*alpha_cdim[21]+alpha_cdim[5]*f[18]+alpha_cdim[6]*f[14]); 
  out[52] += 0.2165063509461096*(f[10]*alpha_cdim[169]+f[1]*alpha_cdim[149]+f[29]*alpha_cdim[148]+f[35]*alpha_cdim[144]+f[13]*alpha_cdim[134]+f[17]*alpha_cdim[133]+f[54]*alpha_cdim[132]+f[38]*alpha_cdim[128]+alpha_cdim[4]*f[56]+f[12]*alpha_cdim[41]+alpha_cdim[0]*f[40]+alpha_cdim[16]*f[37]+alpha_cdim[20]*f[31]+f[3]*alpha_cdim[21]+alpha_cdim[5]*f[19]+alpha_cdim[6]*f[15]); 
  out[53] += 0.2165063509461096*(f[11]*alpha_cdim[169]+f[2]*alpha_cdim[149]+f[30]*alpha_cdim[148]+f[36]*alpha_cdim[144]+f[14]*alpha_cdim[134]+f[18]*alpha_cdim[133]+f[55]*alpha_cdim[132]+f[39]*alpha_cdim[128]+f[12]*alpha_cdim[105]+f[3]*alpha_cdim[85]+f[31]*alpha_cdim[84]+f[37]*alpha_cdim[80]+f[15]*alpha_cdim[70]+f[19]*alpha_cdim[69]+f[56]*alpha_cdim[68]+f[40]*alpha_cdim[64]); 
  out[54] += 0.2165063509461096*(alpha_cdim[0]*f[41]+f[0]*alpha_cdim[41]+alpha_cdim[4]*f[21]+f[4]*alpha_cdim[21]+alpha_cdim[5]*f[20]+f[5]*alpha_cdim[20]+alpha_cdim[6]*f[16]+f[6]*alpha_cdim[16]); 
  out[55] += 0.2165063509461096*(f[0]*alpha_cdim[105]+f[4]*alpha_cdim[85]+f[5]*alpha_cdim[84]+f[6]*alpha_cdim[80]+f[16]*alpha_cdim[70]+f[20]*alpha_cdim[69]+f[21]*alpha_cdim[68]+f[41]*alpha_cdim[64]); 
  out[56] += 0.2165063509461096*(f[0]*alpha_cdim[169]+f[4]*alpha_cdim[149]+f[5]*alpha_cdim[148]+f[6]*alpha_cdim[144]+f[16]*alpha_cdim[134]+f[20]*alpha_cdim[133]+f[21]*alpha_cdim[132]+f[41]*alpha_cdim[128]); 
  out[57] += 0.2165063509461096*(f[32]*alpha_cdim[169]+f[48]*alpha_cdim[149]+f[51]*alpha_cdim[148]+f[7]*alpha_cdim[144]+f[60]*alpha_cdim[134]+f[23]*alpha_cdim[133]+f[26]*alpha_cdim[132]+f[44]*alpha_cdim[128]+f[33]*alpha_cdim[105]+f[49]*alpha_cdim[85]+f[52]*alpha_cdim[84]+f[8]*alpha_cdim[80]+f[61]*alpha_cdim[70]+f[24]*alpha_cdim[69]+f[27]*alpha_cdim[68]+f[45]*alpha_cdim[64]+alpha_cdim[6]*f[62]+alpha_cdim[20]*f[53]+alpha_cdim[21]*f[50]+alpha_cdim[0]*f[46]+f[34]*alpha_cdim[41]+alpha_cdim[4]*f[28]+alpha_cdim[5]*f[25]+f[9]*alpha_cdim[16]); 
  out[58] += 0.2165063509461096*(f[26]*alpha_cdim[169]+f[44]*alpha_cdim[149]+f[7]*alpha_cdim[148]+f[51]*alpha_cdim[144]+f[23]*alpha_cdim[134]+f[60]*alpha_cdim[133]+f[32]*alpha_cdim[132]+f[48]*alpha_cdim[128]+f[27]*alpha_cdim[105]+f[45]*alpha_cdim[85]+f[8]*alpha_cdim[84]+f[52]*alpha_cdim[80]+f[24]*alpha_cdim[70]+f[61]*alpha_cdim[69]+f[33]*alpha_cdim[68]+f[49]*alpha_cdim[64]+alpha_cdim[5]*f[62]+alpha_cdim[16]*f[53]+alpha_cdim[0]*f[50]+alpha_cdim[21]*f[46]+f[28]*alpha_cdim[41]+alpha_cdim[4]*f[34]+alpha_cdim[6]*f[25]+f[9]*alpha_cdim[20]); 
  out[59] += 0.2165063509461096*(f[23]*alpha_cdim[169]+f[7]*alpha_cdim[149]+f[44]*alpha_cdim[148]+f[48]*alpha_cdim[144]+f[26]*alpha_cdim[134]+f[32]*alpha_cdim[133]+f[60]*alpha_cdim[132]+f[51]*alpha_cdim[128]+f[24]*alpha_cdim[105]+f[8]*alpha_cdim[85]+f[45]*alpha_cdim[84]+f[49]*alpha_cdim[80]+f[27]*alpha_cdim[70]+f[33]*alpha_cdim[69]+f[61]*alpha_cdim[68]+f[52]*alpha_cdim[64]+alpha_cdim[4]*f[62]+alpha_cdim[0]*f[53]+alpha_cdim[16]*f[50]+alpha_cdim[20]*f[46]+f[25]*alpha_cdim[41]+alpha_cdim[5]*f[34]+alpha_cdim[6]*f[28]+f[9]*alpha_cdim[21]); 
  out[60] += 0.2165063509461096*(f[1]*alpha_cdim[105]+f[10]*alpha_cdim[85]+f[13]*alpha_cdim[84]+f[17]*alpha_cdim[80]+f[29]*alpha_cdim[70]+f[35]*alpha_cdim[69]+f[38]*alpha_cdim[68]+f[54]*alpha_cdim[64]+alpha_cdim[0]*f[55]+f[2]*alpha_cdim[41]+alpha_cdim[4]*f[39]+alpha_cdim[5]*f[36]+alpha_cdim[6]*f[30]+f[11]*alpha_cdim[21]+f[14]*alpha_cdim[20]+alpha_cdim[16]*f[18]); 
  out[61] += 0.2165063509461096*(f[1]*alpha_cdim[169]+f[10]*alpha_cdim[149]+f[13]*alpha_cdim[148]+f[17]*alpha_cdim[144]+f[29]*alpha_cdim[134]+f[35]*alpha_cdim[133]+f[38]*alpha_cdim[132]+f[54]*alpha_cdim[128]+alpha_cdim[0]*f[56]+f[3]*alpha_cdim[41]+alpha_cdim[4]*f[40]+alpha_cdim[5]*f[37]+alpha_cdim[6]*f[31]+f[12]*alpha_cdim[21]+f[15]*alpha_cdim[20]+alpha_cdim[16]*f[19]); 
  out[62] += 0.2165063509461096*(f[2]*alpha_cdim[169]+f[11]*alpha_cdim[149]+f[14]*alpha_cdim[148]+f[18]*alpha_cdim[144]+f[30]*alpha_cdim[134]+f[36]*alpha_cdim[133]+f[39]*alpha_cdim[132]+f[55]*alpha_cdim[128]+f[3]*alpha_cdim[105]+f[12]*alpha_cdim[85]+f[15]*alpha_cdim[84]+f[19]*alpha_cdim[80]+f[31]*alpha_cdim[70]+f[37]*alpha_cdim[69]+f[40]*alpha_cdim[68]+f[56]*alpha_cdim[64]); 
  out[63] += 0.2165063509461096*(f[7]*alpha_cdim[169]+f[23]*alpha_cdim[149]+f[26]*alpha_cdim[148]+f[32]*alpha_cdim[144]+f[44]*alpha_cdim[134]+f[48]*alpha_cdim[133]+f[51]*alpha_cdim[132]+f[60]*alpha_cdim[128]+f[8]*alpha_cdim[105]+f[24]*alpha_cdim[85]+f[27]*alpha_cdim[84]+f[33]*alpha_cdim[80]+f[45]*alpha_cdim[70]+f[49]*alpha_cdim[69]+f[52]*alpha_cdim[68]+f[61]*alpha_cdim[64]+alpha_cdim[0]*f[62]+alpha_cdim[4]*f[53]+alpha_cdim[5]*f[50]+alpha_cdim[6]*f[46]+f[9]*alpha_cdim[41]+alpha_cdim[16]*f[34]+alpha_cdim[20]*f[28]+alpha_cdim[21]*f[25]); 

  return 3.0*cflFreq_mid; 
} 
