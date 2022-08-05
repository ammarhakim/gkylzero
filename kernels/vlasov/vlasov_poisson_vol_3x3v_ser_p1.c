#include <gkyl_vlasov_kernels.h> 

GKYL_CU_DH double vlasov_poisson_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fac_phi:   potential (scaled by appropriate factors).
  // vecA:      vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[3]/dxv[0]; 
  double w0dx0 = w[3]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 
  double dv1dx1 = dxv[4]/dxv[1]; 
  double w1dx1 = w[4]/dxv[1]; 
  const double dx11 = 2/dxv[1]; 
  double dv2dx2 = dxv[5]/dxv[2]; 
  double w2dx2 = w[5]/dxv[2]; 
  const double dx12 = 2/dxv[2]; 
  const double *phi = &fac_phi[0]; 
  const double dv10 = 2/dxv[3]; 
  const double dv1 = dxv[3], wv1 = w[3]; 
  const double dv11 = 2/dxv[4]; 
  const double dv2 = dxv[4], wv2 = w[4]; 
  const double dv12 = 2/dxv[5]; 
  const double dv3 = dxv[5], wv3 = w[5]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[192]; 
  double alpha_vdim[192]; 

  alpha_cdim[0] = 16.0*w0dx0; 
  alpha_cdim[4] = 4.618802153517007*dv0dx0; 
  alpha_mid += fabs(w0dx0)+0.5*dv0dx0; 

  alpha_cdim[64] = 16.0*w1dx1; 
  alpha_cdim[69] = 4.618802153517007*dv1dx1; 
  alpha_mid += fabs(w1dx1)+0.5*dv1dx1; 

  alpha_cdim[128] = 16.0*w2dx2; 
  alpha_cdim[134] = 4.618802153517007*dv2dx2; 
  alpha_mid += fabs(w2dx2)+0.5*dv2dx2; 

  alpha_vdim[0] = -4.898979485566357*phi[1]*dv10*dx10; 
  alpha_vdim[2] = -4.898979485566357*phi[4]*dv10*dx10; 
  alpha_vdim[3] = -4.898979485566357*phi[5]*dv10*dx10; 
  alpha_vdim[9] = -4.898979485566357*phi[7]*dv10*dx10; 
  alpha_mid += fabs(0.0625*alpha_vdim[0]); 

  alpha_vdim[64] = -4.898979485566357*phi[2]*dv11*dx11; 
  alpha_vdim[65] = -4.898979485566357*phi[4]*dv11*dx11; 
  alpha_vdim[67] = -4.898979485566357*phi[6]*dv11*dx11; 
  alpha_vdim[72] = -4.898979485566357*phi[7]*dv11*dx11; 
  alpha_mid += fabs(0.0625*alpha_vdim[64]); 

  alpha_vdim[128] = -4.898979485566357*phi[3]*dv12*dx12; 
  alpha_vdim[129] = -4.898979485566357*phi[5]*dv12*dx12; 
  alpha_vdim[130] = -4.898979485566357*phi[6]*dv12*dx12; 
  alpha_vdim[135] = -4.898979485566357*phi[7]*dv12*dx12; 
  alpha_mid += fabs(0.0625*alpha_vdim[128]); 

  out[1] += 0.2165063509461096*(alpha_cdim[4]*f[4]+alpha_cdim[0]*f[0]); 
  out[2] += 0.2165063509461096*(f[5]*alpha_cdim[69]+f[0]*alpha_cdim[64]); 
  out[3] += 0.2165063509461096*(f[6]*alpha_cdim[134]+f[0]*alpha_cdim[128]); 
  out[4] += 0.2165063509461096*(alpha_vdim[9]*f[9]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[0]*f[0]); 
  out[5] += 0.2165063509461096*(f[8]*alpha_vdim[72]+f[3]*alpha_vdim[67]+f[1]*alpha_vdim[65]+f[0]*alpha_vdim[64]); 
  out[6] += 0.2165063509461096*(f[7]*alpha_vdim[135]+f[2]*alpha_vdim[130]+f[1]*alpha_vdim[129]+f[0]*alpha_vdim[128]); 
  out[7] += 0.2165063509461096*(f[13]*alpha_cdim[69]+f[1]*alpha_cdim[64]+alpha_cdim[4]*f[11]+alpha_cdim[0]*f[2]); 
  out[8] += 0.2165063509461096*(f[17]*alpha_cdim[134]+f[1]*alpha_cdim[128]+alpha_cdim[4]*f[12]+alpha_cdim[0]*f[3]); 
  out[9] += 0.2165063509461096*(f[18]*alpha_cdim[134]+f[2]*alpha_cdim[128]+f[15]*alpha_cdim[69]+f[3]*alpha_cdim[64]); 
  out[10] += 0.2165063509461096*(alpha_vdim[9]*f[22]+alpha_vdim[3]*f[8]+alpha_vdim[2]*f[7]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]+alpha_vdim[0]*f[1]); 
  out[11] += 0.2165063509461096*(f[16]*alpha_cdim[69]+f[4]*alpha_cdim[64]+alpha_vdim[3]*f[9]+f[3]*alpha_vdim[9]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[12] += 0.2165063509461096*(f[20]*alpha_cdim[134]+f[4]*alpha_cdim[128]+alpha_vdim[2]*f[9]+f[2]*alpha_vdim[9]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[13] += 0.2165063509461096*(f[3]*alpha_vdim[72]+f[8]*alpha_vdim[67]+f[0]*alpha_vdim[65]+f[1]*alpha_vdim[64]+alpha_cdim[4]*f[16]+alpha_cdim[0]*f[5]); 
  out[14] += 0.2165063509461096*(f[22]*alpha_vdim[72]+f[0]*alpha_cdim[69]+f[9]*alpha_vdim[67]+f[7]*alpha_vdim[65]+f[2]*alpha_vdim[64]+f[5]*alpha_cdim[64]); 
  out[15] += 0.2165063509461096*(f[21]*alpha_cdim[134]+f[5]*alpha_cdim[128]+f[1]*alpha_vdim[72]+f[0]*alpha_vdim[67]+f[8]*alpha_vdim[65]+f[3]*alpha_vdim[64]); 
  out[16] += 0.2165063509461096*(f[24]*alpha_vdim[72]+f[12]*alpha_vdim[67]+f[10]*alpha_vdim[65]+f[4]*alpha_vdim[64]+alpha_vdim[9]*f[28]+alpha_vdim[3]*f[15]+alpha_vdim[2]*f[14]+alpha_vdim[0]*f[5]); 
  out[17] += 0.2165063509461096*(f[2]*alpha_vdim[135]+f[7]*alpha_vdim[130]+f[0]*alpha_vdim[129]+f[1]*alpha_vdim[128]+alpha_cdim[4]*f[20]+alpha_cdim[0]*f[6]); 
  out[18] += 0.2165063509461096*(f[1]*alpha_vdim[135]+f[0]*alpha_vdim[130]+f[7]*alpha_vdim[129]+f[2]*alpha_vdim[128]+f[21]*alpha_cdim[69]+f[6]*alpha_cdim[64]); 
  out[19] += 0.2165063509461096*(f[22]*alpha_vdim[135]+f[0]*alpha_cdim[134]+f[9]*alpha_vdim[130]+f[8]*alpha_vdim[129]+f[3]*alpha_vdim[128]+f[6]*alpha_cdim[128]); 
  out[20] += 0.2165063509461096*(f[23]*alpha_vdim[135]+f[11]*alpha_vdim[130]+f[10]*alpha_vdim[129]+f[4]*alpha_vdim[128]+alpha_vdim[9]*f[34]+alpha_vdim[3]*f[19]+alpha_vdim[2]*f[18]+alpha_vdim[0]*f[6]); 
  out[21] += 0.2165063509461096*(f[26]*alpha_vdim[135]+f[14]*alpha_vdim[130]+f[13]*alpha_vdim[129]+f[5]*alpha_vdim[128]+f[33]*alpha_vdim[72]+f[19]*alpha_vdim[67]+f[17]*alpha_vdim[65]+f[6]*alpha_vdim[64]); 
  out[22] += 0.2165063509461096*(f[32]*alpha_cdim[134]+f[7]*alpha_cdim[128]+f[27]*alpha_cdim[69]+f[8]*alpha_cdim[64]+alpha_cdim[4]*f[25]+alpha_cdim[0]*f[9]); 
  out[23] += 0.2165063509461096*(f[29]*alpha_cdim[69]+f[10]*alpha_cdim[64]+alpha_vdim[3]*f[22]+alpha_cdim[0]*f[11]+f[8]*alpha_vdim[9]+alpha_vdim[0]*f[7]+f[2]*alpha_cdim[4]+f[1]*alpha_vdim[2]); 
  out[24] += 0.2165063509461096*(f[35]*alpha_cdim[134]+f[10]*alpha_cdim[128]+alpha_vdim[2]*f[22]+alpha_cdim[0]*f[12]+f[7]*alpha_vdim[9]+alpha_vdim[0]*f[8]+f[3]*alpha_cdim[4]+f[1]*alpha_vdim[3]); 
  out[25] += 0.2165063509461096*(f[36]*alpha_cdim[134]+f[11]*alpha_cdim[128]+f[31]*alpha_cdim[69]+f[12]*alpha_cdim[64]+alpha_vdim[0]*f[9]+f[0]*alpha_vdim[9]+alpha_vdim[2]*f[3]+f[2]*alpha_vdim[3]); 
  out[26] += 0.2165063509461096*(f[9]*alpha_vdim[72]+f[1]*alpha_cdim[69]+f[22]*alpha_vdim[67]+f[2]*alpha_vdim[65]+f[7]*alpha_vdim[64]+f[13]*alpha_cdim[64]+alpha_cdim[4]*f[30]+alpha_cdim[0]*f[14]); 
  out[27] += 0.2165063509461096*(f[38]*alpha_cdim[134]+f[13]*alpha_cdim[128]+f[0]*alpha_vdim[72]+f[1]*alpha_vdim[67]+f[3]*alpha_vdim[65]+f[8]*alpha_vdim[64]+alpha_cdim[4]*f[31]+alpha_cdim[0]*f[15]); 
  out[28] += 0.2165063509461096*(f[39]*alpha_cdim[134]+f[14]*alpha_cdim[128]+f[7]*alpha_vdim[72]+f[3]*alpha_cdim[69]+f[2]*alpha_vdim[67]+f[22]*alpha_vdim[65]+f[9]*alpha_vdim[64]+f[15]*alpha_cdim[64]); 
  out[29] += 0.2165063509461096*(f[12]*alpha_vdim[72]+f[24]*alpha_vdim[67]+f[4]*alpha_vdim[65]+f[10]*alpha_vdim[64]+alpha_vdim[9]*f[43]+alpha_vdim[3]*f[27]+alpha_vdim[2]*f[26]+alpha_cdim[0]*f[16]+alpha_vdim[0]*f[13]+alpha_cdim[4]*f[5]); 
  out[30] += 0.2165063509461096*(f[42]*alpha_vdim[72]+f[4]*alpha_cdim[69]+f[25]*alpha_vdim[67]+f[23]*alpha_vdim[65]+f[11]*alpha_vdim[64]+f[16]*alpha_cdim[64]+alpha_vdim[3]*f[28]+alpha_vdim[9]*f[15]+alpha_vdim[0]*f[14]+alpha_vdim[2]*f[5]); 
  out[31] += 0.2165063509461096*(f[41]*alpha_cdim[134]+f[16]*alpha_cdim[128]+f[10]*alpha_vdim[72]+f[4]*alpha_vdim[67]+f[24]*alpha_vdim[65]+f[12]*alpha_vdim[64]+alpha_vdim[2]*f[28]+alpha_vdim[0]*f[15]+alpha_vdim[9]*f[14]+alpha_vdim[3]*f[5]); 
  out[32] += 0.2165063509461096*(f[0]*alpha_vdim[135]+f[1]*alpha_vdim[130]+f[2]*alpha_vdim[129]+f[7]*alpha_vdim[128]+f[38]*alpha_cdim[69]+f[17]*alpha_cdim[64]+alpha_cdim[4]*f[36]+alpha_cdim[0]*f[18]); 
  out[33] += 0.2165063509461096*(f[9]*alpha_vdim[135]+f[1]*alpha_cdim[134]+f[22]*alpha_vdim[130]+f[3]*alpha_vdim[129]+f[8]*alpha_vdim[128]+f[17]*alpha_cdim[128]+alpha_cdim[4]*f[37]+alpha_cdim[0]*f[19]); 
  out[34] += 0.2165063509461096*(f[8]*alpha_vdim[135]+f[2]*alpha_cdim[134]+f[3]*alpha_vdim[130]+f[22]*alpha_vdim[129]+f[9]*alpha_vdim[128]+f[18]*alpha_cdim[128]+f[40]*alpha_cdim[69]+f[19]*alpha_cdim[64]); 
  out[35] += 0.2165063509461096*(f[11]*alpha_vdim[135]+f[23]*alpha_vdim[130]+f[4]*alpha_vdim[129]+f[10]*alpha_vdim[128]+alpha_vdim[9]*f[47]+alpha_vdim[3]*f[33]+alpha_vdim[2]*f[32]+alpha_cdim[0]*f[20]+alpha_vdim[0]*f[17]+alpha_cdim[4]*f[6]); 
  out[36] += 0.2165063509461096*(f[10]*alpha_vdim[135]+f[4]*alpha_vdim[130]+f[23]*alpha_vdim[129]+f[11]*alpha_vdim[128]+f[41]*alpha_cdim[69]+f[20]*alpha_cdim[64]+alpha_vdim[3]*f[34]+alpha_vdim[9]*f[19]+alpha_vdim[0]*f[18]+alpha_vdim[2]*f[6]); 
  out[37] += 0.2165063509461096*(f[42]*alpha_vdim[135]+f[4]*alpha_cdim[134]+f[25]*alpha_vdim[130]+f[24]*alpha_vdim[129]+f[12]*alpha_vdim[128]+f[20]*alpha_cdim[128]+alpha_vdim[2]*f[34]+alpha_vdim[0]*f[19]+alpha_vdim[9]*f[18]+alpha_vdim[3]*f[6]); 
  out[38] += 0.2165063509461096*(f[14]*alpha_vdim[135]+f[26]*alpha_vdim[130]+f[5]*alpha_vdim[129]+f[13]*alpha_vdim[128]+f[19]*alpha_vdim[72]+f[33]*alpha_vdim[67]+f[6]*alpha_vdim[65]+f[17]*alpha_vdim[64]+alpha_cdim[4]*f[41]+alpha_cdim[0]*f[21]); 
  out[39] += 0.2165063509461096*(f[13]*alpha_vdim[135]+f[5]*alpha_vdim[130]+f[26]*alpha_vdim[129]+f[14]*alpha_vdim[128]+f[47]*alpha_vdim[72]+f[6]*alpha_cdim[69]+f[34]*alpha_vdim[67]+f[32]*alpha_vdim[65]+f[18]*alpha_vdim[64]+f[21]*alpha_cdim[64]); 
  out[40] += 0.2165063509461096*(f[43]*alpha_vdim[135]+f[5]*alpha_cdim[134]+f[28]*alpha_vdim[130]+f[27]*alpha_vdim[129]+f[15]*alpha_vdim[128]+f[21]*alpha_cdim[128]+f[17]*alpha_vdim[72]+f[6]*alpha_vdim[67]+f[33]*alpha_vdim[65]+f[19]*alpha_vdim[64]); 
  out[41] += 0.2165063509461096*(f[44]*alpha_vdim[135]+f[30]*alpha_vdim[130]+f[29]*alpha_vdim[129]+f[16]*alpha_vdim[128]+f[49]*alpha_vdim[72]+f[37]*alpha_vdim[67]+f[35]*alpha_vdim[65]+f[20]*alpha_vdim[64]+alpha_vdim[9]*f[53]+alpha_vdim[3]*f[40]+alpha_vdim[2]*f[39]+alpha_vdim[0]*f[21]); 
  out[42] += 0.2165063509461096*(f[48]*alpha_cdim[134]+f[23]*alpha_cdim[128]+f[45]*alpha_cdim[69]+f[24]*alpha_cdim[64]+alpha_cdim[0]*f[25]+alpha_vdim[0]*f[22]+alpha_cdim[4]*f[9]+f[1]*alpha_vdim[9]+alpha_vdim[2]*f[8]+alpha_vdim[3]*f[7]); 
  out[43] += 0.2165063509461096*(f[51]*alpha_cdim[134]+f[26]*alpha_cdim[128]+f[2]*alpha_vdim[72]+f[8]*alpha_cdim[69]+f[7]*alpha_vdim[67]+f[9]*alpha_vdim[65]+f[22]*alpha_vdim[64]+f[27]*alpha_cdim[64]+alpha_cdim[4]*f[46]+alpha_cdim[0]*f[28]); 
  out[44] += 0.2165063509461096*(f[25]*alpha_vdim[72]+f[10]*alpha_cdim[69]+f[42]*alpha_vdim[67]+f[11]*alpha_vdim[65]+f[23]*alpha_vdim[64]+f[29]*alpha_cdim[64]+alpha_vdim[3]*f[43]+alpha_cdim[0]*f[30]+alpha_vdim[9]*f[27]+alpha_vdim[0]*f[26]+alpha_cdim[4]*f[14]+alpha_vdim[2]*f[13]); 
  out[45] += 0.2165063509461096*(f[54]*alpha_cdim[134]+f[29]*alpha_cdim[128]+f[4]*alpha_vdim[72]+f[10]*alpha_vdim[67]+f[12]*alpha_vdim[65]+f[24]*alpha_vdim[64]+alpha_vdim[2]*f[43]+alpha_cdim[0]*f[31]+alpha_vdim[0]*f[27]+alpha_vdim[9]*f[26]+alpha_cdim[4]*f[15]+alpha_vdim[3]*f[13]); 
  out[46] += 0.2165063509461096*(f[55]*alpha_cdim[134]+f[30]*alpha_cdim[128]+f[23]*alpha_vdim[72]+f[12]*alpha_cdim[69]+f[11]*alpha_vdim[67]+f[42]*alpha_vdim[65]+f[25]*alpha_vdim[64]+f[31]*alpha_cdim[64]+alpha_vdim[0]*f[28]+alpha_vdim[2]*f[15]+alpha_vdim[3]*f[14]+f[5]*alpha_vdim[9]); 
  out[47] += 0.2165063509461096*(f[3]*alpha_vdim[135]+f[7]*alpha_cdim[134]+f[8]*alpha_vdim[130]+f[9]*alpha_vdim[129]+f[22]*alpha_vdim[128]+f[32]*alpha_cdim[128]+f[52]*alpha_cdim[69]+f[33]*alpha_cdim[64]+alpha_cdim[4]*f[50]+alpha_cdim[0]*f[34]); 
  out[48] += 0.2165063509461096*(f[4]*alpha_vdim[135]+f[10]*alpha_vdim[130]+f[11]*alpha_vdim[129]+f[23]*alpha_vdim[128]+f[54]*alpha_cdim[69]+f[35]*alpha_cdim[64]+alpha_vdim[3]*f[47]+alpha_cdim[0]*f[36]+alpha_vdim[9]*f[33]+alpha_vdim[0]*f[32]+alpha_cdim[4]*f[18]+alpha_vdim[2]*f[17]); 
  out[49] += 0.2165063509461096*(f[25]*alpha_vdim[135]+f[10]*alpha_cdim[134]+f[42]*alpha_vdim[130]+f[12]*alpha_vdim[129]+f[24]*alpha_vdim[128]+f[35]*alpha_cdim[128]+alpha_vdim[2]*f[47]+alpha_cdim[0]*f[37]+alpha_vdim[0]*f[33]+alpha_vdim[9]*f[32]+alpha_cdim[4]*f[19]+alpha_vdim[3]*f[17]); 
  out[50] += 0.2165063509461096*(f[24]*alpha_vdim[135]+f[11]*alpha_cdim[134]+f[12]*alpha_vdim[130]+f[42]*alpha_vdim[129]+f[25]*alpha_vdim[128]+f[36]*alpha_cdim[128]+f[56]*alpha_cdim[69]+f[37]*alpha_cdim[64]+alpha_vdim[0]*f[34]+alpha_vdim[2]*f[19]+alpha_vdim[3]*f[18]+f[6]*alpha_vdim[9]); 
  out[51] += 0.2165063509461096*(f[5]*alpha_vdim[135]+f[13]*alpha_vdim[130]+f[14]*alpha_vdim[129]+f[26]*alpha_vdim[128]+f[34]*alpha_vdim[72]+f[17]*alpha_cdim[69]+f[47]*alpha_vdim[67]+f[18]*alpha_vdim[65]+f[32]*alpha_vdim[64]+f[38]*alpha_cdim[64]+alpha_cdim[4]*f[55]+alpha_cdim[0]*f[39]); 
  out[52] += 0.2165063509461096*(f[28]*alpha_vdim[135]+f[13]*alpha_cdim[134]+f[43]*alpha_vdim[130]+f[15]*alpha_vdim[129]+f[27]*alpha_vdim[128]+f[38]*alpha_cdim[128]+f[6]*alpha_vdim[72]+f[17]*alpha_vdim[67]+f[19]*alpha_vdim[65]+f[33]*alpha_vdim[64]+alpha_cdim[4]*f[56]+alpha_cdim[0]*f[40]); 
  out[53] += 0.2165063509461096*(f[27]*alpha_vdim[135]+f[14]*alpha_cdim[134]+f[15]*alpha_vdim[130]+f[43]*alpha_vdim[129]+f[28]*alpha_vdim[128]+f[39]*alpha_cdim[128]+f[32]*alpha_vdim[72]+f[19]*alpha_cdim[69]+f[18]*alpha_vdim[67]+f[47]*alpha_vdim[65]+f[34]*alpha_vdim[64]+f[40]*alpha_cdim[64]); 
  out[54] += 0.2165063509461096*(f[30]*alpha_vdim[135]+f[44]*alpha_vdim[130]+f[16]*alpha_vdim[129]+f[29]*alpha_vdim[128]+f[37]*alpha_vdim[72]+f[49]*alpha_vdim[67]+f[20]*alpha_vdim[65]+f[35]*alpha_vdim[64]+alpha_vdim[9]*f[59]+alpha_vdim[3]*f[52]+alpha_vdim[2]*f[51]+alpha_cdim[0]*f[41]+alpha_vdim[0]*f[38]+alpha_cdim[4]*f[21]); 
  out[55] += 0.2165063509461096*(f[29]*alpha_vdim[135]+f[16]*alpha_vdim[130]+f[44]*alpha_vdim[129]+f[30]*alpha_vdim[128]+f[58]*alpha_vdim[72]+f[20]*alpha_cdim[69]+f[50]*alpha_vdim[67]+f[48]*alpha_vdim[65]+f[36]*alpha_vdim[64]+f[41]*alpha_cdim[64]+alpha_vdim[3]*f[53]+alpha_vdim[9]*f[40]+alpha_vdim[0]*f[39]+alpha_vdim[2]*f[21]); 
  out[56] += 0.2165063509461096*(f[57]*alpha_vdim[135]+f[16]*alpha_cdim[134]+f[46]*alpha_vdim[130]+f[45]*alpha_vdim[129]+f[31]*alpha_vdim[128]+f[41]*alpha_cdim[128]+f[35]*alpha_vdim[72]+f[20]*alpha_vdim[67]+f[49]*alpha_vdim[65]+f[37]*alpha_vdim[64]+alpha_vdim[2]*f[53]+alpha_vdim[0]*f[40]+alpha_vdim[9]*f[39]+alpha_vdim[3]*f[21]); 
  out[57] += 0.2165063509461096*(f[60]*alpha_cdim[134]+f[44]*alpha_cdim[128]+f[11]*alpha_vdim[72]+f[24]*alpha_cdim[69]+f[23]*alpha_vdim[67]+f[25]*alpha_vdim[65]+f[42]*alpha_vdim[64]+f[45]*alpha_cdim[64]+alpha_cdim[0]*f[46]+alpha_vdim[0]*f[43]+alpha_cdim[4]*f[28]+alpha_vdim[2]*f[27]+alpha_vdim[3]*f[26]+alpha_vdim[9]*f[13]); 
  out[58] += 0.2165063509461096*(f[12]*alpha_vdim[135]+f[23]*alpha_cdim[134]+f[24]*alpha_vdim[130]+f[25]*alpha_vdim[129]+f[42]*alpha_vdim[128]+f[48]*alpha_cdim[128]+f[61]*alpha_cdim[69]+f[49]*alpha_cdim[64]+alpha_cdim[0]*f[50]+alpha_vdim[0]*f[47]+alpha_cdim[4]*f[34]+alpha_vdim[2]*f[33]+alpha_vdim[3]*f[32]+alpha_vdim[9]*f[17]); 
  out[59] += 0.2165063509461096*(f[15]*alpha_vdim[135]+f[26]*alpha_cdim[134]+f[27]*alpha_vdim[130]+f[28]*alpha_vdim[129]+f[43]*alpha_vdim[128]+f[51]*alpha_cdim[128]+f[18]*alpha_vdim[72]+f[33]*alpha_cdim[69]+f[32]*alpha_vdim[67]+f[34]*alpha_vdim[65]+f[47]*alpha_vdim[64]+f[52]*alpha_cdim[64]+alpha_cdim[4]*f[62]+alpha_cdim[0]*f[53]); 
  out[60] += 0.2165063509461096*(f[16]*alpha_vdim[135]+f[29]*alpha_vdim[130]+f[30]*alpha_vdim[129]+f[44]*alpha_vdim[128]+f[50]*alpha_vdim[72]+f[35]*alpha_cdim[69]+f[58]*alpha_vdim[67]+f[36]*alpha_vdim[65]+f[48]*alpha_vdim[64]+f[54]*alpha_cdim[64]+alpha_vdim[3]*f[59]+alpha_cdim[0]*f[55]+alpha_vdim[9]*f[52]+alpha_vdim[0]*f[51]+alpha_cdim[4]*f[39]+alpha_vdim[2]*f[38]); 
  out[61] += 0.2165063509461096*(f[46]*alpha_vdim[135]+f[29]*alpha_cdim[134]+f[57]*alpha_vdim[130]+f[31]*alpha_vdim[129]+f[45]*alpha_vdim[128]+f[54]*alpha_cdim[128]+f[20]*alpha_vdim[72]+f[35]*alpha_vdim[67]+f[37]*alpha_vdim[65]+f[49]*alpha_vdim[64]+alpha_vdim[2]*f[59]+alpha_cdim[0]*f[56]+alpha_vdim[0]*f[52]+alpha_vdim[9]*f[51]+alpha_cdim[4]*f[40]+alpha_vdim[3]*f[38]); 
  out[62] += 0.2165063509461096*(f[45]*alpha_vdim[135]+f[30]*alpha_cdim[134]+f[31]*alpha_vdim[130]+f[57]*alpha_vdim[129]+f[46]*alpha_vdim[128]+f[55]*alpha_cdim[128]+f[48]*alpha_vdim[72]+f[37]*alpha_cdim[69]+f[36]*alpha_vdim[67]+f[58]*alpha_vdim[65]+f[50]*alpha_vdim[64]+f[56]*alpha_cdim[64]+alpha_vdim[0]*f[53]+alpha_vdim[2]*f[40]+alpha_vdim[3]*f[39]+alpha_vdim[9]*f[21]); 
  out[63] += 0.2165063509461096*(f[31]*alpha_vdim[135]+f[44]*alpha_cdim[134]+f[45]*alpha_vdim[130]+f[46]*alpha_vdim[129]+f[57]*alpha_vdim[128]+f[60]*alpha_cdim[128]+f[36]*alpha_vdim[72]+f[49]*alpha_cdim[69]+f[48]*alpha_vdim[67]+f[50]*alpha_vdim[65]+f[58]*alpha_vdim[64]+f[61]*alpha_cdim[64]+alpha_cdim[0]*f[62]+alpha_vdim[0]*f[59]+alpha_cdim[4]*f[53]+alpha_vdim[2]*f[52]+alpha_vdim[3]*f[51]+alpha_vdim[9]*f[38]); 

  return alpha_mid; 
} 

