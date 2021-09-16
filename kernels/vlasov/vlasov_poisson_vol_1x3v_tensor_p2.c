#include <gkyl_vlasov_kernels.h> 

GKYL_CU_DH double vlasov_poisson_vol_1x3v_tensor_p2(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fac_phi:   potential (scaled by appropriate factors).
  // vecA:      vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 
  const double *phi = &fac_phi[0]; 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv12 = 2/dxv[3]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[81]; 
  double alpha_vdim[243]; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[2] = 2.309401076758503*dv0dx0; 
  alpha_mid += fabs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = -4.898979485566357*phi[1]*dv10*dx10; 
  alpha_vdim[1] = -10.95445115010332*phi[2]*dv10*dx10; 
  alpha_mid += fabs(0.125*alpha_vdim[0]); 

  alpha_mid += fabs(0.0); 

  alpha_mid += fabs(0.0); 

  out[1] += 0.4330127018922193*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[5] += 0.3872983346207416*(alpha_cdim[2]*f[12]+alpha_vdim[1]*f[11])+0.4330127018922193*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.4330127018922193*(alpha_cdim[2]*f[7]+alpha_cdim[0]*f[3]); 
  out[7] += 0.4330127018922193*(alpha_vdim[1]*f[6]+alpha_vdim[0]*f[3]); 
  out[8] += 0.4330127018922193*(alpha_cdim[2]*f[9]+alpha_cdim[0]*f[4]); 
  out[9] += 0.4330127018922193*(alpha_vdim[1]*f[8]+alpha_vdim[0]*f[4]); 
  out[11] += 0.9682458365518543*(alpha_cdim[2]*f[5]+alpha_cdim[0]*f[1]); 
  out[12] += 0.9682458365518543*(alpha_vdim[1]*f[5]+alpha_vdim[0]*f[2]); 
  out[15] += 0.3872983346207416*(alpha_cdim[2]*f[22]+alpha_vdim[1]*f[21])+0.4330127018922193*(alpha_cdim[0]*f[7]+alpha_vdim[0]*f[6]+(alpha_cdim[2]+alpha_vdim[1])*f[3]); 
  out[16] += 0.3872983346207416*(alpha_cdim[2]*f[26]+alpha_vdim[1]*f[25])+0.4330127018922193*(alpha_cdim[0]*f[9]+alpha_vdim[0]*f[8]+(alpha_cdim[2]+alpha_vdim[1])*f[4]); 
  out[17] += 0.4330127018922193*(alpha_cdim[2]*f[18]+alpha_cdim[0]*f[10]); 
  out[18] += 0.4330127018922193*(alpha_vdim[1]*f[17]+alpha_vdim[0]*f[10]); 
  out[19] += 0.8660254037844386*alpha_cdim[2]*f[20]+0.4330127018922193*alpha_vdim[0]*f[11]+0.9682458365518543*alpha_cdim[0]*f[5]+f[1]*(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1]); 
  out[20] += 0.8660254037844386*alpha_vdim[1]*f[19]+0.4330127018922193*alpha_cdim[0]*f[12]+0.9682458365518543*alpha_vdim[0]*f[5]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[2]; 
  out[21] += 0.9682458365518543*(alpha_cdim[2]*f[15]+alpha_cdim[0]*f[6]); 
  out[22] += 0.9682458365518543*(alpha_vdim[1]*f[15]+alpha_vdim[0]*f[7]); 
  out[23] += 0.4330127018922193*(alpha_cdim[2]*f[24]+alpha_cdim[0]*f[13]); 
  out[24] += 0.4330127018922193*(alpha_vdim[1]*f[23]+alpha_vdim[0]*f[13]); 
  out[25] += 0.9682458365518543*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[8]); 
  out[26] += 0.9682458365518543*(alpha_vdim[1]*f[16]+alpha_vdim[0]*f[9]); 
  out[28] += 0.4330127018922193*(alpha_cdim[2]*f[29]+alpha_cdim[0]*f[14]); 
  out[29] += 0.4330127018922193*(alpha_vdim[1]*f[28]+alpha_vdim[0]*f[14]); 
  out[31] += 0.3872983346207416*(alpha_cdim[2]*f[38]+alpha_vdim[1]*f[37])+0.4330127018922193*(alpha_cdim[0]*f[18]+alpha_vdim[0]*f[17]+(alpha_cdim[2]+alpha_vdim[1])*f[10]); 
  out[32] += 0.8660254037844386*alpha_cdim[2]*f[33]+0.4330127018922193*alpha_vdim[0]*f[21]+0.9682458365518543*alpha_cdim[0]*f[15]+(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])*f[6]; 
  out[33] += 0.8660254037844386*alpha_vdim[1]*f[32]+0.4330127018922193*alpha_cdim[0]*f[22]+0.9682458365518543*alpha_vdim[0]*f[15]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[7]; 
  out[34] += 0.3872983346207416*(alpha_cdim[2]*f[46]+alpha_vdim[1]*f[45])+0.4330127018922193*(alpha_cdim[0]*f[24]+alpha_vdim[0]*f[23]+(alpha_cdim[2]+alpha_vdim[1])*f[13]); 
  out[35] += 0.8660254037844386*alpha_cdim[2]*f[36]+0.4330127018922193*alpha_vdim[0]*f[25]+0.9682458365518543*alpha_cdim[0]*f[16]+(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])*f[8]; 
  out[36] += 0.8660254037844386*alpha_vdim[1]*f[35]+0.4330127018922193*alpha_cdim[0]*f[26]+0.9682458365518543*alpha_vdim[0]*f[16]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[9]; 
  out[37] += 0.9682458365518543*(alpha_cdim[2]*f[31]+alpha_cdim[0]*f[17]); 
  out[38] += 0.9682458365518543*(alpha_vdim[1]*f[31]+alpha_vdim[0]*f[18]); 
  out[39] += 0.4330127018922193*(alpha_cdim[2]*f[40]+alpha_cdim[0]*f[27]); 
  out[40] += 0.4330127018922193*(alpha_vdim[1]*f[39]+alpha_vdim[0]*f[27]); 
  out[41] += 0.3872983346207416*(alpha_cdim[2]*f[48]+alpha_vdim[1]*f[47])+0.4330127018922193*(alpha_cdim[0]*f[29]+alpha_vdim[0]*f[28]+(alpha_cdim[2]+alpha_vdim[1])*f[14]); 
  out[42] += 0.4330127018922193*(alpha_cdim[2]*f[43]+alpha_cdim[0]*f[30]); 
  out[43] += 0.4330127018922193*(alpha_vdim[1]*f[42]+alpha_vdim[0]*f[30]); 
  out[44] += 0.9682458365518543*(alpha_cdim[0]*f[20]+alpha_vdim[0]*f[19])+0.8660254037844386*(alpha_cdim[2]+alpha_vdim[1])*f[5]; 
  out[45] += 0.9682458365518543*(alpha_cdim[2]*f[34]+alpha_cdim[0]*f[23]); 
  out[46] += 0.9682458365518543*(alpha_vdim[1]*f[34]+alpha_vdim[0]*f[24]); 
  out[47] += 0.9682458365518543*(alpha_cdim[2]*f[41]+alpha_cdim[0]*f[28]); 
  out[48] += 0.9682458365518543*(alpha_vdim[1]*f[41]+alpha_vdim[0]*f[29]); 
  out[50] += 0.8660254037844386*alpha_cdim[2]*f[51]+0.4330127018922193*alpha_vdim[0]*f[37]+0.9682458365518543*alpha_cdim[0]*f[31]+(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])*f[17]; 
  out[51] += 0.8660254037844386*alpha_vdim[1]*f[50]+0.4330127018922193*alpha_cdim[0]*f[38]+0.9682458365518543*alpha_vdim[0]*f[31]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[18]; 
  out[52] += 0.3872983346207416*(alpha_cdim[2]*f[59]+alpha_vdim[1]*f[58])+0.4330127018922193*(alpha_cdim[0]*f[40]+alpha_vdim[0]*f[39]+(alpha_cdim[2]+alpha_vdim[1])*f[27]); 
  out[53] += 0.3872983346207416*(alpha_cdim[2]*f[63]+alpha_vdim[1]*f[62])+0.4330127018922193*(alpha_cdim[0]*f[43]+alpha_vdim[0]*f[42]+(alpha_cdim[2]+alpha_vdim[1])*f[30]); 
  out[54] += 0.9682458365518543*(alpha_cdim[0]*f[33]+alpha_vdim[0]*f[32])+0.8660254037844386*(alpha_cdim[2]+alpha_vdim[1])*f[15]; 
  out[55] += 0.8660254037844386*alpha_cdim[2]*f[56]+0.4330127018922193*alpha_vdim[0]*f[45]+0.9682458365518543*alpha_cdim[0]*f[34]+(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])*f[23]; 
  out[56] += 0.8660254037844386*alpha_vdim[1]*f[55]+0.4330127018922193*alpha_cdim[0]*f[46]+0.9682458365518543*alpha_vdim[0]*f[34]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[24]; 
  out[57] += 0.9682458365518543*(alpha_cdim[0]*f[36]+alpha_vdim[0]*f[35])+0.8660254037844386*(alpha_cdim[2]+alpha_vdim[1])*f[16]; 
  out[58] += 0.9682458365518543*(alpha_cdim[2]*f[52]+alpha_cdim[0]*f[39]); 
  out[59] += 0.9682458365518543*(alpha_vdim[1]*f[52]+alpha_vdim[0]*f[40]); 
  out[60] += 0.8660254037844386*alpha_cdim[2]*f[61]+0.4330127018922193*alpha_vdim[0]*f[47]+0.9682458365518543*alpha_cdim[0]*f[41]+(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])*f[28]; 
  out[61] += 0.8660254037844386*alpha_vdim[1]*f[60]+0.4330127018922193*alpha_cdim[0]*f[48]+0.9682458365518543*alpha_vdim[0]*f[41]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[29]; 
  out[62] += 0.9682458365518543*(alpha_cdim[2]*f[53]+alpha_cdim[0]*f[42]); 
  out[63] += 0.9682458365518543*(alpha_vdim[1]*f[53]+alpha_vdim[0]*f[43]); 
  out[64] += 0.4330127018922193*(alpha_cdim[2]*f[65]+alpha_cdim[0]*f[49]); 
  out[65] += 0.4330127018922193*(alpha_vdim[1]*f[64]+alpha_vdim[0]*f[49]); 
  out[66] += 0.9682458365518543*(alpha_cdim[0]*f[51]+alpha_vdim[0]*f[50])+0.8660254037844386*(alpha_cdim[2]+alpha_vdim[1])*f[31]; 
  out[67] += 0.8660254037844386*alpha_cdim[2]*f[68]+0.4330127018922193*alpha_vdim[0]*f[58]+0.9682458365518543*alpha_cdim[0]*f[52]+(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])*f[39]; 
  out[68] += 0.8660254037844386*alpha_vdim[1]*f[67]+0.4330127018922193*alpha_cdim[0]*f[59]+0.9682458365518543*alpha_vdim[0]*f[52]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[40]; 
  out[69] += 0.8660254037844386*alpha_cdim[2]*f[70]+0.4330127018922193*alpha_vdim[0]*f[62]+0.9682458365518543*alpha_cdim[0]*f[53]+(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])*f[42]; 
  out[70] += 0.8660254037844386*alpha_vdim[1]*f[69]+0.4330127018922193*alpha_cdim[0]*f[63]+0.9682458365518543*alpha_vdim[0]*f[53]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[43]; 
  out[71] += 0.3872983346207416*(alpha_cdim[2]*f[75]+alpha_vdim[1]*f[74])+0.4330127018922193*(alpha_cdim[0]*f[65]+alpha_vdim[0]*f[64]+(alpha_cdim[2]+alpha_vdim[1])*f[49]); 
  out[72] += 0.9682458365518543*(alpha_cdim[0]*f[56]+alpha_vdim[0]*f[55])+0.8660254037844386*(alpha_cdim[2]+alpha_vdim[1])*f[34]; 
  out[73] += 0.9682458365518543*(alpha_cdim[0]*f[61]+alpha_vdim[0]*f[60])+0.8660254037844386*(alpha_cdim[2]+alpha_vdim[1])*f[41]; 
  out[74] += 0.9682458365518543*(alpha_cdim[2]*f[71]+alpha_cdim[0]*f[64]); 
  out[75] += 0.9682458365518543*(alpha_vdim[1]*f[71]+alpha_vdim[0]*f[65]); 
  out[76] += 0.9682458365518543*(alpha_cdim[0]*f[68]+alpha_vdim[0]*f[67])+0.8660254037844386*(alpha_cdim[2]+alpha_vdim[1])*f[52]; 
  out[77] += 0.9682458365518543*(alpha_cdim[0]*f[70]+alpha_vdim[0]*f[69])+0.8660254037844386*(alpha_cdim[2]+alpha_vdim[1])*f[53]; 
  out[78] += 0.8660254037844386*alpha_cdim[2]*f[79]+0.4330127018922193*alpha_vdim[0]*f[74]+0.9682458365518543*alpha_cdim[0]*f[71]+(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])*f[64]; 
  out[79] += 0.8660254037844386*alpha_vdim[1]*f[78]+0.4330127018922193*alpha_cdim[0]*f[75]+0.9682458365518543*alpha_vdim[0]*f[71]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[65]; 
  out[80] += 0.9682458365518543*(alpha_cdim[0]*f[79]+alpha_vdim[0]*f[78])+0.8660254037844386*(alpha_cdim[2]+alpha_vdim[1])*f[71]; 

  return alpha_mid; 
} 

