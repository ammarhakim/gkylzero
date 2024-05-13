#include <gkyl_vlasov_poisson_kernels.h> 

GKYL_CU_DH double vlasov_poisson_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out) 
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
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  const double dx11 = 2/dxv[1]; 

  const double *phi = &field[0]; 

  double cflFreq_mid = 0.0; 

  double alpha_cdim[96] = {0.0}; 
  double alpha_vdim[96] = {0.0}; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[3] = 2.309401076758503*dv0dx0; 

  cflFreq_mid += 5.0*(fabs(w0dx0)+0.5*dv0dx0); 

  alpha_cdim[48] = 8.0*w1dx1; 
  alpha_cdim[52] = 2.309401076758503*dv1dx1; 

  cflFreq_mid += 5.0*(fabs(w1dx1)+0.5*dv1dx1); 

  alpha_vdim[0] = -3.464101615137754*phi[1]*dv10*dx10; 
  alpha_vdim[1] = -7.745966692414834*phi[4]*dv10*dx10; 
  alpha_vdim[2] = -3.464101615137754*phi[3]*dv10*dx10; 
  alpha_vdim[5] = -7.745966692414834*phi[6]*dv10*dx10; 
  alpha_vdim[12] = -3.464101615137754*phi[7]*dv10*dx10; 

  cflFreq_mid += 5.0*fabs(0.125*alpha_vdim[0]-0.1397542485937369*alpha_vdim[12]); 

  alpha_vdim[48] = -3.464101615137754*phi[2]*dv11*dx11; 
  alpha_vdim[49] = -3.464101615137754*phi[3]*dv11*dx11; 
  alpha_vdim[50] = -7.745966692414834*phi[5]*dv11*dx11; 
  alpha_vdim[53] = -7.745966692414834*phi[7]*dv11*dx11; 
  alpha_vdim[59] = -3.464101615137754*phi[6]*dv11*dx11; 

  cflFreq_mid += 5.0*fabs(0.125*alpha_vdim[0]-0.1397542485937369*alpha_vdim[11]); 

  out[1] += 0.4330127018922193*(alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha_cdim[4]*f[4]+alpha_cdim[0]*f[0]); 
  out[3] += 0.4330127018922193*(alpha_vdim[12]*f[12]+alpha_vdim[5]*f[5]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.4330127018922193*(alpha_vdim[11]*f[11]+alpha_vdim[5]*f[5]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[5] += 0.4330127018922193*(alpha_cdim[4]*f[8]+alpha_cdim[3]*f[7]+alpha_cdim[0]*(f[2]+f[1])); 
  out[6] += 0.4330127018922193*alpha_vdim[12]*f[20]+0.3872983346207416*(alpha_vdim[5]*f[19]+alpha_cdim[3]*f[13]+alpha_vdim[1]*f[11])+0.4330127018922193*(alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[7] += 0.3872983346207416*(alpha_vdim[5]*f[20]+alpha_vdim[2]*f[12]+f[2]*alpha_vdim[12])+0.4330127018922193*(alpha_cdim[4]*f[10]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_cdim[0]*f[3]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[8] += 0.3872983346207416*(alpha_vdim[5]*f[19]+alpha_vdim[1]*f[11]+f[1]*alpha_vdim[11])+0.4330127018922193*(alpha_cdim[3]*f[10]+alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_cdim[0]*f[4]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[9] += 0.3872983346207416*alpha_vdim[5]*f[20]+0.4330127018922193*alpha_vdim[11]*f[19]+0.3872983346207416*(alpha_cdim[4]*f[14]+alpha_vdim[2]*f[12])+0.4330127018922193*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_cdim[0]*f[4]+f[0]*alpha_cdim[4]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[10] += 0.4330127018922193*(alpha_vdim[12]*f[26]+alpha_vdim[11]*f[21]+alpha_vdim[5]*(f[16]+f[15])+alpha_vdim[2]*f[9]+alpha_vdim[1]*f[8]+alpha_vdim[2]*f[7]+alpha_vdim[1]*f[6]+alpha_vdim[0]*(f[4]+f[3])); 
  out[11] += 0.9682458365518543*(alpha_cdim[3]*f[6]+alpha_cdim[0]*f[1]); 
  out[12] += 0.9682458365518543*(alpha_cdim[4]*f[9]+alpha_cdim[0]*f[2]); 
  out[13] += 0.9682458365518543*(alpha_vdim[12]*f[22]+alpha_vdim[5]*f[15]+alpha_vdim[2]*f[7]+alpha_vdim[1]*f[6]+alpha_vdim[0]*f[3]); 
  out[14] += 0.9682458365518543*(alpha_vdim[11]*f[25]+alpha_vdim[5]*f[16]+alpha_vdim[2]*f[9]+alpha_vdim[1]*f[8]+alpha_vdim[0]*f[4]); 
  out[15] += 0.3872983346207416*(alpha_cdim[3]*f[24]+alpha_vdim[2]*f[20]+alpha_vdim[1]*f[19])+0.4330127018922193*alpha_cdim[4]*f[17]+0.3872983346207416*(alpha_vdim[5]*f[12]+f[5]*alpha_vdim[12]+alpha_vdim[5]*f[11])+0.4330127018922193*(alpha_cdim[0]*(f[7]+f[6])+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+f[2]*(alpha_cdim[3]+alpha_vdim[1])+f[1]*alpha_vdim[2]); 
  out[16] += 0.3872983346207416*(alpha_cdim[4]*f[28]+alpha_vdim[2]*f[20]+alpha_vdim[1]*f[19])+0.4330127018922193*alpha_cdim[3]*f[18]+0.3872983346207416*(alpha_vdim[5]*(f[12]+f[11])+f[5]*alpha_vdim[11])+0.4330127018922193*(alpha_cdim[0]*(f[9]+f[8])+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+f[1]*alpha_cdim[4]+alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[17] += 0.4330127018922193*alpha_vdim[12]*f[36]+0.3872983346207416*(alpha_vdim[5]*(f[35]+f[32])+alpha_cdim[3]*f[27]+alpha_vdim[1]*(f[25]+f[21]))+0.4330127018922193*alpha_vdim[2]*(f[16]+f[15])+0.3872983346207416*f[6]*alpha_vdim[11]+0.4330127018922193*(alpha_cdim[0]*f[10]+alpha_vdim[5]*f[9]+alpha_vdim[0]*f[8]+alpha_vdim[5]*f[7]+alpha_vdim[0]*f[6]+alpha_vdim[1]*(f[4]+f[3])+alpha_cdim[3]*f[4]); 
  out[18] += 0.3872983346207416*alpha_vdim[5]*(f[36]+f[33])+0.4330127018922193*alpha_vdim[11]*f[32]+0.3872983346207416*(alpha_cdim[4]*f[30]+alpha_vdim[2]*(f[26]+f[22]))+0.4330127018922193*alpha_vdim[1]*(f[16]+f[15])+0.3872983346207416*f[9]*alpha_vdim[12]+0.4330127018922193*(alpha_cdim[0]*f[10]+alpha_vdim[0]*f[9]+alpha_vdim[5]*f[8]+alpha_vdim[0]*f[7]+alpha_vdim[5]*f[6]+alpha_vdim[2]*f[4]+f[3]*(alpha_cdim[4]+alpha_vdim[2])); 
  out[19] += 0.4330127018922193*alpha_cdim[4]*f[25]+0.9682458365518543*alpha_cdim[3]*f[15]+alpha_cdim[0]*(0.4330127018922193*f[11]+0.9682458365518543*f[5]); 
  out[20] += 0.4330127018922193*alpha_cdim[3]*f[22]+0.9682458365518543*alpha_cdim[4]*f[16]+alpha_cdim[0]*(0.4330127018922193*f[12]+0.9682458365518543*f[5]); 
  out[21] += 0.8660254037844386*alpha_cdim[3]*f[23]+0.4330127018922193*(alpha_vdim[2]*f[19]+alpha_vdim[0]*f[11])+0.9682458365518543*alpha_cdim[0]*f[6]+0.3872983346207416*alpha_vdim[5]*f[5]+f[1]*(0.9682458365518543*alpha_cdim[3]+0.3872983346207416*alpha_vdim[1]); 
  out[22] += 0.4330127018922193*alpha_vdim[1]*f[20]+0.9682458365518543*alpha_cdim[4]*f[18]+0.276641667586244*alpha_vdim[12]*f[12]+0.4330127018922193*(alpha_vdim[0]*f[12]+f[0]*alpha_vdim[12])+0.9682458365518543*alpha_cdim[0]*f[7]+0.3872983346207416*(alpha_vdim[5]*f[5]+alpha_vdim[2]*f[2]); 
  out[23] += 0.9682458365518543*alpha_vdim[12]*f[33]+0.8660254037844386*(alpha_vdim[5]*f[32]+alpha_vdim[1]*f[21])+0.9682458365518543*alpha_vdim[2]*f[15]+0.4330127018922193*alpha_cdim[0]*f[13]+0.9682458365518543*(alpha_vdim[5]*f[7]+alpha_vdim[0]*f[6])+(0.3872983346207416*alpha_cdim[3]+0.9682458365518543*alpha_vdim[1])*f[3]; 
  out[24] += 0.8660254037844386*alpha_vdim[5]*f[33]+0.4330127018922193*alpha_cdim[4]*f[27]+0.8660254037844386*alpha_vdim[2]*f[22]+0.9682458365518543*alpha_vdim[1]*f[15]+0.4330127018922193*alpha_cdim[0]*f[13]+0.8660254037844386*f[7]*alpha_vdim[12]+0.9682458365518543*(alpha_vdim[0]*f[7]+alpha_vdim[5]*f[6]+alpha_vdim[2]*f[3]); 
  out[25] += 0.4330127018922193*alpha_vdim[2]*f[19]+0.9682458365518543*alpha_cdim[3]*f[17]+0.276641667586244*alpha_vdim[11]*f[11]+0.4330127018922193*(alpha_vdim[0]*f[11]+f[0]*alpha_vdim[11])+0.9682458365518543*alpha_cdim[0]*f[8]+0.3872983346207416*(alpha_vdim[5]*f[5]+alpha_vdim[1]*f[1]); 
  out[26] += 0.8660254037844386*alpha_cdim[4]*f[29]+0.4330127018922193*(alpha_vdim[1]*f[20]+alpha_vdim[0]*f[12])+0.9682458365518543*alpha_cdim[0]*f[9]+0.3872983346207416*alpha_vdim[5]*f[5]+f[2]*(0.9682458365518543*alpha_cdim[4]+0.3872983346207416*alpha_vdim[2]); 
  out[27] += 0.9682458365518543*alpha_vdim[12]*f[38]+alpha_vdim[5]*(0.4330127018922193*f[34]+0.9682458365518543*f[31])+0.4330127018922193*(alpha_vdim[2]*f[24]+alpha_vdim[1]*f[23])+0.9682458365518543*(alpha_vdim[2]*f[18]+alpha_vdim[1]*f[17])+alpha_vdim[0]*(0.4330127018922193*f[13]+0.9682458365518543*f[10]); 
  out[28] += 0.8660254037844386*alpha_vdim[5]*f[35]+0.4330127018922193*alpha_cdim[3]*f[30]+0.8660254037844386*alpha_vdim[1]*f[25]+0.9682458365518543*alpha_vdim[2]*f[16]+0.4330127018922193*alpha_cdim[0]*f[14]+0.8660254037844386*f[8]*alpha_vdim[11]+0.9682458365518543*(alpha_vdim[5]*f[9]+alpha_vdim[0]*f[8]+alpha_vdim[1]*f[4]); 
  out[29] += 0.8660254037844386*alpha_vdim[5]*f[36]+0.9682458365518543*alpha_vdim[11]*f[35]+0.8660254037844386*alpha_vdim[2]*f[26]+0.9682458365518543*alpha_vdim[1]*f[16]+0.4330127018922193*alpha_cdim[0]*f[14]+0.9682458365518543*(alpha_vdim[0]*f[9]+alpha_vdim[5]*f[8])+(0.3872983346207416*alpha_cdim[4]+0.9682458365518543*alpha_vdim[2])*f[4]; 
  out[30] += 0.4330127018922193*alpha_vdim[5]*f[41]+0.9682458365518543*(alpha_vdim[11]*f[37]+alpha_vdim[5]*f[31])+0.4330127018922193*(alpha_vdim[2]*f[29]+alpha_vdim[1]*f[28])+0.9682458365518543*(alpha_vdim[2]*f[18]+alpha_vdim[1]*f[17])+alpha_vdim[0]*(0.4330127018922193*f[14]+0.9682458365518543*f[10]); 
  out[31] += 0.3872983346207416*(alpha_cdim[4]*f[42]+alpha_cdim[3]*f[40]+alpha_vdim[2]*f[36]+alpha_vdim[1]*f[35]+alpha_vdim[2]*f[33]+alpha_vdim[1]*f[32]+alpha_vdim[5]*(f[26]+f[25]+f[22]+f[21]))+0.4330127018922193*alpha_cdim[0]*(f[18]+f[17])+(0.3872983346207416*alpha_vdim[12]+0.4330127018922193*alpha_vdim[0])*f[16]+0.3872983346207416*alpha_vdim[11]*f[15]+0.4330127018922193*(alpha_vdim[0]*f[15]+(alpha_cdim[3]+alpha_vdim[1])*f[9]+alpha_vdim[2]*f[8]+alpha_vdim[1]*f[7]+(alpha_cdim[4]+alpha_vdim[2])*f[6]+(f[4]+f[3])*alpha_vdim[5]); 
  out[32] += 0.4330127018922193*alpha_cdim[4]*f[37]+0.8660254037844386*alpha_cdim[3]*f[34]+0.4330127018922193*alpha_cdim[0]*f[21]+0.3464101615137755*alpha_vdim[5]*f[20]+(0.3872983346207416*alpha_vdim[12]+0.4330127018922193*alpha_vdim[0])*f[19]+0.9682458365518543*alpha_cdim[0]*f[15]+0.4330127018922193*alpha_vdim[2]*f[11]+0.9682458365518543*alpha_cdim[3]*f[5]+0.3872983346207416*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]); 
  out[33] += 0.9682458365518543*alpha_cdim[4]*f[31]+0.4330127018922193*alpha_cdim[0]*f[22]+(0.276641667586244*alpha_vdim[12]+0.4330127018922193*alpha_vdim[0])*f[20]+0.3464101615137755*alpha_vdim[5]*f[19]+0.9682458365518543*alpha_cdim[0]*f[15]+0.4330127018922193*((alpha_cdim[3]+alpha_vdim[1])*f[12]+f[1]*alpha_vdim[12])+0.3872983346207416*(alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]); 
  out[34] += 0.4330127018922193*alpha_cdim[4]*f[39]+0.8660254037844386*(alpha_vdim[2]*f[33]+alpha_vdim[1]*f[32])+0.4330127018922193*alpha_cdim[0]*(f[24]+f[23])+0.8660254037844386*alpha_vdim[5]*(f[22]+f[21])+(0.8660254037844386*alpha_vdim[12]+0.9682458365518543*alpha_vdim[0])*f[15]+0.3872983346207416*alpha_cdim[3]*f[7]+0.9682458365518543*(alpha_vdim[1]*f[7]+alpha_vdim[2]*f[6]+f[3]*alpha_vdim[5]); 
  out[35] += 0.9682458365518543*alpha_cdim[3]*f[31]+0.4330127018922193*alpha_cdim[0]*f[25]+0.3464101615137755*alpha_vdim[5]*f[20]+(0.276641667586244*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[19]+0.9682458365518543*alpha_cdim[0]*f[16]+0.4330127018922193*((alpha_cdim[4]+alpha_vdim[2])*f[11]+f[2]*alpha_vdim[11])+0.3872983346207416*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]); 
  out[36] += 0.8660254037844386*alpha_cdim[4]*f[41]+0.4330127018922193*(alpha_cdim[3]*f[38]+alpha_cdim[0]*f[26])+(0.3872983346207416*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[20]+0.3464101615137755*alpha_vdim[5]*f[19]+0.9682458365518543*alpha_cdim[0]*f[16]+0.4330127018922193*alpha_vdim[1]*f[12]+0.9682458365518543*alpha_cdim[4]*f[5]+0.3872983346207416*(alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]); 
  out[37] += 0.8660254037844386*alpha_cdim[3]*f[39]+0.4330127018922193*(alpha_vdim[2]*(f[35]+f[32])+alpha_vdim[0]*f[25])+(0.276641667586244*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[21]+0.9682458365518543*alpha_cdim[0]*f[17]+0.3872983346207416*alpha_vdim[5]*(f[16]+f[15])+0.4330127018922193*f[3]*alpha_vdim[11]+0.9682458365518543*alpha_cdim[3]*f[8]+0.3872983346207416*alpha_vdim[1]*(f[8]+f[6]); 
  out[38] += 0.8660254037844386*alpha_cdim[4]*f[43]+0.4330127018922193*alpha_vdim[1]*(f[36]+f[33])+0.276641667586244*alpha_vdim[12]*f[26]+0.4330127018922193*alpha_vdim[0]*(f[26]+f[22])+0.9682458365518543*alpha_cdim[0]*f[18]+0.3872983346207416*alpha_vdim[5]*(f[16]+f[15])+0.4330127018922193*f[4]*alpha_vdim[12]+0.3872983346207416*alpha_vdim[2]*f[9]+(0.9682458365518543*alpha_cdim[4]+0.3872983346207416*alpha_vdim[2])*f[7]; 
  out[39] += 0.9682458365518543*alpha_vdim[12]*f[45]+0.8660254037844386*(alpha_vdim[5]*f[44]+alpha_vdim[1]*f[37])+alpha_vdim[2]*(0.4330127018922193*f[34]+0.9682458365518543*f[31])+0.4330127018922193*(alpha_cdim[0]*f[27]+alpha_vdim[5]*f[24])+(0.3872983346207416*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[23]+0.9682458365518543*(alpha_vdim[5]*f[18]+alpha_vdim[0]*f[17])+0.4330127018922193*alpha_vdim[1]*f[13]+(0.3872983346207416*alpha_cdim[3]+0.9682458365518543*alpha_vdim[1])*f[10]; 
  out[40] += 0.8660254037844386*(alpha_vdim[5]*f[45]+alpha_vdim[2]*f[38])+alpha_vdim[1]*(0.4330127018922193*f[34]+0.9682458365518543*f[31])+0.4330127018922193*(alpha_cdim[0]*f[27]+alpha_vdim[0]*f[24]+alpha_vdim[5]*f[23])+0.8660254037844386*alpha_vdim[12]*f[18]+0.9682458365518543*(alpha_vdim[0]*f[18]+alpha_vdim[5]*f[17])+0.4330127018922193*alpha_cdim[4]*f[13]+alpha_vdim[2]*(0.4330127018922193*f[13]+0.9682458365518543*f[10]); 
  out[41] += 0.4330127018922193*alpha_cdim[3]*f[43]+0.8660254037844386*(alpha_vdim[2]*f[36]+alpha_vdim[1]*f[35])+0.4330127018922193*alpha_cdim[0]*(f[29]+f[28])+0.8660254037844386*(alpha_vdim[5]*(f[26]+f[25])+alpha_vdim[11]*f[16])+0.9682458365518543*(alpha_vdim[0]*f[16]+alpha_vdim[1]*f[9])+0.3872983346207416*alpha_cdim[4]*f[8]+0.9682458365518543*(alpha_vdim[2]*f[8]+f[4]*alpha_vdim[5]); 
  out[42] += 0.8660254037844386*alpha_vdim[5]*f[44]+0.4330127018922193*alpha_vdim[2]*f[41]+0.8660254037844386*alpha_vdim[1]*f[37]+0.9682458365518543*alpha_vdim[2]*f[31]+0.4330127018922193*(alpha_cdim[0]*f[30]+alpha_vdim[5]*f[29]+alpha_vdim[0]*f[28])+0.9682458365518543*alpha_vdim[5]*f[18]+(0.8660254037844386*alpha_vdim[11]+0.9682458365518543*alpha_vdim[0])*f[17]+0.4330127018922193*alpha_cdim[3]*f[14]+alpha_vdim[1]*(0.4330127018922193*f[14]+0.9682458365518543*f[10]); 
  out[43] += 0.8660254037844386*alpha_vdim[5]*f[45]+0.9682458365518543*alpha_vdim[11]*f[44]+0.4330127018922193*alpha_vdim[1]*f[41]+0.8660254037844386*alpha_vdim[2]*f[38]+0.9682458365518543*alpha_vdim[1]*f[31]+0.4330127018922193*alpha_cdim[0]*f[30]+0.3872983346207416*alpha_vdim[12]*f[29]+0.4330127018922193*(alpha_vdim[0]*f[29]+alpha_vdim[5]*f[28])+0.9682458365518543*(alpha_vdim[0]*f[18]+alpha_vdim[5]*f[17])+0.4330127018922193*alpha_vdim[2]*f[14]+(0.3872983346207416*alpha_cdim[4]+0.9682458365518543*alpha_vdim[2])*f[10]; 
  out[44] += 0.8660254037844386*alpha_cdim[3]*f[46]+0.4330127018922193*alpha_cdim[0]*f[37]+0.3464101615137755*alpha_vdim[5]*f[36]+(0.3872983346207416*alpha_vdim[12]+0.4330127018922193*alpha_vdim[0])*f[35]+0.3464101615137755*alpha_vdim[5]*f[33]+(0.276641667586244*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[32]+0.9682458365518543*alpha_cdim[0]*f[31]+0.4330127018922193*(alpha_vdim[2]*f[25]+(alpha_cdim[4]+alpha_vdim[2])*f[21])+0.9682458365518543*alpha_cdim[3]*f[16]+0.3872983346207416*alpha_vdim[1]*(f[16]+f[15])+0.4330127018922193*f[7]*alpha_vdim[11]+0.3872983346207416*alpha_vdim[5]*(f[8]+f[6]); 
  out[45] += 0.8660254037844386*alpha_cdim[4]*f[47]+0.4330127018922193*alpha_cdim[0]*f[38]+(0.276641667586244*alpha_vdim[12]+0.4330127018922193*alpha_vdim[0])*f[36]+0.3464101615137755*alpha_vdim[5]*f[35]+(0.3872983346207416*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[33]+0.3464101615137755*alpha_vdim[5]*f[32]+0.9682458365518543*alpha_cdim[0]*f[31]+0.4330127018922193*(alpha_vdim[1]*(f[26]+f[22])+alpha_cdim[3]*f[26])+0.3872983346207416*alpha_vdim[2]*f[16]+(0.9682458365518543*alpha_cdim[4]+0.3872983346207416*alpha_vdim[2])*f[15]+0.4330127018922193*f[8]*alpha_vdim[12]+0.3872983346207416*alpha_vdim[5]*(f[9]+f[7]); 
  out[46] += 0.8660254037844386*(alpha_vdim[2]*f[45]+alpha_vdim[1]*f[44])+0.4330127018922193*alpha_cdim[0]*(f[40]+f[39])+0.8660254037844386*alpha_vdim[5]*(f[38]+f[37])+(0.3872983346207416*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[34]+(0.8660254037844386*alpha_vdim[12]+0.9682458365518543*alpha_vdim[0])*f[31]+0.4330127018922193*(alpha_vdim[1]*f[24]+(alpha_cdim[4]+alpha_vdim[2])*f[23])+0.3872983346207416*alpha_cdim[3]*f[18]+0.9682458365518543*(alpha_vdim[1]*f[18]+alpha_vdim[2]*f[17])+alpha_vdim[5]*(0.4330127018922193*f[13]+0.9682458365518543*f[10]); 
  out[47] += 0.8660254037844386*(alpha_vdim[2]*f[45]+alpha_vdim[1]*f[44])+0.4330127018922193*alpha_cdim[0]*(f[43]+f[42])+(0.3872983346207416*alpha_vdim[12]+0.4330127018922193*alpha_vdim[0])*f[41]+0.8660254037844386*alpha_vdim[5]*(f[38]+f[37])+(0.8660254037844386*alpha_vdim[11]+0.9682458365518543*alpha_vdim[0])*f[31]+0.4330127018922193*((alpha_cdim[3]+alpha_vdim[1])*f[29]+alpha_vdim[2]*f[28])+0.9682458365518543*alpha_vdim[1]*f[18]+(0.3872983346207416*alpha_cdim[4]+0.9682458365518543*alpha_vdim[2])*f[17]+alpha_vdim[5]*(0.4330127018922193*f[14]+0.9682458365518543*f[10]); 

  return cflFreq_mid; 
} 

