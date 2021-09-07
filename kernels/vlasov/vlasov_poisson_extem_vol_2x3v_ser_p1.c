#include <gkyl_vlasov_kernels.h> 

GKYL_CU_DH double vlasov_poisson_extem_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fac_phi:   potential (scaled by appropriate factors).
  // vecA:      vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  const double dx11 = 2/dxv[1]; 
  const double *phi = &fac_phi[0]; 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv11 = 2/dxv[3]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv12 = 2/dxv[4]; 
  const double dv3 = dxv[4], wv3 = w[4]; 

  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[4]; 
  const double *A2 = &vecA[8]; 
  double alpha_mid = 0.0; 
  double alpha_cdim[64]; 
  double alpha_vdim[96]; 

  alpha_cdim[0] = 11.31370849898477*w0dx0; 
  alpha_cdim[3] = 3.265986323710906*dv0dx0; 
  alpha_mid += fabs(w0dx0)+0.5*dv0dx0; 

  alpha_cdim[32] = 11.31370849898477*w1dx1; 
  alpha_cdim[36] = 3.265986323710906*dv1dx1; 
  alpha_mid += fabs(w1dx1)+0.5*dv1dx1; 

  alpha_vdim[0] = dv10*(4.898979485566357*A2[1]*dx10*wv3-4.898979485566357*A0[2]*dx11*wv2+dx10*(4.898979485566357*A1[1]*wv2-4.898979485566357*phi[1])); 
  alpha_vdim[1] = -4.898979485566357*A0[3]*dv10*dx11*wv2; 
  alpha_vdim[2] = dv10*dx10*(4.898979485566357*(A2[3]*wv3+A1[3]*wv2)-4.898979485566357*phi[3]); 
  alpha_vdim[4] = dv10*dv2*(1.414213562373095*A1[1]*dx10-1.414213562373095*A0[2]*dx11); 
  alpha_vdim[5] = 1.414213562373095*A2[1]*dv10*dv3*dx10; 
  alpha_vdim[9] = -1.414213562373095*A0[3]*dv10*dv2*dx11; 
  alpha_vdim[10] = 1.414213562373095*A1[3]*dv10*dv2*dx10; 
  alpha_vdim[13] = 1.414213562373095*A2[3]*dv10*dv3*dx10; 
  alpha_mid += fabs(0.0883883476483184*alpha_vdim[0]); 

  alpha_vdim[32] = dv11*(4.898979485566357*dx11*(A2[2]*wv3+A0[2]*wv1)-4.898979485566357*(A1[1]*dx10*wv1+phi[2]*dx11)); 
  alpha_vdim[33] = dv11*dx11*(4.898979485566357*(A2[3]*wv3+A0[3]*wv1)-4.898979485566357*phi[3]); 
  alpha_vdim[34] = -4.898979485566357*A1[3]*dv11*dx10*wv1; 
  alpha_vdim[35] = dv1*dv11*(1.414213562373095*A0[2]*dx11-1.414213562373095*A1[1]*dx10); 
  alpha_vdim[37] = 1.414213562373095*A2[2]*dv11*dv3*dx11; 
  alpha_vdim[39] = 1.414213562373095*A0[3]*dv1*dv11*dx11; 
  alpha_vdim[40] = -1.414213562373095*A1[3]*dv1*dv11*dx10; 
  alpha_vdim[44] = 1.414213562373095*A2[3]*dv11*dv3*dx11; 
  alpha_mid += fabs(0.0883883476483184*alpha_vdim[32]); 

  alpha_vdim[64] = -4.898979485566357*dv12*(A2[2]*dx11*wv2+A2[1]*dx10*wv1); 
  alpha_vdim[65] = -4.898979485566357*A2[3]*dv12*dx11*wv2; 
  alpha_vdim[66] = -4.898979485566357*A2[3]*dv12*dx10*wv1; 
  alpha_vdim[67] = -1.414213562373095*A2[1]*dv1*dv12*dx10; 
  alpha_vdim[68] = -1.414213562373095*A2[2]*dv12*dv2*dx11; 
  alpha_vdim[72] = -1.414213562373095*A2[3]*dv1*dv12*dx10; 
  alpha_vdim[73] = -1.414213562373095*A2[3]*dv12*dv2*dx11; 
  alpha_mid += fabs(0.0883883476483184*alpha_vdim[64]); 

  out[1] += 0.3061862178478971*(alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.3061862178478971*(f[4]*alpha_cdim[36]+f[0]*alpha_cdim[32]); 
  out[3] += 0.3061862178478971*(alpha_vdim[13]*f[13]+alpha_vdim[10]*f[10]+alpha_vdim[9]*f[9]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.3061862178478971*(f[12]*alpha_vdim[44]+f[8]*alpha_vdim[40]+f[7]*alpha_vdim[39]+f[5]*alpha_vdim[37]+f[3]*alpha_vdim[35]+f[2]*alpha_vdim[34]+f[1]*alpha_vdim[33]+f[0]*alpha_vdim[32]); 
  out[5] += 0.3061862178478971*(f[9]*alpha_vdim[73]+f[8]*alpha_vdim[72]+f[4]*alpha_vdim[68]+f[3]*alpha_vdim[67]+f[2]*alpha_vdim[66]+f[1]*alpha_vdim[65]+f[0]*alpha_vdim[64]); 
  out[6] += 0.3061862178478971*(f[9]*alpha_cdim[36]+f[1]*alpha_cdim[32]+alpha_cdim[3]*f[8]+alpha_cdim[0]*f[2]); 
  out[7] += 0.3061862178478971*(alpha_vdim[13]*f[20]+alpha_vdim[10]*f[17]+alpha_vdim[5]*f[12]+alpha_vdim[4]*f[9]+f[4]*alpha_vdim[9]+alpha_vdim[2]*f[6]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[8] += 0.3061862178478971*(f[11]*alpha_cdim[36]+f[3]*alpha_cdim[32]+alpha_vdim[9]*f[17]+alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13]+alpha_vdim[4]*f[10]+f[4]*alpha_vdim[10]+alpha_vdim[1]*f[6]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[9] += 0.3061862178478971*(f[5]*alpha_vdim[44]+f[16]*alpha_vdim[40]+f[3]*alpha_vdim[39]+f[12]*alpha_vdim[37]+f[7]*alpha_vdim[35]+f[6]*alpha_vdim[34]+f[0]*alpha_vdim[33]+f[1]*alpha_vdim[32]+alpha_cdim[3]*f[11]+alpha_cdim[0]*f[4]); 
  out[10] += 0.3061862178478971*(f[20]*alpha_vdim[44]+f[3]*alpha_vdim[40]+f[16]*alpha_vdim[39]+f[13]*alpha_vdim[37]+f[0]*alpha_cdim[36]+f[8]*alpha_vdim[35]+f[0]*alpha_vdim[34]+f[6]*alpha_vdim[33]+f[2]*alpha_vdim[32]+f[4]*alpha_cdim[32]); 
  out[11] += 0.3061862178478971*(f[21]*alpha_vdim[44]+f[2]*alpha_vdim[40]+f[1]*alpha_vdim[39]+f[14]*alpha_vdim[37]+f[0]*alpha_vdim[35]+f[8]*alpha_vdim[34]+f[7]*alpha_vdim[33]+f[3]*alpha_vdim[32]+alpha_vdim[13]*f[24]+alpha_vdim[5]*f[15]+alpha_vdim[2]*f[10]+f[2]*alpha_vdim[10]+alpha_vdim[1]*f[9]+f[1]*alpha_vdim[9]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]); 
  out[12] += 0.3061862178478971*(f[4]*alpha_vdim[73]+f[16]*alpha_vdim[72]+f[9]*alpha_vdim[68]+f[7]*alpha_vdim[67]+f[6]*alpha_vdim[66]+f[0]*alpha_vdim[65]+f[1]*alpha_vdim[64]+alpha_cdim[3]*f[14]+alpha_cdim[0]*f[5]); 
  out[13] += 0.3061862178478971*(f[17]*alpha_vdim[73]+f[3]*alpha_vdim[72]+f[10]*alpha_vdim[68]+f[8]*alpha_vdim[67]+f[0]*alpha_vdim[66]+f[6]*alpha_vdim[65]+f[2]*alpha_vdim[64]+f[15]*alpha_cdim[36]+f[5]*alpha_cdim[32]); 
  out[14] += 0.3061862178478971*(f[18]*alpha_vdim[73]+f[2]*alpha_vdim[72]+f[11]*alpha_vdim[68]+f[0]*alpha_vdim[67]+f[8]*alpha_vdim[66]+f[7]*alpha_vdim[65]+f[3]*alpha_vdim[64]+alpha_vdim[10]*f[24]+alpha_vdim[9]*f[23]+alpha_vdim[4]*f[15]+alpha_vdim[2]*f[13]+f[2]*alpha_vdim[13]+alpha_vdim[1]*f[12]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]); 
  out[15] += 0.3061862178478971*(f[1]*alpha_vdim[73]+f[19]*alpha_vdim[72]+f[0]*alpha_vdim[68]+f[11]*alpha_vdim[67]+f[10]*alpha_vdim[66]+f[9]*alpha_vdim[65]+f[4]*alpha_vdim[64]+f[1]*alpha_vdim[44]+f[22]*alpha_vdim[40]+f[21]*alpha_vdim[39]+f[0]*alpha_vdim[37]+f[14]*alpha_vdim[35]+f[13]*alpha_vdim[34]+f[12]*alpha_vdim[33]+f[5]*alpha_vdim[32]); 
  out[16] += 0.3061862178478971*(f[18]*alpha_cdim[36]+f[7]*alpha_cdim[32]+alpha_vdim[5]*f[20]+alpha_vdim[4]*f[17]+f[12]*alpha_vdim[13]+alpha_vdim[9]*f[10]+f[9]*alpha_vdim[10]+alpha_cdim[0]*f[8]+alpha_vdim[0]*f[6]+f[2]*(alpha_cdim[3]+alpha_vdim[1])+f[1]*alpha_vdim[2]); 
  out[17] += 0.3061862178478971*(f[13]*alpha_vdim[44]+f[7]*alpha_vdim[40]+f[8]*alpha_vdim[39]+f[20]*alpha_vdim[37]+f[1]*alpha_cdim[36]+f[16]*alpha_vdim[35]+f[1]*alpha_vdim[34]+f[2]*alpha_vdim[33]+f[6]*alpha_vdim[32]+f[9]*alpha_cdim[32]+alpha_cdim[3]*f[19]+alpha_cdim[0]*f[10]); 
  out[18] += 0.3061862178478971*(f[14]*alpha_vdim[44]+f[6]*alpha_vdim[40]+f[0]*alpha_vdim[39]+f[21]*alpha_vdim[37]+f[1]*alpha_vdim[35]+f[16]*alpha_vdim[34]+f[3]*alpha_vdim[33]+f[7]*alpha_vdim[32]+alpha_vdim[13]*f[28]+alpha_vdim[5]*f[23]+alpha_vdim[2]*f[17]+alpha_cdim[0]*f[11]+f[6]*alpha_vdim[10]+alpha_vdim[0]*f[9]+f[0]*alpha_vdim[9]+(alpha_cdim[3]+alpha_vdim[1])*f[4]+f[1]*alpha_vdim[4]); 
  out[19] += 0.3061862178478971*(f[27]*alpha_vdim[44]+f[0]*alpha_vdim[40]+f[6]*alpha_vdim[39]+f[22]*alpha_vdim[37]+f[3]*alpha_cdim[36]+f[2]*alpha_vdim[35]+f[3]*alpha_vdim[34]+f[16]*alpha_vdim[33]+f[8]*alpha_vdim[32]+f[11]*alpha_cdim[32]+alpha_vdim[5]*f[24]+alpha_vdim[1]*f[17]+alpha_vdim[13]*f[15]+alpha_vdim[0]*f[10]+f[0]*alpha_vdim[10]+f[6]*alpha_vdim[9]+alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]); 
  out[20] += 0.3061862178478971*(f[10]*alpha_vdim[73]+f[7]*alpha_vdim[72]+f[17]*alpha_vdim[68]+f[16]*alpha_vdim[67]+f[1]*alpha_vdim[66]+f[2]*alpha_vdim[65]+f[6]*alpha_vdim[64]+f[23]*alpha_cdim[36]+f[12]*alpha_cdim[32]+alpha_cdim[3]*f[22]+alpha_cdim[0]*f[13]); 
  out[21] += 0.3061862178478971*(f[11]*alpha_vdim[73]+f[6]*alpha_vdim[72]+f[18]*alpha_vdim[68]+f[1]*alpha_vdim[67]+f[16]*alpha_vdim[66]+f[3]*alpha_vdim[65]+f[7]*alpha_vdim[64]+alpha_vdim[10]*f[28]+alpha_vdim[4]*f[23]+alpha_vdim[2]*f[20]+alpha_vdim[9]*f[15]+alpha_cdim[0]*f[14]+f[6]*alpha_vdim[13]+alpha_vdim[0]*f[12]+(alpha_cdim[3]+alpha_vdim[1])*f[5]+f[1]*alpha_vdim[5]); 
  out[22] += 0.3061862178478971*(f[26]*alpha_vdim[73]+f[0]*alpha_vdim[72]+f[19]*alpha_vdim[68]+f[2]*alpha_vdim[67]+f[3]*alpha_vdim[66]+f[16]*alpha_vdim[65]+f[8]*alpha_vdim[64]+f[25]*alpha_cdim[36]+f[14]*alpha_cdim[32]+alpha_vdim[9]*f[28]+alpha_vdim[4]*f[24]+alpha_vdim[1]*f[20]+alpha_vdim[10]*f[15]+alpha_vdim[0]*f[13]+f[0]*alpha_vdim[13]+alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]); 
  out[23] += 0.3061862178478971*(f[0]*alpha_vdim[73]+f[26]*alpha_vdim[72]+f[1]*alpha_vdim[68]+f[18]*alpha_vdim[67]+f[17]*alpha_vdim[66]+f[4]*alpha_vdim[65]+f[9]*alpha_vdim[64]+f[0]*alpha_vdim[44]+f[27]*alpha_vdim[40]+f[14]*alpha_vdim[39]+f[1]*alpha_vdim[37]+f[21]*alpha_vdim[35]+f[20]*alpha_vdim[34]+f[5]*alpha_vdim[33]+f[12]*alpha_vdim[32]+alpha_cdim[3]*f[25]+alpha_cdim[0]*f[15]); 
  out[24] += 0.3061862178478971*(f[6]*alpha_vdim[73]+f[11]*alpha_vdim[72]+f[2]*alpha_vdim[68]+f[19]*alpha_vdim[67]+f[4]*alpha_vdim[66]+f[17]*alpha_vdim[65]+f[10]*alpha_vdim[64]+f[6]*alpha_vdim[44]+f[14]*alpha_vdim[40]+f[27]*alpha_vdim[39]+f[2]*alpha_vdim[37]+f[5]*alpha_cdim[36]+f[22]*alpha_vdim[35]+f[5]*alpha_vdim[34]+f[20]*alpha_vdim[33]+f[13]*alpha_vdim[32]+f[15]*alpha_cdim[32]); 
  out[25] += 0.3061862178478971*(f[7]*alpha_vdim[73]+f[10]*alpha_vdim[72]+f[3]*alpha_vdim[68]+f[4]*alpha_vdim[67]+f[19]*alpha_vdim[66]+f[18]*alpha_vdim[65]+f[11]*alpha_vdim[64]+f[7]*alpha_vdim[44]+f[13]*alpha_vdim[40]+f[12]*alpha_vdim[39]+f[3]*alpha_vdim[37]+f[5]*alpha_vdim[35]+f[22]*alpha_vdim[34]+f[21]*alpha_vdim[33]+f[14]*alpha_vdim[32]+alpha_vdim[2]*f[24]+alpha_vdim[1]*f[23]+alpha_vdim[0]*f[15]+alpha_vdim[10]*f[13]+f[10]*alpha_vdim[13]+alpha_vdim[9]*f[12]+alpha_vdim[4]*f[5]+f[4]*alpha_vdim[5]); 
  out[26] += 0.3061862178478971*(f[22]*alpha_vdim[44]+f[1]*alpha_vdim[40]+f[2]*alpha_vdim[39]+f[27]*alpha_vdim[37]+f[7]*alpha_cdim[36]+f[6]*alpha_vdim[35]+f[7]*alpha_vdim[34]+f[8]*alpha_vdim[33]+f[16]*alpha_vdim[32]+f[18]*alpha_cdim[32]+alpha_vdim[5]*f[28]+alpha_vdim[13]*f[23]+alpha_cdim[0]*f[19]+alpha_vdim[0]*f[17]+(alpha_cdim[3]+alpha_vdim[1])*f[10]+f[1]*alpha_vdim[10]+alpha_vdim[2]*f[9]+f[2]*alpha_vdim[9]+alpha_vdim[4]*f[6]); 
  out[27] += 0.3061862178478971*(f[19]*alpha_vdim[73]+f[1]*alpha_vdim[72]+f[26]*alpha_vdim[68]+f[6]*alpha_vdim[67]+f[7]*alpha_vdim[66]+f[8]*alpha_vdim[65]+f[16]*alpha_vdim[64]+f[29]*alpha_cdim[36]+f[21]*alpha_cdim[32]+alpha_vdim[4]*f[28]+alpha_vdim[9]*f[24]+alpha_vdim[10]*f[23]+alpha_cdim[0]*f[22]+alpha_vdim[0]*f[20]+(alpha_cdim[3]+alpha_vdim[1])*f[13]+f[1]*alpha_vdim[13]+alpha_vdim[2]*f[12]+alpha_vdim[5]*f[6]); 
  out[28] += 0.3061862178478971*(f[2]*alpha_vdim[73]+f[18]*alpha_vdim[72]+f[6]*alpha_vdim[68]+f[26]*alpha_vdim[67]+f[9]*alpha_vdim[66]+f[10]*alpha_vdim[65]+f[17]*alpha_vdim[64]+f[2]*alpha_vdim[44]+f[21]*alpha_vdim[40]+f[22]*alpha_vdim[39]+f[6]*alpha_vdim[37]+f[12]*alpha_cdim[36]+f[27]*alpha_vdim[35]+f[12]*alpha_vdim[34]+f[13]*alpha_vdim[33]+f[20]*alpha_vdim[32]+f[23]*alpha_cdim[32]+alpha_cdim[3]*f[30]+alpha_cdim[0]*f[24]); 
  out[29] += 0.3061862178478971*(f[3]*alpha_vdim[73]+f[17]*alpha_vdim[72]+f[7]*alpha_vdim[68]+f[9]*alpha_vdim[67]+f[26]*alpha_vdim[66]+f[11]*alpha_vdim[65]+f[18]*alpha_vdim[64]+f[3]*alpha_vdim[44]+f[20]*alpha_vdim[40]+f[5]*alpha_vdim[39]+f[7]*alpha_vdim[37]+f[12]*alpha_vdim[35]+f[27]*alpha_vdim[34]+f[14]*alpha_vdim[33]+f[21]*alpha_vdim[32]+alpha_vdim[2]*f[28]+alpha_cdim[0]*f[25]+alpha_vdim[0]*f[23]+alpha_vdim[10]*f[20]+alpha_vdim[13]*f[17]+(alpha_cdim[3]+alpha_vdim[1])*f[15]+alpha_vdim[4]*f[12]+alpha_vdim[5]*f[9]+f[5]*alpha_vdim[9]); 
  out[30] += 0.3061862178478971*(f[16]*alpha_vdim[73]+f[4]*alpha_vdim[72]+f[8]*alpha_vdim[68]+f[10]*alpha_vdim[67]+f[11]*alpha_vdim[66]+f[26]*alpha_vdim[65]+f[19]*alpha_vdim[64]+f[16]*alpha_vdim[44]+f[5]*alpha_vdim[40]+f[20]*alpha_vdim[39]+f[8]*alpha_vdim[37]+f[14]*alpha_cdim[36]+f[13]*alpha_vdim[35]+f[14]*alpha_vdim[34]+f[27]*alpha_vdim[33]+f[22]*alpha_vdim[32]+f[25]*alpha_cdim[32]+alpha_vdim[1]*f[28]+alpha_vdim[0]*f[24]+alpha_vdim[9]*f[20]+alpha_vdim[2]*f[15]+alpha_vdim[4]*f[13]+f[4]*alpha_vdim[13]+alpha_vdim[5]*f[10]+f[5]*alpha_vdim[10]); 
  out[31] += 0.3061862178478971*(f[8]*alpha_vdim[73]+f[9]*alpha_vdim[72]+f[16]*alpha_vdim[68]+f[17]*alpha_vdim[67]+f[18]*alpha_vdim[66]+f[19]*alpha_vdim[65]+f[26]*alpha_vdim[64]+f[8]*alpha_vdim[44]+f[12]*alpha_vdim[40]+f[13]*alpha_vdim[39]+f[16]*alpha_vdim[37]+f[21]*alpha_cdim[36]+f[20]*alpha_vdim[35]+f[21]*alpha_vdim[34]+f[22]*alpha_vdim[33]+f[27]*alpha_vdim[32]+f[29]*alpha_cdim[32]+alpha_cdim[0]*f[30]+alpha_vdim[0]*f[28]+(alpha_cdim[3]+alpha_vdim[1])*f[24]+alpha_vdim[2]*f[23]+alpha_vdim[4]*f[20]+alpha_vdim[5]*f[17]+alpha_vdim[9]*f[13]+f[9]*alpha_vdim[13]+alpha_vdim[10]*f[12]); 

  return alpha_mid; 
} 

