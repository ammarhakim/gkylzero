#include <gkyl_vlasov_kernels.h> 

GKYL_CU_DH double vlasov_poisson_extem_vol_1x2v_tensor_p2(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *f, double* GKYL_RESTRICT out) 
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

  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[3]; 
  double alpha_mid = 0.0; 
  double alpha_cdim[27]; 
  double alpha_vdim[54]; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  alpha_mid += fabs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = dv10*dx10*(3.464101615137754*A1[1]*wv2-3.464101615137754*phi[1]); 
  alpha_vdim[1] = dv10*dx10*(7.745966692414834*A1[2]*wv2-7.745966692414834*phi[2]); 
  alpha_vdim[3] = A1[1]*dv10*dv2*dx10; 
  alpha_vdim[5] = 2.23606797749979*A1[2]*dv10*dv2*dx10; 
  alpha_mid += fabs(0.1767766952966368*alpha_vdim[0]); 

  alpha_vdim[27] = -3.464101615137754*A1[1]*dv11*dx10*wv1; 
  alpha_vdim[28] = -7.745966692414834*A1[2]*dv11*dx10*wv1; 
  alpha_vdim[29] = -1.0*A1[1]*dv1*dv11*dx10; 
  alpha_vdim[31] = -2.23606797749979*A1[2]*dv1*dv11*dx10; 
  alpha_mid += fabs(0.1767766952966368*alpha_vdim[27]); 

  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[4]*alpha_vdim[31]+f[2]*alpha_vdim[29]+f[1]*alpha_vdim[28]+f[0]*alpha_vdim[27]); 
  out[4] += 0.5477225575051661*(alpha_vdim[5]*f[13]+alpha_cdim[2]*f[8]+alpha_vdim[1]*f[7])+0.6123724356957944*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.5477225575051661*f[11]*alpha_vdim[31]+0.6123724356957944*(f[2]*alpha_vdim[31]+f[4]*alpha_vdim[29])+0.5477225575051661*f[7]*alpha_vdim[28]+0.6123724356957944*(f[0]*alpha_vdim[28]+f[1]*alpha_vdim[27]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += (0.5477225575051661*f[12]+0.6123724356957944*f[1])*alpha_vdim[31]+0.5477225575051661*f[8]*alpha_vdim[29]+0.6123724356957944*(f[0]*alpha_vdim[29]+f[4]*alpha_vdim[28]+f[2]*alpha_vdim[27])+0.5477225575051661*(alpha_vdim[5]*f[15]+alpha_vdim[3]*f[9])+0.6123724356957944*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 1.369306393762915*(alpha_cdim[2]*f[4]+alpha_cdim[0]*f[1]); 
  out[8] += 1.369306393762915*(alpha_vdim[5]*f[10]+alpha_vdim[3]*f[6]+alpha_vdim[1]*f[4]+alpha_vdim[0]*f[2]); 
  out[9] += 1.369306393762915*(f[10]*alpha_vdim[31]+f[6]*alpha_vdim[29]+f[5]*alpha_vdim[28]+f[3]*alpha_vdim[27]); 
  out[10] += (0.4898979485566357*f[20]+0.5477225575051661*(f[8]+f[7])+0.6123724356957944*f[0])*alpha_vdim[31]+(0.5477225575051661*f[12]+0.6123724356957944*f[1])*alpha_vdim[29]+0.5477225575051661*f[11]*alpha_vdim[28]+0.6123724356957944*(f[2]*alpha_vdim[28]+f[4]*alpha_vdim[27])+0.4898979485566357*alpha_vdim[5]*f[21]+0.5477225575051661*(alpha_vdim[3]*f[15]+alpha_cdim[2]*f[14]+alpha_vdim[1]*f[13]+alpha_vdim[5]*(f[9]+f[7]))+0.6123724356957944*(alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]); 
  out[11] += 0.6123724356957944*alpha_vdim[3]*f[13]+1.224744871391589*alpha_cdim[2]*f[12]+0.6123724356957944*alpha_vdim[0]*f[7]+0.5477225575051661*alpha_vdim[5]*f[5]+1.369306393762915*alpha_cdim[0]*f[4]+f[1]*(1.369306393762915*alpha_cdim[2]+0.5477225575051661*alpha_vdim[1]); 
  out[12] += 1.224744871391589*(alpha_vdim[5]*f[17]+alpha_vdim[1]*f[11])+1.369306393762915*alpha_vdim[3]*f[10]+0.6123724356957944*alpha_cdim[0]*f[8]+1.369306393762915*(alpha_vdim[5]*f[6]+alpha_vdim[0]*f[4])+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[2]; 
  out[13] += 0.5477225575051661*f[4]*alpha_vdim[31]+0.6123724356957944*f[11]*alpha_vdim[29]+0.5477225575051661*f[1]*alpha_vdim[28]+0.6123724356957944*f[7]*alpha_vdim[27]+1.369306393762915*(alpha_cdim[2]*f[10]+alpha_cdim[0]*f[5]); 
  out[14] += 0.5477225575051661*(f[4]*alpha_vdim[31]+f[2]*alpha_vdim[29])+0.6123724356957944*(f[12]*alpha_vdim[28]+f[8]*alpha_vdim[27])+1.224744871391589*(alpha_vdim[5]*f[19]+alpha_vdim[3]*f[16])+1.369306393762915*(alpha_vdim[1]*f[10]+alpha_vdim[0]*f[6]+f[4]*alpha_vdim[5]+f[2]*alpha_vdim[3]); 
  out[15] += 1.224744871391589*f[17]*alpha_vdim[31]+1.369306393762915*(f[6]*alpha_vdim[31]+f[10]*alpha_vdim[29])+1.224744871391589*f[13]*alpha_vdim[28]+1.369306393762915*(f[3]*alpha_vdim[28]+f[5]*alpha_vdim[27])+0.6123724356957944*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[9]); 
  out[16] += (1.224744871391589*f[18]+1.369306393762915*f[5])*alpha_vdim[31]+1.224744871391589*f[14]*alpha_vdim[29]+1.369306393762915*(f[3]*alpha_vdim[29]+f[10]*alpha_vdim[28]+f[6]*alpha_vdim[27])+0.6123724356957944*(alpha_vdim[1]*f[15]+alpha_vdim[0]*f[9])+0.5477225575051661*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]); 
  out[17] += (0.4898979485566357*f[12]+0.5477225575051661*f[1])*alpha_vdim[31]+(0.5477225575051661*f[20]+0.6123724356957944*f[7])*alpha_vdim[29]+0.5477225575051661*f[4]*alpha_vdim[28]+0.6123724356957944*f[11]*alpha_vdim[27]+0.5477225575051661*alpha_vdim[3]*f[21]+1.224744871391589*alpha_cdim[2]*f[18]+0.4898979485566357*alpha_vdim[5]*f[15]+0.6123724356957944*alpha_vdim[0]*f[13]+1.369306393762915*alpha_cdim[0]*f[10]+0.6123724356957944*alpha_vdim[3]*f[7]+1.369306393762915*alpha_cdim[2]*f[5]+0.5477225575051661*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]); 
  out[18] += 0.4898979485566357*f[11]*alpha_vdim[31]+0.5477225575051661*(f[2]*alpha_vdim[31]+f[4]*alpha_vdim[29]+f[20]*alpha_vdim[28])+0.6123724356957944*(f[8]*alpha_vdim[28]+f[12]*alpha_vdim[27])+1.095445115010332*alpha_vdim[5]*f[24]+1.224744871391589*(alpha_vdim[3]*f[19]+alpha_vdim[1]*f[17]+alpha_vdim[5]*f[16])+0.6123724356957944*alpha_cdim[0]*f[14]+1.224744871391589*alpha_vdim[5]*f[11]+1.369306393762915*alpha_vdim[0]*f[10]+0.5477225575051661*alpha_cdim[2]*f[6]+1.369306393762915*(alpha_vdim[1]*f[6]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]); 
  out[19] += (1.095445115010332*f[23]+1.224744871391589*(f[14]+f[13])+1.369306393762915*f[3])*alpha_vdim[31]+(1.224744871391589*f[18]+1.369306393762915*f[5])*alpha_vdim[29]+1.224744871391589*f[17]*alpha_vdim[28]+1.369306393762915*(f[6]*alpha_vdim[28]+f[10]*alpha_vdim[27])+0.5477225575051661*(alpha_cdim[2]*f[22]+alpha_vdim[1]*f[21])+0.6123724356957944*(alpha_cdim[0]*f[16]+alpha_vdim[0]*f[15])+0.4898979485566357*alpha_vdim[5]*f[13]+0.6123724356957944*(alpha_cdim[2]+alpha_vdim[1])*f[9]+0.5477225575051661*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]); 
  out[20] += 1.369306393762915*(alpha_vdim[3]*f[17]+alpha_cdim[0]*f[12]+alpha_vdim[0]*f[11])+1.224744871391589*(alpha_vdim[5]*f[10]+(alpha_cdim[2]+alpha_vdim[1])*f[4]); 
  out[21] += 1.224744871391589*f[10]*alpha_vdim[31]+1.369306393762915*f[17]*alpha_vdim[29]+1.224744871391589*f[5]*alpha_vdim[28]+1.369306393762915*(f[13]*alpha_vdim[27]+alpha_cdim[2]*f[19]+alpha_cdim[0]*f[15]); 
  out[22] += 1.224744871391589*(f[10]*alpha_vdim[31]+f[6]*alpha_vdim[29])+1.369306393762915*(f[18]*alpha_vdim[28]+f[14]*alpha_vdim[27]+alpha_vdim[1]*f[19]+alpha_vdim[0]*f[16])+1.224744871391589*(alpha_vdim[5]*f[10]+alpha_vdim[3]*f[6]); 
  out[23] += 0.4898979485566357*f[4]*alpha_vdim[31]+0.5477225575051661*(f[11]*alpha_vdim[29]+f[12]*alpha_vdim[28])+0.6123724356957944*f[20]*alpha_vdim[27]+1.224744871391589*alpha_vdim[3]*f[24]+1.095445115010332*alpha_vdim[5]*f[19]+1.369306393762915*(alpha_cdim[0]*f[18]+alpha_vdim[0]*f[17]+alpha_vdim[3]*f[11])+1.224744871391589*((alpha_cdim[2]+alpha_vdim[1])*f[10]+f[4]*alpha_vdim[5]); 
  out[24] += (1.095445115010332*f[18]+1.224744871391589*f[5])*alpha_vdim[31]+(1.224744871391589*f[23]+1.369306393762915*f[13])*alpha_vdim[29]+1.224744871391589*f[10]*alpha_vdim[28]+1.369306393762915*f[17]*alpha_vdim[27]+1.224744871391589*alpha_cdim[2]*f[25]+0.6123724356957944*alpha_vdim[0]*f[21]+1.369306393762915*(alpha_cdim[0]*f[19]+alpha_cdim[2]*f[15])+0.5477225575051661*(alpha_vdim[1]*f[15]+alpha_vdim[3]*f[13])+0.4898979485566357*alpha_vdim[5]*f[5]; 
  out[25] += 1.095445115010332*f[17]*alpha_vdim[31]+1.224744871391589*(f[6]*alpha_vdim[31]+f[10]*alpha_vdim[29]+f[23]*alpha_vdim[28])+1.369306393762915*(f[14]*alpha_vdim[28]+f[18]*alpha_vdim[27])+1.224744871391589*alpha_vdim[1]*f[24]+0.6123724356957944*alpha_cdim[0]*f[22]+1.369306393762915*alpha_vdim[0]*f[19]+1.095445115010332*alpha_vdim[5]*f[17]+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[16]+1.224744871391589*(alpha_vdim[3]*f[10]+alpha_vdim[5]*f[6]); 
  out[26] += 1.095445115010332*f[10]*alpha_vdim[31]+1.224744871391589*(f[17]*alpha_vdim[29]+f[18]*alpha_vdim[28])+1.369306393762915*(f[23]*alpha_vdim[27]+alpha_cdim[0]*f[25]+alpha_vdim[0]*f[24])+1.224744871391589*((alpha_cdim[2]+alpha_vdim[1])*f[19]+alpha_vdim[3]*f[17])+1.095445115010332*alpha_vdim[5]*f[10]; 

  return alpha_mid; 
} 

