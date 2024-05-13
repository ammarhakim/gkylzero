#include <gkyl_vlasov_poisson_kernels.h> 

GKYL_CU_DH double vlasov_poisson_extem_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:     potentials, including external (scaled by appropriate factors).
  // f:         Input distribution function.
  // out:       Incremented output.

  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *phi = &field[0]; 

  const double *A0 = &field[3]; 
  const double *A1 = &field[6]; 

  double cflFreq_mid = 0.0; 

  double alpha_cdim[20]; 
  double alpha_vdim[40]; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  cflFreq_mid += 5.0*(fabs(w0dx0)+0.5*dv0dx0); 

  alpha_vdim[0] = dv10*dx10*(3.464101615137754*A1[1]*wv2-3.464101615137754*phi[1]); 
  alpha_vdim[1] = dv10*dx10*(7.745966692414834*A1[2]*wv2-7.745966692414834*phi[2]); 
  alpha_vdim[3] = A1[1]*dv10*dv2*dx10; 
  alpha_vdim[5] = 2.23606797749979*A1[2]*dv10*dv2*dx10; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[0]); 

  alpha_vdim[20] = -3.464101615137754*A1[1]*dv11*dx10*wv1; 
  alpha_vdim[21] = -7.745966692414834*A1[2]*dv11*dx10*wv1; 
  alpha_vdim[22] = -1.0*A1[1]*dv1*dv11*dx10; 
  alpha_vdim[24] = -2.23606797749979*A1[2]*dv1*dv11*dx10; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[0]); 

  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.5477225575051661*(alpha_vdim[5]*f[13]+alpha_cdim[2]*f[8]+alpha_vdim[1]*f[7])+0.6123724356957944*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.5477225575051661*(alpha_vdim[4]*f[11]+alpha_vdim[1]*f[7])+0.6123724356957944*(alpha_cdim[2]*f[6]+alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]+alpha_cdim[0]*f[3]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.5477225575051661*(alpha_vdim[5]*f[15]+alpha_vdim[4]*f[12]+alpha_vdim[3]*f[9]+alpha_vdim[2]*f[8])+0.6123724356957944*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[7] += 1.369306393762915*(alpha_cdim[2]*f[4]+alpha_cdim[0]*f[1]); 
  out[8] += 1.369306393762915*(alpha_vdim[5]*f[10]+alpha_vdim[3]*f[6]+alpha_vdim[1]*f[4]+alpha_vdim[0]*f[2]); 
  out[9] += 1.369306393762915*(alpha_vdim[4]*f[10]+alpha_vdim[2]*f[6]+alpha_vdim[1]*f[5]+alpha_vdim[0]*f[3]); 
  out[10] += 0.5477225575051661*(alpha_vdim[3]*f[15]+alpha_cdim[2]*f[14]+alpha_vdim[1]*f[13]+alpha_vdim[2]*f[12]+alpha_vdim[1]*f[11]+alpha_vdim[5]*f[9]+alpha_vdim[4]*f[8]+(alpha_vdim[5]+alpha_vdim[4])*f[7])+0.6123724356957944*(alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]+alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 
  out[11] += 0.6123724356957944*alpha_vdim[3]*f[13]+1.224744871391589*alpha_cdim[2]*f[12]+0.6123724356957944*alpha_vdim[0]*f[7]+0.5477225575051661*alpha_vdim[5]*f[5]+1.369306393762915*alpha_cdim[0]*f[4]+f[1]*(1.369306393762915*alpha_cdim[2]+0.5477225575051661*alpha_vdim[1]); 
  out[12] += 1.224744871391589*(alpha_vdim[5]*f[17]+alpha_vdim[1]*f[11])+1.369306393762915*alpha_vdim[3]*f[10]+0.6123724356957944*alpha_cdim[0]*f[8]+1.369306393762915*(alpha_vdim[5]*f[6]+alpha_vdim[0]*f[4])+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[2]; 
  out[13] += 0.6123724356957944*alpha_vdim[2]*f[11]+1.369306393762915*alpha_cdim[2]*f[10]+0.6123724356957944*alpha_vdim[0]*f[7]+1.369306393762915*alpha_cdim[0]*f[5]+0.5477225575051661*(alpha_vdim[4]*f[4]+alpha_vdim[1]*f[1]); 
  out[14] += 1.224744871391589*(alpha_vdim[5]*f[19]+alpha_vdim[3]*f[16])+alpha_vdim[1]*(0.6123724356957944*f[12]+1.369306393762915*f[10])+alpha_vdim[0]*(0.6123724356957944*f[8]+1.369306393762915*f[6])+f[4]*(1.369306393762915*alpha_vdim[5]+0.5477225575051661*alpha_vdim[4])+f[2]*(1.369306393762915*alpha_vdim[3]+0.5477225575051661*alpha_vdim[2]); 
  out[15] += 1.224744871391589*alpha_vdim[4]*f[17]+0.6123724356957944*alpha_cdim[2]*f[16]+1.224744871391589*alpha_vdim[1]*f[13]+1.369306393762915*alpha_vdim[2]*f[10]+0.6123724356957944*alpha_cdim[0]*f[9]+1.369306393762915*(alpha_vdim[4]*f[6]+alpha_vdim[0]*f[5]+alpha_vdim[1]*f[3]); 
  out[16] += 1.224744871391589*alpha_vdim[4]*f[18]+0.6123724356957944*alpha_vdim[1]*f[15]+1.224744871391589*alpha_vdim[2]*f[14]+1.369306393762915*alpha_vdim[1]*f[10]+alpha_vdim[0]*(0.6123724356957944*f[9]+1.369306393762915*f[6])+(0.5477225575051661*alpha_vdim[5]+1.369306393762915*alpha_vdim[4])*f[5]+(0.5477225575051661*alpha_vdim[3]+1.369306393762915*alpha_vdim[2])*f[3]; 
  out[17] += 1.224744871391589*alpha_cdim[2]*f[18]+0.4898979485566357*alpha_vdim[5]*f[15]+0.6123724356957944*alpha_vdim[0]*f[13]+0.4898979485566357*alpha_vdim[4]*f[12]+0.6123724356957944*alpha_vdim[0]*f[11]+1.369306393762915*alpha_cdim[0]*f[10]+0.6123724356957944*(alpha_vdim[3]+alpha_vdim[2])*f[7]+1.369306393762915*alpha_cdim[2]*f[5]+0.5477225575051661*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4]); 
  out[18] += 1.224744871391589*(alpha_vdim[3]*f[19]+alpha_vdim[1]*f[17]+alpha_vdim[5]*f[16])+0.6123724356957944*(alpha_cdim[0]*f[14]+alpha_vdim[0]*f[12])+(1.224744871391589*alpha_vdim[5]+0.4898979485566357*alpha_vdim[4])*f[11]+1.369306393762915*alpha_vdim[0]*f[10]+0.6123724356957944*alpha_vdim[1]*f[8]+0.5477225575051661*alpha_cdim[2]*f[6]+1.369306393762915*(alpha_vdim[1]*f[6]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4])+0.5477225575051661*(alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]); 
  out[19] += 1.224744871391589*(alpha_vdim[2]*f[18]+alpha_vdim[1]*f[17])+0.6123724356957944*(alpha_cdim[0]*f[16]+alpha_vdim[0]*f[15])+1.224744871391589*alpha_vdim[4]*f[14]+(0.4898979485566357*alpha_vdim[5]+1.224744871391589*alpha_vdim[4])*f[13]+1.369306393762915*alpha_vdim[0]*f[10]+0.6123724356957944*alpha_cdim[2]*f[9]+alpha_vdim[1]*(0.6123724356957944*f[9]+1.369306393762915*f[6])+(0.5477225575051661*alpha_vdim[3]+1.369306393762915*alpha_vdim[2])*f[5]+f[3]*(0.5477225575051661*alpha_vdim[5]+1.369306393762915*alpha_vdim[4]); 

  return cflFreq_mid; 
} 

