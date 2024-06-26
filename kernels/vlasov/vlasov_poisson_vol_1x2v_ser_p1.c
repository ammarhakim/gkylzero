#include <gkyl_vlasov_kernels.h> 

GKYL_CU_DH double vlasov_poisson_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:     potential (scaled by appropriate factors).
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dx10 = 2/dxv[0]; 
  const double *phi = &field[0]; 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[16]; 
  double alpha_vdim[32]; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  cflFreq_mid += 3.0*(fabs(w0dx0)+0.5*dv0dx0); 

  alpha_vdim[0] = -3.464101615137754*phi[1]*dv10*dx10; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[0]); 

  cflFreq_mid += 5.0*fabs(0.0); 

  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*alpha_vdim[0]*f[0]; 
  out[4] += 0.5477225575051661*alpha_cdim[2]*f[8]+0.6123724356957944*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]); 
  out[5] += 0.6123724356957944*(alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*alpha_vdim[0]*f[3]; 
  out[7] += 0.5477225575051661*alpha_cdim[2]*f[10]+0.6123724356957944*(alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+alpha_cdim[2]*f[3]); 
  out[8] += 1.369306393762915*alpha_vdim[0]*f[2]; 
  out[9] += 0.6123724356957944*alpha_cdim[0]*f[8]+1.369306393762915*alpha_vdim[0]*f[4]+0.5477225575051661*alpha_cdim[2]*f[2]; 
  out[10] += 1.369306393762915*alpha_vdim[0]*f[6]; 
  out[11] += 0.6123724356957944*alpha_cdim[0]*f[10]+1.369306393762915*alpha_vdim[0]*f[7]+0.5477225575051661*alpha_cdim[2]*f[6]; 
  out[13] += 0.6123724356957944*(alpha_cdim[2]*f[14]+alpha_cdim[0]*f[12]); 
  out[14] += 0.6123724356957944*alpha_vdim[0]*f[12]; 
  out[15] += 0.6123724356957944*(alpha_cdim[0]*f[14]+alpha_vdim[0]*f[13]+alpha_cdim[2]*f[12]); 

  return cflFreq_mid; 
} 

