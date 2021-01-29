#include <gkyl_vlasov_kernels.h> 
gkyl_real vlasov_vol_1x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *qmem, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields.
  // f:         Input distribution function.
  // out:       Incremented output.
  gkyl_real dv0dx0 = dxv[1]/dxv[0]; 
  gkyl_real w0dx0 = w[1]/dxv[0]; 
  const gkyl_real dv10 = 2/dxv[1]; 
  const gkyl_real *E0 = &qmem[0]; 
  const gkyl_real dv1 = dxv[1], wv1 = w[1]; 
  const gkyl_real dv11 = 2/dxv[2]; 
  const gkyl_real *E1 = &qmem[2]; 
  const gkyl_real dv2 = dxv[2], wv2 = w[2]; 

  const gkyl_real *B2 = &qmem[10]; 
  gkyl_real alpha_mid = 0.0; 
  gkyl_real alpha_cdim[8]; 
  gkyl_real alpha_vdim[16]; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  alpha_mid += fabs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = 2.0*dv10*(B2[0]*wv2+E0[0]); 
  alpha_vdim[1] = 2.0*dv10*(B2[1]*wv2+E0[1]); 
  alpha_vdim[3] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha_vdim[5] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha_mid += fabs(0.1767766952966368*alpha_vdim[0]); 

  alpha_vdim[8] = dv11*(2.0*E1[0]-2.0*B2[0]*wv1); 
  alpha_vdim[9] = dv11*(2.0*E1[1]-2.0*B2[1]*wv1); 
  alpha_vdim[10] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha_vdim[12] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha_mid += fabs(0.1767766952966368*alpha_vdim[8]); 

  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[4]*alpha_vdim[12]+f[2]*alpha_vdim[10]+f[1]*alpha_vdim[9]+f[0]*alpha_vdim[8]); 
  out[4] += 0.6123724356957944*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.6123724356957944*(f[2]*alpha_vdim[12]+f[4]*alpha_vdim[10]+f[0]*alpha_vdim[9]+f[1]*alpha_vdim[8]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*(f[1]*alpha_vdim[12]+f[0]*alpha_vdim[10]+f[4]*alpha_vdim[9]+f[2]*alpha_vdim[8]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 0.6123724356957944*(f[0]*alpha_vdim[12]+f[1]*alpha_vdim[10]+f[2]*alpha_vdim[9]+f[4]*alpha_vdim[8]+alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]); 

  return alpha_mid; 
} 
