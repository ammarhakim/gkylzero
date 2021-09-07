#include <gkyl_vlasov_kernels.h> 

GKYL_CU_DH double vlasov_poisson_extem_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *f, double* GKYL_RESTRICT out) 
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
  const double *A1 = &vecA[2]; 
  double alpha_mid = 0.0; 
  double alpha_cdim[8]; 
  double alpha_vdim[16]; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  alpha_mid += fabs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = dv10*dx10*(3.464101615137754*A1[1]*wv2-3.464101615137754*phi[1]); 
  alpha_vdim[3] = A1[1]*dv10*dv2*dx10; 
  alpha_mid += fabs(0.1767766952966368*alpha_vdim[0]); 

  alpha_vdim[8] = -3.464101615137754*A1[1]*dv11*dx10*wv1; 
  alpha_vdim[10] = -1.0*A1[1]*dv1*dv11*dx10; 
  alpha_mid += fabs(0.1767766952966368*alpha_vdim[8]); 

  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[3]*f[3]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[2]*alpha_vdim[10]+f[0]*alpha_vdim[8]); 
  out[4] += 0.6123724356957944*(alpha_vdim[3]*f[5]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]); 
  out[5] += 0.6123724356957944*(f[4]*alpha_vdim[10]+f[1]*alpha_vdim[8]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*(f[0]*alpha_vdim[10]+f[2]*alpha_vdim[8]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 0.6123724356957944*(f[1]*alpha_vdim[10]+f[4]*alpha_vdim[8]+alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+alpha_cdim[2]*f[3]+f[1]*alpha_vdim[3]); 

  return alpha_mid; 
} 

