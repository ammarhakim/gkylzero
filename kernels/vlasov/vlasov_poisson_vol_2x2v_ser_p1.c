#include <gkyl_vlasov_kernels.h> 

GKYL_CU_DH double vlasov_poisson_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *f, double* GKYL_RESTRICT out) 
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

  double alpha_mid = 0.0; 
  double alpha_cdim[32]; 
  double alpha_vdim[32]; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[3] = 2.309401076758503*dv0dx0; 
  alpha_mid += fabs(w0dx0)+0.5*dv0dx0; 

  alpha_cdim[16] = 8.0*w1dx1; 
  alpha_cdim[20] = 2.309401076758503*dv1dx1; 
  alpha_mid += fabs(w1dx1)+0.5*dv1dx1; 

  alpha_vdim[0] = -3.464101615137754*phi[1]*dv10*dx10; 
  alpha_vdim[2] = -3.464101615137754*phi[3]*dv10*dx10; 
  alpha_mid += fabs(0.125*alpha_vdim[0]); 

  alpha_vdim[16] = -3.464101615137754*phi[2]*dv11*dx11; 
  alpha_vdim[17] = -3.464101615137754*phi[3]*dv11*dx11; 
  alpha_mid += fabs(0.125*alpha_vdim[16]); 

  out[1] += 0.4330127018922193*(alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(f[4]*alpha_cdim[20]+f[0]*alpha_cdim[16]); 
  out[3] += 0.4330127018922193*(alpha_vdim[2]*f[2]+alpha_vdim[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[1]*alpha_vdim[17]+f[0]*alpha_vdim[16]); 
  out[5] += 0.4330127018922193*(f[8]*alpha_cdim[20]+f[1]*alpha_cdim[16]+alpha_cdim[3]*f[7]+alpha_cdim[0]*f[2]); 
  out[6] += 0.4330127018922193*(alpha_vdim[2]*f[5]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_vdim[0]*f[1]); 
  out[7] += 0.4330127018922193*(f[10]*alpha_cdim[20]+f[3]*alpha_cdim[16]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[8] += 0.4330127018922193*(f[0]*alpha_vdim[17]+f[1]*alpha_vdim[16]+alpha_cdim[3]*f[10]+alpha_cdim[0]*f[4]); 
  out[9] += 0.4330127018922193*(f[0]*alpha_cdim[20]+f[5]*alpha_vdim[17]+f[2]*alpha_vdim[16]+f[4]*alpha_cdim[16]); 
  out[10] += 0.4330127018922193*(f[6]*alpha_vdim[17]+f[3]*alpha_vdim[16]+alpha_vdim[2]*f[9]+alpha_vdim[0]*f[4]); 
  out[11] += 0.4330127018922193*(f[13]*alpha_cdim[20]+f[6]*alpha_cdim[16]+alpha_cdim[0]*f[7]+alpha_vdim[0]*f[5]+f[2]*alpha_cdim[3]+f[1]*alpha_vdim[2]); 
  out[12] += 0.4330127018922193*(f[1]*alpha_cdim[20]+f[2]*alpha_vdim[17]+f[5]*alpha_vdim[16]+f[8]*alpha_cdim[16]+alpha_cdim[3]*f[14]+alpha_cdim[0]*f[9]); 
  out[13] += 0.4330127018922193*(f[3]*alpha_vdim[17]+f[6]*alpha_vdim[16]+alpha_vdim[2]*f[12]+alpha_cdim[0]*f[10]+alpha_vdim[0]*f[8]+alpha_cdim[3]*f[4]); 
  out[14] += 0.4330127018922193*(f[3]*alpha_cdim[20]+f[11]*alpha_vdim[17]+f[7]*alpha_vdim[16]+f[10]*alpha_cdim[16]+alpha_vdim[0]*f[9]+alpha_vdim[2]*f[4]); 
  out[15] += 0.4330127018922193*(f[6]*alpha_cdim[20]+f[7]*alpha_vdim[17]+f[11]*alpha_vdim[16]+f[13]*alpha_cdim[16]+alpha_cdim[0]*f[14]+alpha_vdim[0]*f[12]+alpha_cdim[3]*f[9]+alpha_vdim[2]*f[8]); 

  return alpha_mid; 
} 

