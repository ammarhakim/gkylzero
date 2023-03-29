#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_sr_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // qmem:      q/m*EM fields.
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx10 = 2/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &qmem[0]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &qmem[2]; 
  const double *p1_over_gamma = &p_over_gamma[4]; 

  const double *B2 = &qmem[10]; 
  double cflFreq_mid = 0.0; 
  double alpha_cdim[8] = {0.0}; 
  double alpha_vdim[16] = {0.0}; 

  alpha_cdim[0] = 1.414213562373095*p0_over_gamma[0]*dx10; 
  alpha_cdim[2] = 1.414213562373095*p0_over_gamma[1]*dx10; 
  alpha_cdim[3] = 1.414213562373095*p0_over_gamma[2]*dx10; 
  alpha_cdim[6] = 1.414213562373095*p0_over_gamma[3]*dx10; 
  cflFreq_mid += 3.0*fabs(0.1767766952966368*alpha_cdim[0]); 

  alpha_vdim[0] = (B2[0]*p1_over_gamma[0]+2.0*E0[0])*dv10; 
  alpha_vdim[1] = (2.0*E0[1]+p1_over_gamma[0]*B2[1])*dv10; 
  alpha_vdim[2] = B2[0]*p1_over_gamma[1]*dv10; 
  alpha_vdim[3] = B2[0]*p1_over_gamma[2]*dv10; 
  alpha_vdim[4] = B2[1]*p1_over_gamma[1]*dv10; 
  alpha_vdim[5] = B2[1]*p1_over_gamma[2]*dv10; 
  alpha_vdim[6] = B2[0]*p1_over_gamma[3]*dv10; 
  alpha_vdim[7] = B2[1]*p1_over_gamma[3]*dv10; 
  cflFreq_mid += 3.0*fabs(0.1767766952966368*alpha_vdim[0]); 

  alpha_vdim[8] = (2.0*E1[0]-1.0*B2[0]*p0_over_gamma[0])*dv11; 
  alpha_vdim[9] = (2.0*E1[1]-1.0*p0_over_gamma[0]*B2[1])*dv11; 
  alpha_vdim[10] = -1.0*B2[0]*p0_over_gamma[1]*dv11; 
  alpha_vdim[11] = -1.0*B2[0]*p0_over_gamma[2]*dv11; 
  alpha_vdim[12] = -1.0*B2[1]*p0_over_gamma[1]*dv11; 
  alpha_vdim[13] = -1.0*B2[1]*p0_over_gamma[2]*dv11; 
  alpha_vdim[14] = -1.0*B2[0]*p0_over_gamma[3]*dv11; 
  alpha_vdim[15] = -1.0*B2[1]*p0_over_gamma[3]*dv11; 
  cflFreq_mid += 3.0*fabs(0.1767766952966368*alpha_vdim[8]); 

  out[1] += 0.6123724356957944*(alpha_cdim[6]*f[6]+alpha_cdim[3]*f[3]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[7]*alpha_vdim[15]+f[6]*alpha_vdim[14]+f[5]*alpha_vdim[13]+f[4]*alpha_vdim[12]+f[3]*alpha_vdim[11]+f[2]*alpha_vdim[10]+f[1]*alpha_vdim[9]+f[0]*alpha_vdim[8]); 
  out[4] += 0.6123724356957944*(alpha_vdim[6]*f[7]+f[6]*(alpha_vdim[7]+alpha_cdim[3])+f[3]*alpha_cdim[6]+alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_vdim[2]*f[4]+f[2]*(alpha_vdim[4]+alpha_cdim[0])+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.6123724356957944*(f[6]*alpha_vdim[15]+f[7]*alpha_vdim[14]+f[3]*alpha_vdim[13]+f[2]*alpha_vdim[12]+f[5]*alpha_vdim[11]+f[4]*alpha_vdim[10]+f[0]*alpha_vdim[9]+f[1]*alpha_vdim[8]+alpha_cdim[2]*f[6]+f[2]*alpha_cdim[6]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]); 
  out[6] += 0.6123724356957944*(f[5]*alpha_vdim[15]+f[3]*alpha_vdim[14]+f[7]*alpha_vdim[13]+f[1]*alpha_vdim[12]+f[6]*alpha_vdim[11]+f[0]*alpha_vdim[10]+f[4]*alpha_vdim[9]+f[2]*alpha_vdim[8]+alpha_vdim[4]*f[7]+f[4]*alpha_vdim[7]+alpha_vdim[2]*f[6]+f[2]*alpha_vdim[6]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 0.6123724356957944*(f[3]*alpha_vdim[15]+f[5]*alpha_vdim[14]+f[6]*alpha_vdim[13]+f[0]*alpha_vdim[12]+f[7]*alpha_vdim[11]+f[1]*alpha_vdim[10]+f[2]*alpha_vdim[9]+f[4]*alpha_vdim[8]+alpha_vdim[2]*f[7]+f[2]*alpha_vdim[7]+(alpha_vdim[4]+alpha_cdim[0])*f[6]+f[4]*alpha_vdim[6]+f[0]*alpha_cdim[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]+f[2]*alpha_cdim[3]); 

  return cflFreq_mid; 
} 
