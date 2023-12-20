#include <gkyl_vlasov_poisson_kernels.h> 

GKYL_CU_DH double vlasov_poisson_nonuniformv_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vcoord, const double *field, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // vcoord     Discrete (DG) velocity coordinate.
  // field:     potential (scaled by appropriate factors).
  // f:         Input distribution function.
  // out:       Incremented output.
  double rdx = 1./dxv[0];

  const double rdx2 = 2.*rdx; 
  const double *phi = &field[0]; 
  const double rdvx2 = 2./dxv[1]; 
  const double rdvy2 = 2./dxv[2]; 

  double cflFreq = 0.0; 
  double alpha_cdim[16]; 
  double alpha_vdim[32]; 

  alpha_cdim[0] = 2.828427124746191*vcoord[0]*rdx; 
  alpha_cdim[2] = 2.828427124746191*vcoord[1]*rdx; 
  alpha_cdim[8] = 2.828427124746191*vcoord[4]*rdx; 
  cflFreq += 3.0*rdx*fmax(fabs(1.118033988749895*vcoord[4]-0.8660254037844386*vcoord[1]+0.5*vcoord[0]),fabs(1.118033988749895*vcoord[4]+0.8660254037844386*vcoord[1]+0.5*vcoord[0])); 

  alpha_vdim[0] = -3.464101615137754*phi[1]*rdvx2*rdx2; 
  cflFreq += 5.0*fabs(0.1767766952966368*alpha_vdim[0]); 

  cflFreq += 5.0*fabs(0.0); 

  out[1] += 0.6123724356957944*(alpha_cdim[8]*f[8]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*alpha_vdim[0]*f[0]; 
  out[4] += 0.5477225575051661*(alpha_cdim[2]*f[8]+f[2]*alpha_cdim[8])+0.6123724356957944*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]); 
  out[5] += 0.6123724356957944*(alpha_cdim[8]*f[10]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*alpha_vdim[0]*f[3]; 
  out[7] += 0.5477225575051661*(alpha_cdim[2]*f[10]+f[6]*alpha_cdim[8])+0.6123724356957944*(alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+alpha_cdim[2]*f[3]); 
  out[8] += 1.369306393762915*alpha_vdim[0]*f[2]; 
  out[9] += 0.3912303982179757*alpha_cdim[8]*f[8]+0.6123724356957944*(alpha_cdim[0]*f[8]+f[0]*alpha_cdim[8])+1.369306393762915*alpha_vdim[0]*f[4]+0.5477225575051661*alpha_cdim[2]*f[2]; 
  out[10] += 1.369306393762915*alpha_vdim[0]*f[6]; 
  out[11] += 0.3912303982179757*alpha_cdim[8]*f[10]+0.6123724356957944*(alpha_cdim[0]*f[10]+f[3]*alpha_cdim[8])+1.369306393762915*alpha_vdim[0]*f[7]+0.5477225575051661*alpha_cdim[2]*f[6]; 
  out[13] += 0.6123724356957944*(alpha_cdim[2]*f[14]+alpha_cdim[0]*f[12]); 
  out[14] += 0.6123724356957944*alpha_vdim[0]*f[12]; 
  out[15] += 0.5477225575051661*alpha_cdim[8]*f[14]+0.6123724356957944*(alpha_cdim[0]*f[14]+alpha_vdim[0]*f[13]+alpha_cdim[2]*f[12]); 

  return cflFreq; 
} 


GKYL_CU_DH double vlasov_poisson_extem_nonuniformv_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vcoord, const double *field, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // vcoord     Discrete (DG) velocity coordinate.
  // field:     potentials, including external (scaled by appropriate factors).
  // f:         Input distribution function.
  // out:       Incremented output.
  double rdx = 1.0/dxv[0];

  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double rdx2 = 2.*rdx; 
  const double *phi = &field[0]; 
  const double rdvx2 = 2./dxv[1]; 
  const double rdvy2 = 2./dxv[2]; 

  const double *A0 = &field[2]; 
  const double *A1 = &field[4]; 
  double cflFreq = 0.0; 
  double alpha_cdim[16]; 
  double alpha_vdim[32]; 

  alpha_cdim[0] = 2.828427124746191*vcoord[0]*rdx; 
  alpha_cdim[2] = 2.828427124746191*vcoord[1]*rdx; 
  alpha_cdim[8] = 2.828427124746191*vcoord[4]*rdx; 
  cflFreq += 3.0*rdx*fmax(fabs(1.118033988749895*vcoord[4]-0.8660254037844386*vcoord[1]+0.5*vcoord[0]),fabs(1.118033988749895*vcoord[4]+0.8660254037844386*vcoord[1]+0.5*vcoord[0])); 

  alpha_vdim[0] = (1.732050807568877*A1[1]*vcoord[8]-3.464101615137754*phi[1])*rdvx2*rdx2; 
  alpha_vdim[3] = 1.732050807568877*A1[1]*vcoord[10]*rdvx2*rdx2; 
  alpha_vdim[12] = 1.732050807568877*A1[1]*vcoord[13]*rdvx2*rdx2; 
  cflFreq += 5.0*fabs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*alpha_vdim[12]); 

  alpha_vdim[16] = -1.732050807568877*vcoord[0]*A1[1]*rdvy2*rdx2; 
  alpha_vdim[18] = -1.732050807568877*A1[1]*vcoord[1]*rdvy2*rdx2; 
  alpha_vdim[24] = -1.732050807568877*A1[1]*vcoord[4]*rdvy2*rdx2; 
  cflFreq += 5.0*fabs(0.1767766952966368*alpha_vdim[16]-0.1976423537605236*alpha_vdim[24]); 

  out[1] += 0.6123724356957944*(alpha_cdim[8]*f[8]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[12]*f[12]+alpha_vdim[3]*f[3]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[8]*alpha_vdim[24]+f[2]*alpha_vdim[18]+f[0]*alpha_vdim[16]); 
  out[4] += 0.6123724356957944*alpha_vdim[12]*f[13]+0.5477225575051661*(alpha_cdim[2]*f[8]+f[2]*alpha_cdim[8])+0.6123724356957944*(alpha_vdim[3]*f[5]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]); 
  out[5] += 0.6123724356957944*(f[9]*alpha_vdim[24]+f[4]*alpha_vdim[18]+f[1]*alpha_vdim[16]+alpha_cdim[8]*f[10]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.5477225575051661*(f[2]*alpha_vdim[24]+f[8]*alpha_vdim[18])+0.6123724356957944*(f[0]*alpha_vdim[18]+f[2]*alpha_vdim[16])+0.5477225575051661*(alpha_vdim[3]*f[12]+f[3]*alpha_vdim[12])+0.6123724356957944*(alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 0.5477225575051661*(f[4]*alpha_vdim[24]+f[9]*alpha_vdim[18])+0.6123724356957944*(f[1]*alpha_vdim[18]+f[4]*alpha_vdim[16])+0.5477225575051661*(alpha_vdim[3]*f[13]+f[5]*alpha_vdim[12]+alpha_cdim[2]*f[10]+f[6]*alpha_cdim[8])+0.6123724356957944*(alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+alpha_cdim[2]*f[3]+f[1]*alpha_vdim[3]); 
  out[8] += 1.369306393762915*(alpha_vdim[12]*f[14]+alpha_vdim[3]*f[6]+alpha_vdim[0]*f[2]); 
  out[9] += 1.369306393762915*alpha_vdim[12]*f[15]+0.3912303982179757*alpha_cdim[8]*f[8]+0.6123724356957944*(alpha_cdim[0]*f[8]+f[0]*alpha_cdim[8])+1.369306393762915*(alpha_vdim[3]*f[7]+alpha_vdim[0]*f[4])+0.5477225575051661*alpha_cdim[2]*f[2]; 
  out[10] += (0.3912303982179757*f[8]+0.6123724356957944*f[0])*alpha_vdim[24]+0.5477225575051661*f[2]*alpha_vdim[18]+0.6123724356957944*f[8]*alpha_vdim[16]+1.224744871391589*(alpha_vdim[3]*f[14]+f[6]*alpha_vdim[12])+1.369306393762915*(alpha_vdim[0]*f[6]+f[2]*alpha_vdim[3]); 
  out[11] += (0.3912303982179757*f[9]+0.6123724356957944*f[1])*alpha_vdim[24]+0.5477225575051661*f[4]*alpha_vdim[18]+0.6123724356957944*f[9]*alpha_vdim[16]+1.224744871391589*(alpha_vdim[3]*f[15]+f[7]*alpha_vdim[12])+0.3912303982179757*alpha_cdim[8]*f[10]+0.6123724356957944*(alpha_cdim[0]*f[10]+f[3]*alpha_cdim[8])+1.369306393762915*alpha_vdim[0]*f[7]+0.5477225575051661*alpha_cdim[2]*f[6]+1.369306393762915*alpha_vdim[3]*f[4]; 
  out[12] += 1.369306393762915*(f[10]*alpha_vdim[24]+f[6]*alpha_vdim[18]+f[3]*alpha_vdim[16]); 
  out[13] += 1.369306393762915*(f[11]*alpha_vdim[24]+f[7]*alpha_vdim[18]+f[5]*alpha_vdim[16])+0.6123724356957944*(alpha_cdim[2]*f[14]+alpha_cdim[0]*f[12]); 
  out[14] += 1.224744871391589*(f[6]*alpha_vdim[24]+f[10]*alpha_vdim[18])+1.369306393762915*(f[3]*alpha_vdim[18]+f[6]*alpha_vdim[16])+0.3912303982179757*alpha_vdim[12]*f[12]+0.6123724356957944*(alpha_vdim[0]*f[12]+f[0]*alpha_vdim[12])+0.5477225575051661*alpha_vdim[3]*f[3]; 
  out[15] += 1.224744871391589*(f[7]*alpha_vdim[24]+f[11]*alpha_vdim[18])+1.369306393762915*(f[5]*alpha_vdim[18]+f[7]*alpha_vdim[16])+(0.5477225575051661*alpha_cdim[8]+0.6123724356957944*alpha_cdim[0])*f[14]+0.3912303982179757*alpha_vdim[12]*f[13]+0.6123724356957944*(alpha_vdim[0]*f[13]+alpha_cdim[2]*f[12]+f[1]*alpha_vdim[12])+0.5477225575051661*alpha_vdim[3]*f[5]; 

  return cflFreq; 
} 

