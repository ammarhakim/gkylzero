#include <gkyl_vlasov_poisson_kernels.h> 

GKYL_CU_DH double vlasov_poisson_nonuniformv_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *vcoord, const double *field, const double *f, double* GKYL_RESTRICT out) 
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
  double alpha_cdim[20]; 
  double alpha_vdim[40]; 

  alpha_cdim[0] = 2.828427124746191*vcoord[0]*rdx; 
  alpha_cdim[2] = 2.828427124746191*vcoord[1]*rdx; 
  alpha_cdim[8] = 2.828427124746191*vcoord[4]*rdx; 
  cflFreq += 5.0*rdx*fmax(fabs(1.118033988749895*vcoord[4]-0.8660254037844386*vcoord[1]+0.5*vcoord[0]),fabs(1.118033988749895*vcoord[4]+0.8660254037844386*vcoord[1]+0.5*vcoord[0])); 

  alpha_vdim[0] = -3.464101615137754*phi[1]*rdvx2*rdx2; 
  alpha_vdim[1] = -7.745966692414834*phi[2]*rdvx2*rdx2; 
  cflFreq += 5.0*fabs(0.1767766952966368*alpha_vdim[0]); 

  cflFreq += 5.0*fabs(0.0); 

  out[1] += 0.6123724356957944*(alpha_cdim[8]*f[8]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.5477225575051661*(alpha_cdim[2]*f[8]+f[2]*alpha_cdim[8]+alpha_vdim[1]*f[7])+0.6123724356957944*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.6123724356957944*(alpha_cdim[8]*f[14]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*(alpha_vdim[1]*f[5]+alpha_vdim[0]*f[3]); 
  out[7] += 1.369306393762915*(alpha_cdim[8]*f[12]+alpha_cdim[2]*f[4]+alpha_cdim[0]*f[1]); 
  out[8] += 1.369306393762915*(alpha_vdim[1]*f[4]+alpha_vdim[0]*f[2]); 
  out[10] += 0.5477225575051661*(alpha_cdim[2]*f[14]+alpha_vdim[1]*f[13]+f[6]*alpha_cdim[8])+0.6123724356957944*(alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]); 
  out[11] += 1.224744871391589*(alpha_cdim[2]*f[12]+f[4]*alpha_cdim[8])+0.6123724356957944*alpha_vdim[0]*f[7]+1.369306393762915*alpha_cdim[0]*f[4]+f[1]*(1.369306393762915*alpha_cdim[2]+0.5477225575051661*alpha_vdim[1]); 
  out[12] += 1.224744871391589*alpha_vdim[1]*f[11]+0.3912303982179757*alpha_cdim[8]*f[8]+0.6123724356957944*(alpha_cdim[0]*f[8]+f[0]*alpha_cdim[8])+1.369306393762915*alpha_vdim[0]*f[4]+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[2]; 
  out[13] += 1.369306393762915*(alpha_cdim[8]*f[18]+alpha_cdim[2]*f[10]+alpha_cdim[0]*f[5]); 
  out[14] += 1.369306393762915*(alpha_vdim[1]*f[10]+alpha_vdim[0]*f[6]); 
  out[15] += 0.6123724356957944*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[9]); 
  out[16] += 0.6123724356957944*(alpha_vdim[1]*f[15]+alpha_vdim[0]*f[9]); 
  out[17] += 1.224744871391589*alpha_cdim[2]*f[18]+0.6123724356957944*alpha_vdim[0]*f[13]+(1.224744871391589*alpha_cdim[8]+1.369306393762915*alpha_cdim[0])*f[10]+(1.369306393762915*alpha_cdim[2]+0.5477225575051661*alpha_vdim[1])*f[5]; 
  out[18] += 1.224744871391589*alpha_vdim[1]*f[17]+(0.3912303982179757*alpha_cdim[8]+0.6123724356957944*alpha_cdim[0])*f[14]+1.369306393762915*alpha_vdim[0]*f[10]+0.6123724356957944*f[3]*alpha_cdim[8]+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[6]; 
  out[19] += 0.5477225575051661*alpha_cdim[8]*f[16]+0.6123724356957944*(alpha_cdim[0]*f[16]+alpha_vdim[0]*f[15]+(alpha_cdim[2]+alpha_vdim[1])*f[9]); 

  return cflFreq; 
} 


GKYL_CU_DH double vlasov_poisson_extem_nonuniformv_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *vcoord, const double *field, const double *f, double* GKYL_RESTRICT out) 
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

  const double *A0 = &field[3]; 
  const double *A1 = &field[6]; 
  double cflFreq = 0.0; 
  double alpha_cdim[20]; 
  double alpha_vdim[40]; 

  alpha_cdim[0] = 2.828427124746191*vcoord[0]*rdx; 
  alpha_cdim[2] = 2.828427124746191*vcoord[1]*rdx; 
  alpha_cdim[8] = 2.828427124746191*vcoord[4]*rdx; 
  cflFreq += 5.0*rdx*fmax(fabs(1.118033988749895*vcoord[4]-0.8660254037844386*vcoord[1]+0.5*vcoord[0]),fabs(1.118033988749895*vcoord[4]+0.8660254037844386*vcoord[1]+0.5*vcoord[0])); 

  alpha_vdim[0] = (1.732050807568877*A1[1]*vcoord[8]-3.464101615137754*phi[1])*rdvx2*rdx2; 
  alpha_vdim[1] = (3.872983346207417*A1[2]*vcoord[8]-7.745966692414834*phi[2])*rdvx2*rdx2; 
  alpha_vdim[3] = 1.732050807568877*A1[1]*vcoord[10]*rdvx2*rdx2; 
  alpha_vdim[5] = 3.872983346207417*A1[2]*vcoord[10]*rdvx2*rdx2; 
  alpha_vdim[9] = 1.732050807568877*A1[1]*vcoord[13]*rdvx2*rdx2; 
  alpha_vdim[15] = 3.872983346207417*A1[2]*vcoord[13]*rdvx2*rdx2; 
  cflFreq += 5.0*fabs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*alpha_vdim[9]); 

  alpha_vdim[20] = -1.732050807568877*vcoord[0]*A1[1]*rdvy2*rdx2; 
  alpha_vdim[21] = -3.872983346207417*vcoord[0]*A1[2]*rdvy2*rdx2; 
  alpha_vdim[22] = -1.732050807568877*A1[1]*vcoord[1]*rdvy2*rdx2; 
  alpha_vdim[24] = -3.872983346207417*vcoord[1]*A1[2]*rdvy2*rdx2; 
  alpha_vdim[28] = -1.732050807568877*A1[1]*vcoord[4]*rdvy2*rdx2; 
  alpha_vdim[32] = -3.872983346207417*A1[2]*vcoord[4]*rdvy2*rdx2; 
  cflFreq += 5.0*fabs(0.1767766952966368*alpha_vdim[20]-0.1976423537605236*alpha_vdim[28]); 

  out[1] += 0.6123724356957944*(alpha_cdim[8]*f[8]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[15]*f[15]+alpha_vdim[9]*f[9]+alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[12]*alpha_vdim[32]+f[8]*alpha_vdim[28]+f[4]*alpha_vdim[24]+f[2]*alpha_vdim[22]+f[1]*alpha_vdim[21]+f[0]*alpha_vdim[20]); 
  out[4] += 0.6123724356957944*(alpha_vdim[9]*f[15]+f[9]*alpha_vdim[15])+0.5477225575051661*(alpha_vdim[5]*f[13]+alpha_cdim[2]*f[8]+f[2]*alpha_cdim[8]+alpha_vdim[1]*f[7])+0.6123724356957944*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.6123724356957944*(f[8]*alpha_vdim[32]+f[12]*alpha_vdim[28])+0.5477225575051661*f[11]*alpha_vdim[24]+0.6123724356957944*(f[2]*alpha_vdim[24]+f[4]*alpha_vdim[22])+0.5477225575051661*f[7]*alpha_vdim[21]+0.6123724356957944*(f[0]*alpha_vdim[21]+f[1]*alpha_vdim[20]+alpha_cdim[8]*f[14]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.5477225575051661*(f[4]*alpha_vdim[32]+f[2]*alpha_vdim[28])+(0.5477225575051661*f[12]+0.6123724356957944*f[1])*alpha_vdim[24]+0.5477225575051661*f[8]*alpha_vdim[22]+0.6123724356957944*(f[0]*alpha_vdim[22]+f[4]*alpha_vdim[21]+f[2]*alpha_vdim[20])+0.5477225575051661*(alpha_vdim[5]*f[15]+f[5]*alpha_vdim[15]+alpha_vdim[3]*f[9]+f[3]*alpha_vdim[9])+0.6123724356957944*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 1.369306393762915*(alpha_cdim[8]*f[12]+alpha_cdim[2]*f[4]+alpha_cdim[0]*f[1]); 
  out[8] += 1.369306393762915*(alpha_vdim[15]*f[19]+alpha_vdim[9]*f[16]+alpha_vdim[5]*f[10]+alpha_vdim[3]*f[6]+alpha_vdim[1]*f[4]+alpha_vdim[0]*f[2]); 
  out[9] += 1.369306393762915*(f[18]*alpha_vdim[32]+f[14]*alpha_vdim[28]+f[10]*alpha_vdim[24]+f[6]*alpha_vdim[22]+f[5]*alpha_vdim[21]+f[3]*alpha_vdim[20]); 
  out[10] += 0.4898979485566357*f[11]*alpha_vdim[32]+0.5477225575051661*(f[2]*alpha_vdim[32]+f[4]*alpha_vdim[28])+(0.5477225575051661*(f[8]+f[7])+0.6123724356957944*f[0])*alpha_vdim[24]+(0.5477225575051661*f[12]+0.6123724356957944*f[1])*alpha_vdim[22]+0.5477225575051661*f[11]*alpha_vdim[21]+0.6123724356957944*(f[2]*alpha_vdim[21]+f[4]*alpha_vdim[20])+0.5477225575051661*alpha_vdim[3]*f[15]+0.4898979485566357*f[13]*alpha_vdim[15]+0.5477225575051661*(f[3]*alpha_vdim[15]+alpha_cdim[2]*f[14]+alpha_vdim[1]*f[13]+alpha_vdim[5]*f[9]+f[5]*alpha_vdim[9]+f[6]*alpha_cdim[8]+alpha_vdim[5]*f[7])+0.6123724356957944*(alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]); 
  out[11] += 0.5477225575051661*alpha_vdim[15]*f[15]+0.6123724356957944*alpha_vdim[3]*f[13]+1.224744871391589*(alpha_cdim[2]*f[12]+f[4]*alpha_cdim[8])+0.6123724356957944*alpha_vdim[0]*f[7]+0.5477225575051661*alpha_vdim[5]*f[5]+1.369306393762915*alpha_cdim[0]*f[4]+f[1]*(1.369306393762915*alpha_cdim[2]+0.5477225575051661*alpha_vdim[1]); 
  out[12] += 1.369306393762915*alpha_vdim[9]*f[19]+1.224744871391589*alpha_vdim[5]*f[17]+1.369306393762915*alpha_vdim[15]*f[16]+1.224744871391589*alpha_vdim[1]*f[11]+1.369306393762915*alpha_vdim[3]*f[10]+0.3912303982179757*alpha_cdim[8]*f[8]+0.6123724356957944*(alpha_cdim[0]*f[8]+f[0]*alpha_cdim[8])+1.369306393762915*(alpha_vdim[5]*f[6]+alpha_vdim[0]*f[4])+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[2]; 
  out[13] += 0.5477225575051661*(f[12]*alpha_vdim[32]+f[4]*alpha_vdim[24])+0.6123724356957944*f[11]*alpha_vdim[22]+0.5477225575051661*f[1]*alpha_vdim[21]+0.6123724356957944*f[7]*alpha_vdim[20]+1.369306393762915*(alpha_cdim[8]*f[18]+alpha_cdim[2]*f[10]+alpha_cdim[0]*f[5]); 
  out[14] += (0.3912303982179757*f[12]+0.6123724356957944*f[1])*alpha_vdim[32]+(0.3912303982179757*f[8]+0.6123724356957944*f[0])*alpha_vdim[28]+0.5477225575051661*(f[4]*alpha_vdim[24]+f[2]*alpha_vdim[22])+0.6123724356957944*(f[12]*alpha_vdim[21]+f[8]*alpha_vdim[20])+1.224744871391589*(alpha_vdim[5]*f[19]+alpha_vdim[3]*f[16])+f[10]*(1.224744871391589*alpha_vdim[15]+1.369306393762915*alpha_vdim[1])+1.224744871391589*f[6]*alpha_vdim[9]+1.369306393762915*(alpha_vdim[0]*f[6]+f[4]*alpha_vdim[5]+f[2]*alpha_vdim[3]); 
  out[15] += 1.369306393762915*(f[14]*alpha_vdim[32]+f[18]*alpha_vdim[28])+1.224744871391589*f[17]*alpha_vdim[24]+1.369306393762915*(f[6]*alpha_vdim[24]+f[10]*alpha_vdim[22])+1.224744871391589*f[13]*alpha_vdim[21]+1.369306393762915*(f[3]*alpha_vdim[21]+f[5]*alpha_vdim[20])+0.6123724356957944*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[9]); 
  out[16] += 1.224744871391589*(f[10]*alpha_vdim[32]+f[6]*alpha_vdim[28])+(1.224744871391589*f[18]+1.369306393762915*f[5])*alpha_vdim[24]+1.224744871391589*f[14]*alpha_vdim[22]+1.369306393762915*(f[3]*alpha_vdim[22]+f[10]*alpha_vdim[21]+f[6]*alpha_vdim[20])+0.3912303982179757*alpha_vdim[15]*f[15]+0.6123724356957944*(alpha_vdim[1]*f[15]+f[1]*alpha_vdim[15])+0.3912303982179757*alpha_vdim[9]*f[9]+0.6123724356957944*(alpha_vdim[0]*f[9]+f[0]*alpha_vdim[9])+0.5477225575051661*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]); 
  out[17] += 0.4898979485566357*f[4]*alpha_vdim[32]+0.5477225575051661*f[11]*alpha_vdim[28]+(0.4898979485566357*f[12]+0.5477225575051661*f[1])*alpha_vdim[24]+0.6123724356957944*f[7]*alpha_vdim[22]+0.5477225575051661*f[4]*alpha_vdim[21]+0.6123724356957944*f[11]*alpha_vdim[20]+1.224744871391589*alpha_cdim[2]*f[18]+0.4898979485566357*(alpha_vdim[5]*f[15]+f[5]*alpha_vdim[15])+(0.5477225575051661*alpha_vdim[9]+0.6123724356957944*alpha_vdim[0])*f[13]+(1.224744871391589*alpha_cdim[8]+1.369306393762915*alpha_cdim[0])*f[10]+0.6123724356957944*alpha_vdim[3]*f[7]+1.369306393762915*alpha_cdim[2]*f[5]+0.5477225575051661*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]); 
  out[18] += (0.3912303982179757*f[8]+0.5477225575051661*f[7]+0.6123724356957944*f[0])*alpha_vdim[32]+(0.3912303982179757*f[12]+0.6123724356957944*f[1])*alpha_vdim[28]+0.4898979485566357*f[11]*alpha_vdim[24]+0.5477225575051661*(f[2]*alpha_vdim[24]+f[4]*alpha_vdim[22])+0.6123724356957944*(f[8]*alpha_vdim[21]+f[12]*alpha_vdim[20])+1.224744871391589*alpha_vdim[3]*f[19]+1.095445115010332*alpha_vdim[15]*f[17]+1.224744871391589*(alpha_vdim[1]*f[17]+alpha_vdim[5]*f[16]+f[6]*alpha_vdim[15])+(0.3912303982179757*alpha_cdim[8]+0.6123724356957944*alpha_cdim[0])*f[14]+1.224744871391589*alpha_vdim[5]*f[11]+(1.224744871391589*alpha_vdim[9]+1.369306393762915*alpha_vdim[0])*f[10]+0.6123724356957944*f[3]*alpha_cdim[8]+0.5477225575051661*alpha_cdim[2]*f[6]+1.369306393762915*(alpha_vdim[1]*f[6]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]); 
  out[19] += 1.095445115010332*f[17]*alpha_vdim[32]+1.224744871391589*(f[6]*alpha_vdim[32]+f[10]*alpha_vdim[28])+(1.224744871391589*(f[14]+f[13])+1.369306393762915*f[3])*alpha_vdim[24]+(1.224744871391589*f[18]+1.369306393762915*f[5])*alpha_vdim[22]+1.224744871391589*f[17]*alpha_vdim[21]+1.369306393762915*(f[6]*alpha_vdim[21]+f[10]*alpha_vdim[20])+(0.5477225575051661*alpha_cdim[8]+0.6123724356957944*alpha_cdim[0])*f[16]+(0.3912303982179757*alpha_vdim[9]+0.6123724356957944*alpha_vdim[0])*f[15]+(0.3912303982179757*f[9]+0.5477225575051661*f[7]+0.6123724356957944*f[0])*alpha_vdim[15]+0.4898979485566357*alpha_vdim[5]*f[13]+0.6123724356957944*((alpha_cdim[2]+alpha_vdim[1])*f[9]+f[1]*alpha_vdim[9])+0.5477225575051661*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]); 

  return cflFreq; 
} 

