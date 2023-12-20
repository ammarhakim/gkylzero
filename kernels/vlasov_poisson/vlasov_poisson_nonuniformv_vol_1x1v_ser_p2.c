#include <gkyl_vlasov_poisson_kernels.h> 

GKYL_CU_DH double vlasov_poisson_nonuniformv_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *vcoord, const double *field, const double *f, double* GKYL_RESTRICT out) 
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

  double cflFreq = 0.0; 
  double alpha_cdim[8]; 
  double alpha_vdim[8]; 

  alpha_cdim[0] = 2.828427124746191*vcoord[0]*rdx; 
  alpha_cdim[2] = 2.828427124746191*vcoord[1]*rdx; 
  alpha_cdim[5] = 2.828427124746191*vcoord[2]*rdx; 
  cflFreq += 5.0*rdx*fmax(fabs(1.58113883008419*vcoord[2]-1.224744871391589*vcoord[1]+0.7071067811865475*vcoord[0]),fabs(1.58113883008419*vcoord[2]+1.224744871391589*vcoord[1]+0.7071067811865475*vcoord[0])); 

  alpha_vdim[0] = -2.449489742783178*phi[1]*rdvx2*rdx2; 
  alpha_vdim[1] = -5.477225575051662*phi[2]*rdvx2*rdx2; 
  cflFreq += 5.0*fabs(0.25*alpha_vdim[0]); 

  out[1] += 0.8660254037844386*(alpha_cdim[5]*f[5]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.8660254037844386*(alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.7745966692414833*(alpha_cdim[2]*f[5]+f[2]*alpha_cdim[5]+alpha_vdim[1]*f[4])+0.8660254037844386*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[4] += 1.936491673103709*(alpha_cdim[5]*f[7]+alpha_cdim[2]*f[3]+alpha_cdim[0]*f[1]); 
  out[5] += 1.936491673103709*(alpha_vdim[1]*f[3]+alpha_vdim[0]*f[2]); 
  out[6] += 1.732050807568877*(alpha_cdim[2]*f[7]+f[3]*alpha_cdim[5])+0.8660254037844386*alpha_vdim[0]*f[4]+1.936491673103709*alpha_cdim[0]*f[3]+f[1]*(1.936491673103709*alpha_cdim[2]+0.7745966692414833*alpha_vdim[1]); 
  out[7] += 1.732050807568877*alpha_vdim[1]*f[6]+0.5532833351724881*alpha_cdim[5]*f[5]+0.8660254037844386*(alpha_cdim[0]*f[5]+f[0]*alpha_cdim[5])+1.936491673103709*alpha_vdim[0]*f[3]+(0.7745966692414833*alpha_cdim[2]+1.936491673103709*alpha_vdim[1])*f[2]; 

  return cflFreq; 
} 


GKYL_CU_DH double vlasov_poisson_extem_nonuniformv_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *vcoord, const double *field, const double *f, double* GKYL_RESTRICT out) 
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

  double cflFreq = 0.0; 
  double alpha_cdim[8]; 
  double alpha_vdim[8]; 

  alpha_cdim[0] = 2.828427124746191*vcoord[0]*rdx; 
  alpha_cdim[2] = 2.828427124746191*vcoord[1]*rdx; 
  alpha_cdim[5] = 2.828427124746191*vcoord[2]*rdx; 
  cflFreq += 5.0*rdx*fmax(fabs(1.58113883008419*vcoord[2]-1.224744871391589*vcoord[1]+0.7071067811865475*vcoord[0]),fabs(1.58113883008419*vcoord[2]+1.224744871391589*vcoord[1]+0.7071067811865475*vcoord[0])); 

  alpha_vdim[0] = -2.449489742783178*phi[1]*rdvx2*rdx2; 
  alpha_vdim[1] = -5.477225575051662*phi[2]*rdvx2*rdx2; 
  cflFreq += 5.0*fabs(0.25*alpha_vdim[0]); 

  out[1] += 0.8660254037844386*(alpha_cdim[5]*f[5]+alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.8660254037844386*(alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.7745966692414833*(alpha_cdim[2]*f[5]+f[2]*alpha_cdim[5]+alpha_vdim[1]*f[4])+0.8660254037844386*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[4] += 1.936491673103709*(alpha_cdim[5]*f[7]+alpha_cdim[2]*f[3]+alpha_cdim[0]*f[1]); 
  out[5] += 1.936491673103709*(alpha_vdim[1]*f[3]+alpha_vdim[0]*f[2]); 
  out[6] += 1.732050807568877*(alpha_cdim[2]*f[7]+f[3]*alpha_cdim[5])+0.8660254037844386*alpha_vdim[0]*f[4]+1.936491673103709*alpha_cdim[0]*f[3]+f[1]*(1.936491673103709*alpha_cdim[2]+0.7745966692414833*alpha_vdim[1]); 
  out[7] += 1.732050807568877*alpha_vdim[1]*f[6]+0.5532833351724881*alpha_cdim[5]*f[5]+0.8660254037844386*(alpha_cdim[0]*f[5]+f[0]*alpha_cdim[5])+1.936491673103709*alpha_vdim[0]*f[3]+(0.7745966692414833*alpha_cdim[2]+1.936491673103709*alpha_vdim[1])*f[2]; 

  return cflFreq; 
} 

