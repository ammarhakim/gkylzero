#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_pkpm_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *u_i, const double *p_ij, const double *bvar, const double *rho_inv_b, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:       bulk flow velocity (ux, uy, uz).
  // p_ij:      pressure tensor (P_xx, P_xy, P_xz, P_yy, P_yz, P_zz).
  // bvar:      magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // rho_inv_b: b_i/rho (for pressure force 1/rho * b . div(P)).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx0 = 2.0/dxv[0]; 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[3]; 
  const double *uz = &u_i[6]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[3]; 
  const double *Pxz = &p_ij[6]; 
  const double *Pyy = &p_ij[9]; 
  const double *Pyz = &p_ij[12]; 
  const double *Pzz = &p_ij[15]; 

  const double *bx = &bvar[0]; 
  const double *by = &bvar[3]; 
  const double *bz = &bvar[6]; 
  const double *bxbx = &bvar[9]; 
  const double *bxby = &bvar[12]; 
  const double *bxbz = &bvar[15]; 
  const double *byby = &bvar[18]; 
  const double *bybz = &bvar[21]; 
  const double *bzbz = &bvar[24]; 

  const double *rho_inv_bx = &rho_inv_b[0]; 
  const double *rho_inv_by = &rho_inv_b[3]; 
  const double *rho_inv_bz = &rho_inv_b[6]; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[8] = {0.0}; 
  double alpha_vdim[8] = {0.0}; 

  alpha_cdim[0] = 1.414213562373095*dx0*(bx[0]*wvpar+ux[0]); 
  alpha_cdim[1] = 1.414213562373095*dx0*(bx[1]*wvpar+ux[1]); 
  alpha_cdim[2] = 0.408248290463863*bx[0]*dvpar*dx0; 
  alpha_cdim[3] = 0.408248290463863*bx[1]*dvpar*dx0; 
  alpha_cdim[4] = 1.414213562373095*dx0*(bx[2]*wvpar+ux[2]); 
  alpha_cdim[6] = 0.408248290463863*bx[2]*dvpar*dx0; 
  cflFreq_mid += 5.0*fabs(0.25*alpha_cdim[0]-0.2795084971874737*alpha_cdim[4]); 

  alpha_vdim[0] = ((-3.872983346207417*bxbz[1]*uz[2])-3.872983346207417*bxby[1]*uy[2]-3.872983346207417*bxbx[1]*ux[2]-1.732050807568877*bxbz[0]*uz[1]-1.732050807568877*bxby[0]*uy[1]-1.732050807568877*bxbx[0]*ux[1])*dv1par*dx0*wvpar+(3.872983346207417*rho_inv_bz[1]*Pxz[2]+3.872983346207417*rho_inv_by[1]*Pxy[2]+3.872983346207417*rho_inv_bx[1]*Pxx[2]+1.732050807568877*rho_inv_bz[0]*Pxz[1]+1.732050807568877*rho_inv_by[0]*Pxy[1]+1.732050807568877*rho_inv_bx[0]*Pxx[1])*dv1par*dx0; 
  alpha_vdim[1] = ((-3.464101615137754*bxbz[2]*uz[2])-3.872983346207417*bxbz[0]*uz[2]-3.464101615137754*bxby[2]*uy[2]-3.872983346207417*bxby[0]*uy[2]-3.464101615137754*bxbx[2]*ux[2]-3.872983346207417*bxbx[0]*ux[2]-1.732050807568877*bxbz[1]*uz[1]-1.732050807568877*bxby[1]*uy[1]-1.732050807568877*bxbx[1]*ux[1])*dv1par*dx0*wvpar+(3.464101615137754*Pxz[2]*rho_inv_bz[2]+3.464101615137754*Pxy[2]*rho_inv_by[2]+3.464101615137754*Pxx[2]*rho_inv_bx[2]+3.872983346207417*rho_inv_bz[0]*Pxz[2]+3.872983346207417*rho_inv_by[0]*Pxy[2]+3.872983346207417*rho_inv_bx[0]*Pxx[2]+1.732050807568877*Pxz[1]*rho_inv_bz[1]+1.732050807568877*Pxy[1]*rho_inv_by[1]+1.732050807568877*Pxx[1]*rho_inv_bx[1])*dv1par*dx0; 
  alpha_vdim[2] = ((-1.118033988749895*bxbz[1]*uz[2])-1.118033988749895*bxby[1]*uy[2]-1.118033988749895*bxbx[1]*ux[2]-0.5*bxbz[0]*uz[1]-0.5*bxby[0]*uy[1]-0.5*bxbx[0]*ux[1])*dv1par*dvpar*dx0; 
  alpha_vdim[3] = ((-1.0*bxbz[2]*uz[2])-1.118033988749895*bxbz[0]*uz[2]-1.0*bxby[2]*uy[2]-1.118033988749895*bxby[0]*uy[2]-1.0*bxbx[2]*ux[2]-1.118033988749895*bxbx[0]*ux[2]-0.5*bxbz[1]*uz[1]-0.5*bxby[1]*uy[1]-0.5*bxbx[1]*ux[1])*dv1par*dvpar*dx0; 
  alpha_vdim[4] = ((-3.464101615137754*bxbz[1]*uz[2])-3.464101615137754*bxby[1]*uy[2]-3.464101615137754*bxbx[1]*ux[2]-1.732050807568877*uz[1]*bxbz[2]-1.732050807568877*uy[1]*bxby[2]-1.732050807568877*ux[1]*bxbx[2])*dv1par*dx0*wvpar+(1.732050807568877*Pxz[1]*rho_inv_bz[2]+1.732050807568877*Pxy[1]*rho_inv_by[2]+1.732050807568877*Pxx[1]*rho_inv_bx[2]+3.464101615137754*rho_inv_bz[1]*Pxz[2]+3.464101615137754*rho_inv_by[1]*Pxy[2]+3.464101615137754*rho_inv_bx[1]*Pxx[2])*dv1par*dx0; 
  alpha_vdim[6] = ((-1.0*bxbz[1]*uz[2])-1.0*bxby[1]*uy[2]-1.0*bxbx[1]*ux[2]-0.5000000000000001*uz[1]*bxbz[2]-0.5000000000000001*uy[1]*bxby[2]-0.5000000000000001*ux[1]*bxbx[2])*dv1par*dvpar*dx0; 
  cflFreq_mid += 5.0*fabs(0.25*alpha_vdim[0]-0.2795084971874737*alpha_vdim[4]); 

  out[1] += 0.8660254037844386*(alpha_cdim[6]*f[6]+alpha_cdim[4]*f[4]+alpha_cdim[3]*f[3]+alpha_cdim[2]*f[2]+alpha_cdim[1]*f[1]+alpha_cdim[0]*f[0]); 
  out[2] += 0.8660254037844386*(alpha_vdim[6]*f[6]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.7745966692414833*alpha_cdim[3]*f[7]+0.8660254037844386*alpha_cdim[4]*f[6]+0.7745966692414833*(alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6])+0.8660254037844386*f[4]*alpha_cdim[6]+0.7745966692414833*(alpha_cdim[2]*f[5]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4])+0.8660254037844386*((alpha_vdim[2]+alpha_cdim[1])*f[3]+f[2]*alpha_vdim[3]+f[1]*alpha_cdim[3]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[4] += 1.732050807568877*(alpha_cdim[3]*f[6]+f[3]*alpha_cdim[6]+alpha_cdim[1]*f[4]+f[1]*alpha_cdim[4])+1.936491673103709*(alpha_cdim[2]*f[3]+f[2]*alpha_cdim[3]+alpha_cdim[0]*f[1]+f[0]*alpha_cdim[1]); 
  out[5] += 1.732050807568877*alpha_vdim[3]*f[7]+1.936491673103709*(alpha_vdim[4]*f[6]+f[4]*alpha_vdim[6])+1.732050807568877*alpha_vdim[2]*f[5]+1.936491673103709*(alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[6] += (1.549193338482967*alpha_cdim[6]+1.732050807568877*alpha_cdim[2])*f[7]+(0.5532833351724881*alpha_vdim[6]+0.8660254037844386*alpha_vdim[2]+1.732050807568877*alpha_cdim[1])*f[6]+0.8660254037844386*f[2]*alpha_vdim[6]+1.732050807568877*(f[1]*alpha_cdim[6]+alpha_cdim[3]*f[5])+(0.5532833351724881*alpha_vdim[4]+1.732050807568877*alpha_cdim[3])*f[4]+0.8660254037844386*(alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4])+f[3]*(1.732050807568877*alpha_cdim[4]+0.7745966692414833*alpha_vdim[3])+1.936491673103709*(alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_cdim[1]*f[2])+f[1]*(1.936491673103709*alpha_cdim[2]+0.7745966692414833*alpha_vdim[1]); 
  out[7] += (1.549193338482967*alpha_vdim[6]+1.732050807568877*alpha_vdim[2]+0.8660254037844386*alpha_cdim[1])*f[7]+0.7745966692414833*alpha_cdim[6]*f[6]+1.732050807568877*(alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6])+(1.732050807568877*alpha_vdim[3]+0.8660254037844386*alpha_cdim[0])*f[5]+1.732050807568877*alpha_vdim[3]*f[4]+f[3]*(1.732050807568877*alpha_vdim[4]+0.7745966692414833*alpha_cdim[3])+1.936491673103709*(alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3])+0.7745966692414833*alpha_cdim[2]*f[2]+1.936491673103709*(alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 

  return cflFreq_mid; 
} 
