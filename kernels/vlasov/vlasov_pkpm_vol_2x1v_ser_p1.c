#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_pkpm_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *u_i, const double *p_ij, const double *bvar, const double *rho_inv_b, const double *f, double* GKYL_RESTRICT out) 
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
  const double dx1 = 2.0/dxv[1]; 
  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[4]; 
  const double *uz = &u_i[8]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[4]; 
  const double *Pxz = &p_ij[8]; 
  const double *Pyy = &p_ij[12]; 
  const double *Pyz = &p_ij[16]; 
  const double *Pzz = &p_ij[20]; 

  const double *bx = &bvar[0]; 
  const double *by = &bvar[4]; 
  const double *bz = &bvar[8]; 
  const double *bxbx = &bvar[12]; 
  const double *bxby = &bvar[16]; 
  const double *bxbz = &bvar[20]; 
  const double *byby = &bvar[24]; 
  const double *bybz = &bvar[28]; 
  const double *bzbz = &bvar[32]; 

  const double *rho_inv_bx = &rho_inv_b[0]; 
  const double *rho_inv_by = &rho_inv_b[4]; 
  const double *rho_inv_bz = &rho_inv_b[8]; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[24] = {0.0}; 
  double alpha_vdim[12] = {0.0}; 

  alpha_cdim[0] = 1.414213562373095*dx0*(bx[0]*wvpar+ux[0]); 
  alpha_cdim[1] = 1.414213562373095*dx0*(bx[1]*wvpar+ux[1]); 
  alpha_cdim[2] = 1.414213562373095*dx0*(bx[2]*wvpar+ux[2]); 
  alpha_cdim[3] = 0.408248290463863*bx[0]*dvpar*dx0; 
  alpha_cdim[4] = 1.414213562373095*dx0*(bx[3]*wvpar+ux[3]); 
  alpha_cdim[5] = 0.408248290463863*bx[1]*dvpar*dx0; 
  alpha_cdim[6] = 0.408248290463863*bx[2]*dvpar*dx0; 
  alpha_cdim[7] = 0.408248290463863*bx[3]*dvpar*dx0; 
  cflFreq_mid += 3.0*fabs(0.1767766952966368*alpha_cdim[0]); 

  alpha_cdim[12] = 1.414213562373095*dx1*(by[0]*wvpar+uy[0]); 
  alpha_cdim[13] = 1.414213562373095*dx1*(by[1]*wvpar+uy[1]); 
  alpha_cdim[14] = 1.414213562373095*dx1*(by[2]*wvpar+uy[2]); 
  alpha_cdim[15] = 0.408248290463863*by[0]*dvpar*dx1; 
  alpha_cdim[16] = 1.414213562373095*dx1*(by[3]*wvpar+uy[3]); 
  alpha_cdim[17] = 0.408248290463863*by[1]*dvpar*dx1; 
  alpha_cdim[18] = 0.408248290463863*by[2]*dvpar*dx1; 
  alpha_cdim[19] = 0.408248290463863*by[3]*dvpar*dx1; 
  cflFreq_mid += 3.0*fabs(0.1767766952966368*alpha_cdim[12]); 

  alpha_vdim[0] = ((-1.224744871391589*bybz[1]*uz[3])-1.224744871391589*byby[1]*uy[3]-1.224744871391589*bxby[1]*ux[3]-1.224744871391589*bybz[0]*uz[2]-1.224744871391589*byby[0]*uy[2]-1.224744871391589*bxby[0]*ux[2])*dv1par*dx1*wvpar+((-1.224744871391589*bxbz[2]*uz[3])-1.224744871391589*bxby[2]*uy[3]-1.224744871391589*bxbx[2]*ux[3]-1.224744871391589*bxbz[0]*uz[1]-1.224744871391589*bxby[0]*uy[1]-1.224744871391589*bxbx[0]*ux[1])*dv1par*dx0*wvpar+(1.224744871391589*rho_inv_bz[1]*Pyz[3]+1.224744871391589*rho_inv_by[1]*Pyy[3]+1.224744871391589*rho_inv_bx[1]*Pxy[3]+1.224744871391589*rho_inv_bz[0]*Pyz[2]+1.224744871391589*rho_inv_by[0]*Pyy[2]+1.224744871391589*rho_inv_bx[0]*Pxy[2])*dv1par*dx1+(1.224744871391589*rho_inv_bz[2]*Pxz[3]+1.224744871391589*rho_inv_by[2]*Pxy[3]+1.224744871391589*rho_inv_bx[2]*Pxx[3]+1.224744871391589*rho_inv_bz[0]*Pxz[1]+1.224744871391589*rho_inv_by[0]*Pxy[1]+1.224744871391589*rho_inv_bx[0]*Pxx[1])*dv1par*dx0; 
  alpha_vdim[1] = ((-1.224744871391589*bybz[0]*uz[3])-1.224744871391589*byby[0]*uy[3]-1.224744871391589*bxby[0]*ux[3]-1.224744871391589*bybz[1]*uz[2]-1.224744871391589*byby[1]*uy[2]-1.224744871391589*bxby[1]*ux[2])*dv1par*dx1*wvpar+((-1.224744871391589*bxbz[3]*uz[3])-1.224744871391589*bxby[3]*uy[3]-1.224744871391589*bxbx[3]*ux[3]-1.224744871391589*bxbz[1]*uz[1]-1.224744871391589*bxby[1]*uy[1]-1.224744871391589*bxbx[1]*ux[1])*dv1par*dx0*wvpar+(1.224744871391589*rho_inv_bz[0]*Pyz[3]+1.224744871391589*rho_inv_by[0]*Pyy[3]+1.224744871391589*rho_inv_bx[0]*Pxy[3]+1.224744871391589*rho_inv_bz[1]*Pyz[2]+1.224744871391589*rho_inv_by[1]*Pyy[2]+1.224744871391589*rho_inv_bx[1]*Pxy[2])*dv1par*dx1+(1.224744871391589*Pxz[3]*rho_inv_bz[3]+1.224744871391589*Pxy[3]*rho_inv_by[3]+1.224744871391589*Pxx[3]*rho_inv_bx[3]+1.224744871391589*Pxz[1]*rho_inv_bz[1]+1.224744871391589*Pxy[1]*rho_inv_by[1]+1.224744871391589*Pxx[1]*rho_inv_bx[1])*dv1par*dx0; 
  alpha_vdim[2] = ((-1.224744871391589*bybz[3]*uz[3])-1.224744871391589*byby[3]*uy[3]-1.224744871391589*bxby[3]*ux[3]-1.224744871391589*bybz[2]*uz[2]-1.224744871391589*byby[2]*uy[2]-1.224744871391589*bxby[2]*ux[2])*dv1par*dx1*wvpar+((-1.224744871391589*bxbz[0]*uz[3])-1.224744871391589*bxby[0]*uy[3]-1.224744871391589*bxbx[0]*ux[3]-1.224744871391589*uz[1]*bxbz[2]-1.224744871391589*uy[1]*bxby[2]-1.224744871391589*ux[1]*bxbx[2])*dv1par*dx0*wvpar+(1.224744871391589*Pyz[3]*rho_inv_bz[3]+1.224744871391589*Pyy[3]*rho_inv_by[3]+1.224744871391589*Pxy[3]*rho_inv_bx[3]+1.224744871391589*Pyz[2]*rho_inv_bz[2]+1.224744871391589*Pyy[2]*rho_inv_by[2]+1.224744871391589*Pxy[2]*rho_inv_bx[2])*dv1par*dx1+(1.224744871391589*rho_inv_bz[0]*Pxz[3]+1.224744871391589*rho_inv_by[0]*Pxy[3]+1.224744871391589*rho_inv_bx[0]*Pxx[3]+1.224744871391589*Pxz[1]*rho_inv_bz[2]+1.224744871391589*Pxy[1]*rho_inv_by[2]+1.224744871391589*Pxx[1]*rho_inv_bx[2])*dv1par*dx0; 
  alpha_vdim[3] = ((-0.3535533905932737*bybz[1]*uz[3])-0.3535533905932737*byby[1]*uy[3]-0.3535533905932737*bxby[1]*ux[3]-0.3535533905932737*bybz[0]*uz[2]-0.3535533905932737*byby[0]*uy[2]-0.3535533905932737*bxby[0]*ux[2])*dv1par*dvpar*dx1+((-0.3535533905932737*bxbz[2]*uz[3])-0.3535533905932737*bxby[2]*uy[3]-0.3535533905932737*bxbx[2]*ux[3]-0.3535533905932737*bxbz[0]*uz[1]-0.3535533905932737*bxby[0]*uy[1]-0.3535533905932737*bxbx[0]*ux[1])*dv1par*dvpar*dx0; 
  alpha_vdim[4] = ((-1.224744871391589*bybz[2]*uz[3])-1.224744871391589*byby[2]*uy[3]-1.224744871391589*bxby[2]*ux[3]-1.224744871391589*uz[2]*bybz[3]-1.224744871391589*uy[2]*byby[3]-1.224744871391589*ux[2]*bxby[3])*dv1par*dx1*wvpar+((-1.224744871391589*bxbz[1]*uz[3])-1.224744871391589*bxby[1]*uy[3]-1.224744871391589*bxbx[1]*ux[3]-1.224744871391589*uz[1]*bxbz[3]-1.224744871391589*uy[1]*bxby[3]-1.224744871391589*ux[1]*bxbx[3])*dv1par*dx0*wvpar+(1.224744871391589*Pyz[2]*rho_inv_bz[3]+1.224744871391589*Pyy[2]*rho_inv_by[3]+1.224744871391589*Pxy[2]*rho_inv_bx[3]+1.224744871391589*rho_inv_bz[2]*Pyz[3]+1.224744871391589*rho_inv_by[2]*Pyy[3]+1.224744871391589*rho_inv_bx[2]*Pxy[3])*dv1par*dx1+(1.224744871391589*Pxz[1]*rho_inv_bz[3]+1.224744871391589*Pxy[1]*rho_inv_by[3]+1.224744871391589*Pxx[1]*rho_inv_bx[3]+1.224744871391589*rho_inv_bz[1]*Pxz[3]+1.224744871391589*rho_inv_by[1]*Pxy[3]+1.224744871391589*rho_inv_bx[1]*Pxx[3])*dv1par*dx0; 
  alpha_vdim[5] = ((-0.3535533905932737*bybz[0]*uz[3])-0.3535533905932737*byby[0]*uy[3]-0.3535533905932737*bxby[0]*ux[3]-0.3535533905932737*bybz[1]*uz[2]-0.3535533905932737*byby[1]*uy[2]-0.3535533905932737*bxby[1]*ux[2])*dv1par*dvpar*dx1+((-0.3535533905932737*bxbz[3]*uz[3])-0.3535533905932737*bxby[3]*uy[3]-0.3535533905932737*bxbx[3]*ux[3]-0.3535533905932737*bxbz[1]*uz[1]-0.3535533905932737*bxby[1]*uy[1]-0.3535533905932737*bxbx[1]*ux[1])*dv1par*dvpar*dx0; 
  alpha_vdim[6] = ((-0.3535533905932737*bybz[3]*uz[3])-0.3535533905932737*byby[3]*uy[3]-0.3535533905932737*bxby[3]*ux[3]-0.3535533905932737*bybz[2]*uz[2]-0.3535533905932737*byby[2]*uy[2]-0.3535533905932737*bxby[2]*ux[2])*dv1par*dvpar*dx1+((-0.3535533905932737*bxbz[0]*uz[3])-0.3535533905932737*bxby[0]*uy[3]-0.3535533905932737*bxbx[0]*ux[3]-0.3535533905932737*uz[1]*bxbz[2]-0.3535533905932737*uy[1]*bxby[2]-0.3535533905932737*ux[1]*bxbx[2])*dv1par*dvpar*dx0; 
  alpha_vdim[7] = ((-0.3535533905932737*bybz[2]*uz[3])-0.3535533905932737*byby[2]*uy[3]-0.3535533905932737*bxby[2]*ux[3]-0.3535533905932737*uz[2]*bybz[3]-0.3535533905932737*uy[2]*byby[3]-0.3535533905932737*ux[2]*bxby[3])*dv1par*dvpar*dx1+((-0.3535533905932737*bxbz[1]*uz[3])-0.3535533905932737*bxby[1]*uy[3]-0.3535533905932737*bxbx[1]*ux[3]-0.3535533905932737*uz[1]*bxbz[3]-0.3535533905932737*uy[1]*bxby[3]-0.3535533905932737*ux[1]*bxbx[3])*dv1par*dvpar*dx0; 
  cflFreq_mid += 5.0*fabs(0.1767766952966368*alpha_vdim[0]); 

  out[1] += 0.6123724356957944*(alpha_cdim[7]*f[7]+alpha_cdim[6]*f[6]+alpha_cdim[5]*f[5]+alpha_cdim[4]*f[4]+alpha_cdim[3]*f[3]+alpha_cdim[2]*f[2]+alpha_cdim[1]*f[1]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(f[7]*alpha_cdim[19]+f[6]*alpha_cdim[18]+f[5]*alpha_cdim[17]+f[4]*alpha_cdim[16]+f[3]*alpha_cdim[15]+f[2]*alpha_cdim[14]+f[1]*alpha_cdim[13]+f[0]*alpha_cdim[12]); 
  out[3] += 0.6123724356957944*(alpha_vdim[7]*f[7]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.6123724356957944*(f[6]*alpha_cdim[19]+f[7]*alpha_cdim[18]+f[3]*alpha_cdim[17]+f[2]*alpha_cdim[16]+f[5]*alpha_cdim[15]+f[4]*alpha_cdim[14]+f[0]*alpha_cdim[13]+f[1]*alpha_cdim[12]+alpha_cdim[5]*f[7]+f[5]*alpha_cdim[7]+alpha_cdim[3]*f[6]+f[3]*alpha_cdim[6]+alpha_cdim[1]*f[4]+f[1]*alpha_cdim[4]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]); 
  out[5] += 0.5477225575051661*(alpha_cdim[7]*f[11]+alpha_cdim[6]*f[10]+alpha_cdim[5]*f[9]+alpha_cdim[3]*f[8])+0.6123724356957944*((alpha_vdim[6]+alpha_cdim[4])*f[7]+f[6]*alpha_vdim[7]+f[4]*alpha_cdim[7]+alpha_cdim[2]*f[6]+f[2]*alpha_cdim[6]+(alpha_vdim[3]+alpha_cdim[1])*f[5]+f[3]*alpha_vdim[5]+f[1]*alpha_cdim[5]+alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += (0.5477225575051661*f[11]+0.6123724356957944*f[4])*alpha_cdim[19]+(0.5477225575051661*f[10]+0.6123724356957944*f[2])*alpha_cdim[18]+0.5477225575051661*f[9]*alpha_cdim[17]+0.6123724356957944*(f[1]*alpha_cdim[17]+f[7]*alpha_cdim[16])+0.5477225575051661*f[8]*alpha_cdim[15]+0.6123724356957944*(f[0]*alpha_cdim[15]+f[6]*alpha_cdim[14]+f[5]*alpha_cdim[13]+f[3]*alpha_cdim[12]+alpha_vdim[5]*f[7]+f[5]*alpha_vdim[7]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[7] += (0.5477225575051661*f[10]+0.6123724356957944*f[2])*alpha_cdim[19]+(0.5477225575051661*f[11]+0.6123724356957944*f[4])*alpha_cdim[18]+0.5477225575051661*f[8]*alpha_cdim[17]+0.6123724356957944*(f[0]*alpha_cdim[17]+f[6]*alpha_cdim[16])+0.5477225575051661*f[9]*alpha_cdim[15]+0.6123724356957944*(f[1]*alpha_cdim[15]+f[7]*alpha_cdim[14]+f[3]*alpha_cdim[13]+f[5]*alpha_cdim[12])+0.5477225575051661*(alpha_cdim[5]*f[11]+alpha_cdim[3]*f[10]+alpha_cdim[7]*f[9]+alpha_cdim[6]*f[8])+0.6123724356957944*((alpha_vdim[3]+alpha_cdim[1])*f[7]+f[3]*alpha_vdim[7]+f[1]*alpha_cdim[7]+(alpha_vdim[5]+alpha_cdim[0])*f[6]+f[5]*alpha_vdim[6]+f[0]*alpha_cdim[6]+alpha_cdim[4]*f[5]+f[4]*(alpha_cdim[5]+alpha_vdim[0])+f[0]*alpha_vdim[4]+alpha_cdim[2]*f[3]+f[2]*(alpha_cdim[3]+alpha_vdim[1])+f[1]*alpha_vdim[2]); 
  out[8] += 1.224744871391589*(alpha_vdim[7]*f[11]+alpha_vdim[6]*f[10]+alpha_vdim[5]*f[9]+alpha_vdim[3]*f[8])+1.369306393762915*(alpha_vdim[4]*f[7]+f[4]*alpha_vdim[7]+alpha_vdim[2]*f[6]+f[2]*alpha_vdim[6]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[9] += (1.224744871391589*alpha_vdim[6]+0.6123724356957944*alpha_cdim[4])*f[11]+(1.224744871391589*alpha_vdim[7]+0.6123724356957944*alpha_cdim[2])*f[10]+(1.224744871391589*alpha_vdim[3]+0.6123724356957944*alpha_cdim[1])*f[9]+(1.224744871391589*alpha_vdim[5]+0.6123724356957944*alpha_cdim[0])*f[8]+0.5477225575051661*alpha_cdim[7]*f[7]+1.369306393762915*(alpha_vdim[2]*f[7]+f[2]*alpha_vdim[7])+0.5477225575051661*alpha_cdim[6]*f[6]+1.369306393762915*(alpha_vdim[4]*f[6]+f[4]*alpha_vdim[6])+0.5477225575051661*alpha_cdim[5]*f[5]+1.369306393762915*(alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5])+0.5477225575051661*alpha_cdim[3]*f[3]+1.369306393762915*(alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]); 
  out[10] += 0.5477225575051661*(f[7]*alpha_cdim[19]+f[6]*alpha_cdim[18]+f[5]*alpha_cdim[17])+0.6123724356957944*f[11]*alpha_cdim[16]+0.5477225575051661*f[3]*alpha_cdim[15]+0.6123724356957944*(f[10]*alpha_cdim[14]+f[9]*alpha_cdim[13]+f[8]*alpha_cdim[12])+1.224744871391589*(alpha_vdim[5]*f[11]+alpha_vdim[3]*f[10]+alpha_vdim[7]*f[9]+alpha_vdim[6]*f[8])+1.369306393762915*(alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7]+alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+alpha_vdim[4]*f[5]+f[4]*alpha_vdim[5]+alpha_vdim[2]*f[3]+f[2]*alpha_vdim[3]); 
  out[11] += 0.5477225575051661*(f[6]*alpha_cdim[19]+f[7]*alpha_cdim[18]+f[3]*alpha_cdim[17])+0.6123724356957944*f[10]*alpha_cdim[16]+0.5477225575051661*f[5]*alpha_cdim[15]+0.6123724356957944*(f[11]*alpha_cdim[14]+f[8]*alpha_cdim[13]+f[9]*alpha_cdim[12])+(1.224744871391589*alpha_vdim[3]+0.6123724356957944*alpha_cdim[1])*f[11]+(1.224744871391589*alpha_vdim[5]+0.6123724356957944*alpha_cdim[0])*f[10]+(1.224744871391589*alpha_vdim[6]+0.6123724356957944*alpha_cdim[4])*f[9]+(1.224744871391589*alpha_vdim[7]+0.6123724356957944*alpha_cdim[2])*f[8]+0.5477225575051661*alpha_cdim[5]*f[7]+1.369306393762915*(alpha_vdim[0]*f[7]+f[0]*alpha_vdim[7])+0.5477225575051661*(f[5]*alpha_cdim[7]+alpha_cdim[3]*f[6])+1.369306393762915*(alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6])+0.5477225575051661*f[3]*alpha_cdim[6]+1.369306393762915*(alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]+f[3]*alpha_vdim[4]); 

  return cflFreq_mid; 
} 
