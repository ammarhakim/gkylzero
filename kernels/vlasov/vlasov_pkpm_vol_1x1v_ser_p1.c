#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_pkpm_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *bvar, const double *rho_inv_b, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:       bulk flow velocity (ux, uy, uz).
  // div_p:     divergence of the pressure tensor.
  // bvar:      magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // rho_inv_b: b_i/rho (for pressure force 1/rho * b . div(P)).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx0 = 2.0/dxv[0]; 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[2]; 
  const double *uz = &u_i[4]; 

  const double *div_p_x = &div_p[0]; 
  const double *div_p_y = &div_p[2]; 
  const double *div_p_z = &div_p[4]; 

  const double *bx = &bvar[0]; 
  const double *by = &bvar[2]; 
  const double *bz = &bvar[4]; 
  const double *bxbx = &bvar[6]; 
  const double *bxby = &bvar[8]; 
  const double *bxbz = &bvar[10]; 
  const double *byby = &bvar[12]; 
  const double *bybz = &bvar[14]; 
  const double *bzbz = &bvar[16]; 

  const double *rho_inv_bx = &rho_inv_b[0]; 
  const double *rho_inv_by = &rho_inv_b[2]; 
  const double *rho_inv_bz = &rho_inv_b[4]; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[6] = {0.0}; 
  double alpha_vdim[6] = {0.0}; 

  alpha_cdim[0] = 1.414213562373095*dx0*(bx[0]*wvpar+ux[0]); 
  alpha_cdim[1] = 1.414213562373095*dx0*(bx[1]*wvpar+ux[1]); 
  alpha_cdim[2] = 0.408248290463863*bx[0]*dvpar*dx0; 
  alpha_cdim[3] = 0.408248290463863*bx[1]*dvpar*dx0; 
  cflFreq_mid += 3.0*fabs(0.25*alpha_cdim[0]); 

  alpha_vdim[0] = ((-1.732050807568877*bxbz[0]*uz[1])-1.732050807568877*bxby[0]*uy[1]-1.732050807568877*bxbx[0]*ux[1])*dv1par*dx0*wvpar+(div_p_z[1]*rho_inv_bz[1]+div_p_y[1]*rho_inv_by[1]+div_p_x[1]*rho_inv_bx[1]+div_p_z[0]*rho_inv_bz[0]+div_p_y[0]*rho_inv_by[0]+div_p_x[0]*rho_inv_bx[0])*dv1par; 
  alpha_vdim[1] = ((-1.732050807568877*bxbz[1]*uz[1])-1.732050807568877*bxby[1]*uy[1]-1.732050807568877*bxbx[1]*ux[1])*dv1par*dx0*wvpar+(div_p_z[0]*rho_inv_bz[1]+div_p_y[0]*rho_inv_by[1]+div_p_x[0]*rho_inv_bx[1]+rho_inv_bz[0]*div_p_z[1]+rho_inv_by[0]*div_p_y[1]+rho_inv_bx[0]*div_p_x[1])*dv1par; 
  alpha_vdim[2] = ((-0.5*bxbz[0]*uz[1])-0.5*bxby[0]*uy[1]-0.5*bxbx[0]*ux[1])*dv1par*dvpar*dx0; 
  alpha_vdim[3] = ((-0.5*bxbz[1]*uz[1])-0.5*bxby[1]*uy[1]-0.5*bxbx[1]*ux[1])*dv1par*dvpar*dx0; 
  cflFreq_mid += 5.0*fabs(0.25*alpha_vdim[0]); 

  out[1] += 0.8660254037844386*(alpha_cdim[3]*f[3]+alpha_cdim[2]*f[2]+alpha_cdim[1]*f[1]+alpha_cdim[0]*f[0]); 
  out[2] += 0.8660254037844386*(alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.7745966692414833*(alpha_cdim[3]*f[5]+alpha_cdim[2]*f[4])+0.8660254037844386*((alpha_vdim[2]+alpha_cdim[1])*f[3]+f[2]*alpha_vdim[3]+f[1]*alpha_cdim[3]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[4] += 1.732050807568877*(alpha_vdim[3]*f[5]+alpha_vdim[2]*f[4])+1.936491673103709*(alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[5] += (1.732050807568877*alpha_vdim[2]+0.8660254037844386*alpha_cdim[1])*f[5]+(1.732050807568877*alpha_vdim[3]+0.8660254037844386*alpha_cdim[0])*f[4]+0.7745966692414833*alpha_cdim[3]*f[3]+1.936491673103709*(alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3])+0.7745966692414833*alpha_cdim[2]*f[2]+1.936491673103709*(alpha_vdim[1]*f[2]+f[1]*alpha_vdim[2]); 

  return cflFreq_mid; 
} 
