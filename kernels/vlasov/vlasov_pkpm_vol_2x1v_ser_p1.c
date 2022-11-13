#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_pkpm_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *bvar, const double *u_i, const double *bb_grad_u, const double *p_force, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // bvar:      magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // u_i:       flow velocity  // p_force:   total pressure force = 1/rho (b . div(P) + p_perp div(b)) for Euler PKPM.
  // bb_grad_u: bb : grad(u).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx0 = 2.0/dxv[0]; 
  const double dx1 = 2.0/dxv[1]; 
  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[4]; 
  const double *uz = &u_i[8]; 
  const double *bx = &bvar[0]; 
  const double *by = &bvar[4]; 
  const double *bz = &bvar[8]; 

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

  alpha_vdim[0] = 1.414213562373095*p_force[0]*dv1par-1.414213562373095*bb_grad_u[0]*dv1par*wvpar; 
  alpha_vdim[1] = 1.414213562373095*p_force[1]*dv1par-1.414213562373095*bb_grad_u[1]*dv1par*wvpar; 
  alpha_vdim[2] = 1.414213562373095*p_force[2]*dv1par-1.414213562373095*bb_grad_u[2]*dv1par*wvpar; 
  alpha_vdim[3] = -0.408248290463863*bb_grad_u[0]*dv1par*dvpar; 
  alpha_vdim[4] = 1.414213562373095*p_force[3]*dv1par-1.414213562373095*bb_grad_u[3]*dv1par*wvpar; 
  alpha_vdim[5] = -0.408248290463863*bb_grad_u[1]*dv1par*dvpar; 
  alpha_vdim[6] = -0.408248290463863*bb_grad_u[2]*dv1par*dvpar; 
  alpha_vdim[7] = -0.408248290463863*bb_grad_u[3]*dv1par*dvpar; 
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
