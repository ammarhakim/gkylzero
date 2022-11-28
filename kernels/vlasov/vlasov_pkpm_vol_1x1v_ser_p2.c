#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_pkpm_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *bvar, const double *u_i, const double *bb_grad_u, const double *p_force, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // bvar:      magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // u_i:       flow velocity  // p_force:   total pressure force = 1/rho (b . div(P) + p_perp div(b)) for Euler PKPM.
  // bb_grad_u: bb : grad(u).
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dx0 = 2.0/dxv[0]; 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[3]; 
  const double *uz = &u_i[6]; 
  const double *bx = &bvar[0]; 
  const double *by = &bvar[3]; 
  const double *bz = &bvar[6]; 

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

  alpha_vdim[0] = 1.414213562373095*p_force[0]*dv1par-1.414213562373095*bb_grad_u[0]*dv1par*wvpar; 
  alpha_vdim[1] = 1.414213562373095*p_force[1]*dv1par-1.414213562373095*bb_grad_u[1]*dv1par*wvpar; 
  alpha_vdim[2] = -0.408248290463863*bb_grad_u[0]*dv1par*dvpar; 
  alpha_vdim[3] = -0.408248290463863*bb_grad_u[1]*dv1par*dvpar; 
  alpha_vdim[4] = 1.414213562373095*p_force[2]*dv1par-1.414213562373095*bb_grad_u[2]*dv1par*wvpar; 
  alpha_vdim[6] = -0.408248290463863*bb_grad_u[2]*dv1par*dvpar; 
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
