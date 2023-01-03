#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_pkpm_vol_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *bvar, const double *u_i, 
  const double *bb_grad_u, const double *p_force, 
  const double *div_b, const double *p_perp_source, 
  const double *p_perp_div_b, const double *g_dist_source, 
  const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates.
  // dxv[NDIM]:     Cell spacing.
  // bvar:          magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // u_i:           flow velocity  // p_force:       total pressure force; for F_0 = 1/rho (div(p_parallel b_hat), for G_1 = 1/rho (div(p_parallel b_hat) + 3*T_perp/m*div(b).
  // bb_grad_u:     bb : grad(u).
  // div_b:         divergence of the magnetic field unit vector.
  // p_perp_source: perpendicular pressure source bb : grad(u) - div(u) - nu + nu rho vth^2/p_perp.
  // p_perp_div_b:  p_perp/rho*div(b) = T_perp/m*div(b).
  // g_dist_source: [-F_1, -2 T_perp/m*(F_0 - F_2)].
  // f:             Input distribution function [F_0, T_perp/m G = T_perp/m (F_0 - F_1)].
  // out:           Incremented output.
  const double dx0 = 2.0/dxv[0]; 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[2]; 
  const double *uz = &u_i[4]; 
  const double *bx = &bvar[0]; 
  const double *by = &bvar[2]; 
  const double *bz = &bvar[4]; 

  const double *p_force_F_0 = &p_force[0]; 
  const double *p_force_G_1 = &p_force[2]; 
  const double *F_0 = &f[0]; 
  const double *G_1 = &f[6]; 
  const double *F_0_source = &g_dist_source[0]; 
  const double *G_1_source = &g_dist_source[6]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[6]; 
  double cflFreq_mid = 0.0; 
  double alpha_cdim[6] = {0.0}; 
  double alpha_vdim_F_0[6] = {0.0}; 

  double alpha_vdim_G_1[6] = {0.0}; 

  double alpha_div_b[6] = {0.0}; 

  double alpha_G_1_source[6] = {0.0}; 

  alpha_cdim[0] = 1.414213562373095*dx0*(bx[0]*wvpar+ux[0]); 
  alpha_cdim[1] = 1.414213562373095*dx0*(bx[1]*wvpar+ux[1]); 
  alpha_cdim[2] = 0.408248290463863*bx[0]*dvpar*dx0; 
  alpha_cdim[3] = 0.408248290463863*bx[1]*dvpar*dx0; 
  cflFreq_mid += 3.0*fabs(0.25*alpha_cdim[0]); 

  alpha_vdim_F_0[0] = 1.414213562373095*p_force_F_0[0]*dv1par-1.414213562373095*bb_grad_u[0]*dv1par*wvpar; 
  alpha_vdim_F_0[1] = 1.414213562373095*p_force_F_0[1]*dv1par-1.414213562373095*bb_grad_u[1]*dv1par*wvpar; 
  alpha_vdim_F_0[2] = -0.408248290463863*bb_grad_u[0]*dv1par*dvpar; 
  alpha_vdim_F_0[3] = -0.408248290463863*bb_grad_u[1]*dv1par*dvpar; 

  alpha_vdim_G_1[0] = 1.414213562373095*p_force_G_1[0]*dv1par-1.414213562373095*bb_grad_u[0]*dv1par*wvpar; 
  alpha_vdim_G_1[1] = 1.414213562373095*p_force_G_1[1]*dv1par-1.414213562373095*bb_grad_u[1]*dv1par*wvpar; 
  alpha_vdim_G_1[2] = -0.408248290463863*bb_grad_u[0]*dv1par*dvpar; 
  alpha_vdim_G_1[3] = -0.408248290463863*bb_grad_u[1]*dv1par*dvpar; 

  cflFreq_mid += 5.0*fabs(0.25*alpha_vdim_G_1[0]); 

  alpha_div_b[0] = -1.414213562373095*div_b[0]*dv1par; 
  alpha_div_b[1] = -1.414213562373095*div_b[1]*dv1par; 

  cflFreq_mid += 5.0*fabs(0.3535533905932737*p_perp_div_b[0]); 

  alpha_G_1_source[0] = 1.414213562373095*p_perp_source[0]-1.414213562373095*div_b[0]*wvpar; 
  alpha_G_1_source[1] = 1.414213562373095*p_perp_source[1]-1.414213562373095*div_b[1]*wvpar; 
  alpha_G_1_source[2] = -0.408248290463863*div_b[0]*dvpar; 
  alpha_G_1_source[3] = -0.408248290463863*div_b[1]*dvpar; 

  out_F_0[1] += 0.8660254037844386*F_0[3]*alpha_cdim[3]+0.8660254037844386*F_0[2]*alpha_cdim[2]+0.8660254037844386*F_0[1]*alpha_cdim[1]+0.8660254037844386*F_0[0]*alpha_cdim[0]; 
  out_F_0[2] += 0.8660254037844386*F_0[3]*alpha_vdim_F_0[3]+0.8660254037844386*F_0[2]*alpha_vdim_F_0[2]+0.8660254037844386*F_0[1]*alpha_vdim_F_0[1]+0.8660254037844386*F_0_source[1]*alpha_div_b[1]+0.8660254037844386*F_0[0]*alpha_vdim_F_0[0]+0.8660254037844386*F_0_source[0]*alpha_div_b[0]; 
  out_F_0[3] += 0.7745966692414834*alpha_cdim[3]*F_0[5]+0.7745966692414833*alpha_cdim[2]*F_0[4]+0.8660254037844386*F_0[2]*alpha_vdim_F_0[3]+0.8660254037844386*F_0[1]*alpha_cdim[3]+0.8660254037844386*alpha_vdim_F_0[2]*F_0[3]+0.8660254037844386*alpha_cdim[1]*F_0[3]+0.8660254037844386*F_0[0]*alpha_cdim[2]+0.8660254037844386*alpha_cdim[0]*F_0[2]+0.8660254037844386*F_0[0]*alpha_vdim_F_0[1]+0.8660254037844386*F_0_source[0]*alpha_div_b[1]+0.8660254037844386*alpha_div_b[0]*F_0_source[1]+0.8660254037844386*alpha_vdim_F_0[0]*F_0[1]; 
  out_F_0[4] += 1.732050807568877*alpha_vdim_F_0[3]*F_0[5]+1.732050807568877*alpha_vdim_F_0[2]*F_0[4]+1.936491673103709*F_0[1]*alpha_vdim_F_0[3]+1.936491673103709*alpha_div_b[1]*F_0_source[3]+1.936491673103709*alpha_vdim_F_0[1]*F_0[3]+1.936491673103709*F_0[0]*alpha_vdim_F_0[2]+1.936491673103709*alpha_div_b[0]*F_0_source[2]+1.936491673103709*alpha_vdim_F_0[0]*F_0[2]; 
  out_F_0[5] += 1.732050807568877*alpha_vdim_F_0[2]*F_0[5]+0.8660254037844386*alpha_cdim[1]*F_0[5]+1.732050807568877*alpha_vdim_F_0[3]*F_0[4]+0.8660254037844387*alpha_cdim[0]*F_0[4]+1.936491673103709*F_0[0]*alpha_vdim_F_0[3]+0.7745966692414834*F_0[3]*alpha_cdim[3]+1.936491673103709*alpha_div_b[0]*F_0_source[3]+1.936491673103709*alpha_vdim_F_0[0]*F_0[3]+1.936491673103709*F_0[1]*alpha_vdim_F_0[2]+0.7745966692414834*F_0[2]*alpha_cdim[2]+1.936491673103709*alpha_div_b[1]*F_0_source[2]+1.936491673103709*alpha_vdim_F_0[1]*F_0[2]; 
  out_G_1[0] += 0.5*G_1[3]*alpha_G_1_source[3]+0.5*G_1[2]*alpha_G_1_source[2]+0.5*G_1[1]*alpha_G_1_source[1]+0.5*G_1[0]*alpha_G_1_source[0]; 
  out_G_1[1] += 0.8660254037844386*G_1[3]*alpha_cdim[3]+0.5*G_1[2]*alpha_G_1_source[3]+0.5*alpha_G_1_source[2]*G_1[3]+0.8660254037844386*G_1[2]*alpha_cdim[2]+0.8660254037844386*G_1[1]*alpha_cdim[1]+0.5*G_1[0]*alpha_G_1_source[1]+0.5*alpha_G_1_source[0]*G_1[1]+0.8660254037844386*G_1[0]*alpha_cdim[0]; 
  out_G_1[2] += 0.447213595499958*alpha_G_1_source[3]*G_1[5]+0.4472135954999579*alpha_G_1_source[2]*G_1[4]+0.8660254037844386*G_1[3]*alpha_vdim_G_1[3]+0.5*G_1[1]*alpha_G_1_source[3]+0.5*alpha_G_1_source[1]*G_1[3]+0.8660254037844386*G_1[2]*alpha_vdim_G_1[2]+0.5*G_1[0]*alpha_G_1_source[2]+0.5*alpha_G_1_source[0]*G_1[2]+0.8660254037844386*G_1[1]*alpha_vdim_G_1[1]+0.8660254037844386*G_1_source[1]*alpha_div_b[1]+0.8660254037844386*G_1[0]*alpha_vdim_G_1[0]+0.8660254037844386*G_1_source[0]*alpha_div_b[0]; 
  out_G_1[3] += 0.7745966692414834*alpha_cdim[3]*G_1[5]+0.447213595499958*alpha_G_1_source[2]*G_1[5]+0.4472135954999579*alpha_G_1_source[3]*G_1[4]+0.7745966692414833*alpha_cdim[2]*G_1[4]+0.8660254037844386*G_1[2]*alpha_vdim_G_1[3]+0.8660254037844386*G_1[1]*alpha_cdim[3]+0.5*G_1[0]*alpha_G_1_source[3]+0.8660254037844386*alpha_vdim_G_1[2]*G_1[3]+0.8660254037844386*alpha_cdim[1]*G_1[3]+0.5*alpha_G_1_source[0]*G_1[3]+0.8660254037844386*G_1[0]*alpha_cdim[2]+0.5*G_1[1]*alpha_G_1_source[2]+0.5*alpha_G_1_source[1]*G_1[2]+0.8660254037844386*alpha_cdim[0]*G_1[2]+0.8660254037844386*G_1[0]*alpha_vdim_G_1[1]+0.8660254037844386*G_1_source[0]*alpha_div_b[1]+0.8660254037844386*alpha_div_b[0]*G_1_source[1]+0.8660254037844386*alpha_vdim_G_1[0]*G_1[1]; 
  out_G_1[4] += 1.732050807568877*alpha_vdim_G_1[3]*G_1[5]+0.5000000000000001*alpha_G_1_source[1]*G_1[5]+1.732050807568877*alpha_vdim_G_1[2]*G_1[4]+0.5*alpha_G_1_source[0]*G_1[4]+1.936491673103709*G_1[1]*alpha_vdim_G_1[3]+0.4472135954999579*G_1[3]*alpha_G_1_source[3]+1.936491673103709*alpha_div_b[1]*G_1_source[3]+1.936491673103709*alpha_vdim_G_1[1]*G_1[3]+1.936491673103709*G_1[0]*alpha_vdim_G_1[2]+0.4472135954999579*G_1[2]*alpha_G_1_source[2]+1.936491673103709*alpha_div_b[0]*G_1_source[2]+1.936491673103709*alpha_vdim_G_1[0]*G_1[2]; 
  out_G_1[5] += 1.732050807568877*alpha_vdim_G_1[2]*G_1[5]+0.8660254037844386*alpha_cdim[1]*G_1[5]+0.5*alpha_G_1_source[0]*G_1[5]+1.732050807568877*alpha_vdim_G_1[3]*G_1[4]+0.5000000000000001*alpha_G_1_source[1]*G_1[4]+0.8660254037844387*alpha_cdim[0]*G_1[4]+1.936491673103709*G_1[0]*alpha_vdim_G_1[3]+0.7745966692414834*G_1[3]*alpha_cdim[3]+0.447213595499958*G_1[2]*alpha_G_1_source[3]+1.936491673103709*alpha_div_b[0]*G_1_source[3]+0.447213595499958*alpha_G_1_source[2]*G_1[3]+1.936491673103709*alpha_vdim_G_1[0]*G_1[3]+1.936491673103709*G_1[1]*alpha_vdim_G_1[2]+0.7745966692414834*G_1[2]*alpha_cdim[2]+1.936491673103709*alpha_div_b[1]*G_1_source[2]+1.936491673103709*alpha_vdim_G_1[1]*G_1[2]; 

  return cflFreq_mid; 
} 
