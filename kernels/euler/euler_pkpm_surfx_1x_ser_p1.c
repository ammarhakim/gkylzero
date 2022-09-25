#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p1_surfx1_eval_quad.h> 
GKYL_CU_DH void euler_pkpm_surfx_1x_ser_p1(const double *w, const double *dxv, 
  const double *u_i, const double *p_ijl, const double *p_ijc, const double *p_ijr,
  const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i: Input bulk velocity (ux,uy,uz) in cell being updated (ASSUMED TO BE CONTINUOUS).
  // p_ijl/p_ijc/p_ijr: Pressure tensor in left/center/right cells.
  // statevecl/statevecc/statevecr: [rho ux, rho uy, rho uz, p_perp], Fluid input state vector in left/center/right cells.
  // out: Incremented output.

  const double dx1 = 2.0/dxv[0]; 

  const double *rhoux_l = &statevecl[0]; 
  const double *rhouy_l = &statevecl[2]; 
  const double *rhouz_l = &statevecl[4]; 
  const double *p_perp_l = &statevecl[6]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[2]; 
  const double *rhouz_c = &statevecc[4]; 
  const double *p_perp_c = &statevecc[6]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[2]; 
  const double *rhouz_r = &statevecr[4]; 
  const double *p_perp_r = &statevecr[6]; 

  const double *ux_c = &u_i[0]; 
  const double *uy_c = &u_i[2]; 
  const double *uz_c = &u_i[4]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxx_r = &p_ijr[0]; 

  const double *Pxy_l = &p_ijl[2]; 
  const double *Pxy_c = &p_ijc[2]; 
  const double *Pxy_r = &p_ijr[2]; 

  const double *Pxz_l = &p_ijl[4]; 
  const double *Pxz_c = &p_ijc[4]; 
  const double *Pxz_r = &p_ijr[4]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[2]; 
  double *outrhou2 = &out[4]; 
  double *outp_perp = &out[6]; 

  double ux_c_l = ser_1x_p1_surfx1_eval_quad_node_0_l(ux_c); 
  double ux_c_r = ser_1x_p1_surfx1_eval_quad_node_0_r(ux_c); 

  double Ghat_rhoux_l = 0.0; 
  double Ghat_rhoux_r = 0.0; 
  double Ghat_rhouy_l = 0.0; 
  double Ghat_rhouy_r = 0.0; 
  double Ghat_rhouz_l = 0.0; 
  double Ghat_rhouz_r = 0.0; 
  double Ghat_p_perp_l = 0.0; 
  double Ghat_p_perp_r = 0.0; 
  if (ux_c_l > 0) { 
  Ghat_rhoux_l = (1.224744871391589*rhoux_l[1]+0.7071067811865475*rhoux_l[0])*ux_c_l+0.6123724356957945*Pxx_l[1]-0.6123724356957945*Pxx_c[1]+0.3535533905932737*Pxx_l[0]+0.3535533905932737*Pxx_c[0]; 
  Ghat_rhouy_l = (1.224744871391589*rhouy_l[1]+0.7071067811865475*rhouy_l[0])*ux_c_l+0.6123724356957945*Pxy_l[1]-0.6123724356957945*Pxy_c[1]+0.3535533905932737*Pxy_l[0]+0.3535533905932737*Pxy_c[0]; 
  Ghat_rhouz_l = (1.224744871391589*rhouz_l[1]+0.7071067811865475*rhouz_l[0])*ux_c_l+0.6123724356957945*Pxz_l[1]-0.6123724356957945*Pxz_c[1]+0.3535533905932737*Pxz_l[0]+0.3535533905932737*Pxz_c[0]; 
  Ghat_p_perp_l = (1.224744871391589*p_perp_l[1]+0.7071067811865475*p_perp_l[0])*ux_c_l; 
  } else { 
  Ghat_rhoux_l = (0.7071067811865475*rhoux_c[0]-1.224744871391589*rhoux_c[1])*ux_c_l+0.6123724356957945*Pxx_l[1]-0.6123724356957945*Pxx_c[1]+0.3535533905932737*Pxx_l[0]+0.3535533905932737*Pxx_c[0]; 
  Ghat_rhouy_l = (0.7071067811865475*rhouy_c[0]-1.224744871391589*rhouy_c[1])*ux_c_l+0.6123724356957945*Pxy_l[1]-0.6123724356957945*Pxy_c[1]+0.3535533905932737*Pxy_l[0]+0.3535533905932737*Pxy_c[0]; 
  Ghat_rhouz_l = (0.7071067811865475*rhouz_c[0]-1.224744871391589*rhouz_c[1])*ux_c_l+0.6123724356957945*Pxz_l[1]-0.6123724356957945*Pxz_c[1]+0.3535533905932737*Pxz_l[0]+0.3535533905932737*Pxz_c[0]; 
  Ghat_p_perp_l = (0.7071067811865475*p_perp_c[0]-1.224744871391589*p_perp_c[1])*ux_c_l; 
  } 
  if (ux_c_r > 0) { 
  Ghat_rhoux_r = (1.224744871391589*rhoux_c[1]+0.7071067811865475*rhoux_c[0])*ux_c_r-0.6123724356957945*Pxx_r[1]+0.6123724356957945*Pxx_c[1]+0.3535533905932737*Pxx_r[0]+0.3535533905932737*Pxx_c[0]; 
  Ghat_rhouy_r = (1.224744871391589*rhouy_c[1]+0.7071067811865475*rhouy_c[0])*ux_c_r-0.6123724356957945*Pxy_r[1]+0.6123724356957945*Pxy_c[1]+0.3535533905932737*Pxy_r[0]+0.3535533905932737*Pxy_c[0]; 
  Ghat_rhouz_r = (1.224744871391589*rhouz_c[1]+0.7071067811865475*rhouz_c[0])*ux_c_r-0.6123724356957945*Pxz_r[1]+0.6123724356957945*Pxz_c[1]+0.3535533905932737*Pxz_r[0]+0.3535533905932737*Pxz_c[0]; 
  Ghat_p_perp_r = (1.224744871391589*p_perp_c[1]+0.7071067811865475*p_perp_c[0])*ux_c_r; 
  } else { 
  Ghat_rhoux_r = (0.7071067811865475*rhoux_r[0]-1.224744871391589*rhoux_r[1])*ux_c_r-0.6123724356957945*Pxx_r[1]+0.6123724356957945*Pxx_c[1]+0.3535533905932737*Pxx_r[0]+0.3535533905932737*Pxx_c[0]; 
  Ghat_rhouy_r = (0.7071067811865475*rhouy_r[0]-1.224744871391589*rhouy_r[1])*ux_c_r-0.6123724356957945*Pxy_r[1]+0.6123724356957945*Pxy_c[1]+0.3535533905932737*Pxy_r[0]+0.3535533905932737*Pxy_c[0]; 
  Ghat_rhouz_r = (0.7071067811865475*rhouz_r[0]-1.224744871391589*rhouz_r[1])*ux_c_r-0.6123724356957945*Pxz_r[1]+0.6123724356957945*Pxz_c[1]+0.3535533905932737*Pxz_r[0]+0.3535533905932737*Pxz_c[0]; 
  Ghat_p_perp_r = (0.7071067811865475*p_perp_r[0]-1.224744871391589*p_perp_r[1])*ux_c_r; 
  } 

  outrhou0[0] += 0.7071067811865475*Ghat_rhoux_l*dx1-0.7071067811865475*Ghat_rhoux_r*dx1; 
  outrhou0[1] += (-1.224744871391589*Ghat_rhoux_r*dx1)-1.224744871391589*Ghat_rhoux_l*dx1; 

  outrhou1[0] += 0.7071067811865475*Ghat_rhouy_l*dx1-0.7071067811865475*Ghat_rhouy_r*dx1; 
  outrhou1[1] += (-1.224744871391589*Ghat_rhouy_r*dx1)-1.224744871391589*Ghat_rhouy_l*dx1; 

  outrhou2[0] += 0.7071067811865475*Ghat_rhouz_l*dx1-0.7071067811865475*Ghat_rhouz_r*dx1; 
  outrhou2[1] += (-1.224744871391589*Ghat_rhouz_r*dx1)-1.224744871391589*Ghat_rhouz_l*dx1; 

  outp_perp[0] += 0.7071067811865475*Ghat_p_perp_l*dx1-0.7071067811865475*Ghat_p_perp_r*dx1; 
  outp_perp[1] += (-1.224744871391589*Ghat_p_perp_r*dx1)-1.224744871391589*Ghat_p_perp_l*dx1; 

} 
