#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p2_surfx1_eval_quad.h> 
GKYL_CU_DH void euler_pkpm_surfx_1x_ser_p2(const double *w, const double *dxv, 
  const double *u_il, const double *u_ic, const double *u_ir,
  const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                       Cell-center coordinates.
  // dxv[NDIM]:                     Cell spacing.
  // u_il/u_ic/u_ir:                Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // statevecl/statevecc/statevecr: [rho ux, rho uy, rho uz, p_perp], Fluid input state vector in left/center/right cells.
  // out: Incremented output.

  const double dx1 = 2.0/dxv[0]; 

  const double *rhoux_l = &statevecl[0]; 
  const double *rhouy_l = &statevecl[3]; 
  const double *rhouz_l = &statevecl[6]; 
  const double *p_perp_l = &statevecl[9]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[3]; 
  const double *rhouz_c = &statevecc[6]; 
  const double *p_perp_c = &statevecc[9]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[3]; 
  const double *rhouz_r = &statevecr[6]; 
  const double *p_perp_r = &statevecr[9]; 

  const double *ux_l = &u_il[0]; 
  const double *ux_c = &u_ic[0]; 
  const double *ux_r = &u_ir[0]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[3]; 
  double *outrhou2 = &out[6]; 
  double *outp_perp = &out[9]; 

  double ux_l_r = ser_1x_p2_surfx1_eval_quad_node_0_r(ux_l); 
  double ux_c_l = ser_1x_p2_surfx1_eval_quad_node_0_l(ux_c); 
  double ux_c_r = ser_1x_p2_surfx1_eval_quad_node_0_r(ux_c); 
  double ux_r_l = ser_1x_p2_surfx1_eval_quad_node_0_l(ux_r); 

  double ux_max_l = fmax(fabs(ux_l_r), fabs(ux_c_l)); 
  double ux_max_r = fmax(fabs(ux_c_r), fabs(ux_r_l)); 

  double Ghat_rhoux_l = 0.7905694150420948*rhoux_l[2]*ux_max_l-0.7905694150420948*rhoux_c[2]*ux_max_l+0.6123724356957945*rhoux_l[1]*ux_max_l+0.6123724356957945*rhoux_c[1]*ux_max_l+0.3535533905932737*rhoux_l[0]*ux_max_l-0.3535533905932737*rhoux_c[0]*ux_max_l+0.7905694150420948*rhoux_l[2]*ux_l_r+0.6123724356957945*rhoux_l[1]*ux_l_r+0.3535533905932737*rhoux_l[0]*ux_l_r+0.7905694150420948*rhoux_c[2]*ux_c_l-0.6123724356957945*rhoux_c[1]*ux_c_l+0.3535533905932737*rhoux_c[0]*ux_c_l; 
  double Ghat_rhoux_r = 0.7905694150420948*rhoux_r[2]*ux_r_l-0.6123724356957945*rhoux_r[1]*ux_r_l+0.3535533905932737*rhoux_r[0]*ux_r_l-0.7905694150420948*rhoux_r[2]*ux_max_r+0.7905694150420948*rhoux_c[2]*ux_max_r+0.6123724356957945*rhoux_r[1]*ux_max_r+0.6123724356957945*rhoux_c[1]*ux_max_r-0.3535533905932737*rhoux_r[0]*ux_max_r+0.3535533905932737*rhoux_c[0]*ux_max_r+0.7905694150420948*rhoux_c[2]*ux_c_r+0.6123724356957945*rhoux_c[1]*ux_c_r+0.3535533905932737*rhoux_c[0]*ux_c_r; 
  double Ghat_rhouy_l = 0.7905694150420948*rhouy_l[2]*ux_max_l-0.7905694150420948*rhouy_c[2]*ux_max_l+0.6123724356957945*rhouy_l[1]*ux_max_l+0.6123724356957945*rhouy_c[1]*ux_max_l+0.3535533905932737*rhouy_l[0]*ux_max_l-0.3535533905932737*rhouy_c[0]*ux_max_l+0.7905694150420948*rhouy_l[2]*ux_l_r+0.6123724356957945*rhouy_l[1]*ux_l_r+0.3535533905932737*rhouy_l[0]*ux_l_r+0.7905694150420948*rhouy_c[2]*ux_c_l-0.6123724356957945*rhouy_c[1]*ux_c_l+0.3535533905932737*rhouy_c[0]*ux_c_l; 
  double Ghat_rhouy_r = 0.7905694150420948*rhouy_r[2]*ux_r_l-0.6123724356957945*rhouy_r[1]*ux_r_l+0.3535533905932737*rhouy_r[0]*ux_r_l-0.7905694150420948*rhouy_r[2]*ux_max_r+0.7905694150420948*rhouy_c[2]*ux_max_r+0.6123724356957945*rhouy_r[1]*ux_max_r+0.6123724356957945*rhouy_c[1]*ux_max_r-0.3535533905932737*rhouy_r[0]*ux_max_r+0.3535533905932737*rhouy_c[0]*ux_max_r+0.7905694150420948*rhouy_c[2]*ux_c_r+0.6123724356957945*rhouy_c[1]*ux_c_r+0.3535533905932737*rhouy_c[0]*ux_c_r; 
  double Ghat_rhouz_l = 0.7905694150420948*rhouz_l[2]*ux_max_l-0.7905694150420948*rhouz_c[2]*ux_max_l+0.6123724356957945*rhouz_l[1]*ux_max_l+0.6123724356957945*rhouz_c[1]*ux_max_l+0.3535533905932737*rhouz_l[0]*ux_max_l-0.3535533905932737*rhouz_c[0]*ux_max_l+0.7905694150420948*rhouz_l[2]*ux_l_r+0.6123724356957945*rhouz_l[1]*ux_l_r+0.3535533905932737*rhouz_l[0]*ux_l_r+0.7905694150420948*rhouz_c[2]*ux_c_l-0.6123724356957945*rhouz_c[1]*ux_c_l+0.3535533905932737*rhouz_c[0]*ux_c_l; 
  double Ghat_rhouz_r = 0.7905694150420948*rhouz_r[2]*ux_r_l-0.6123724356957945*rhouz_r[1]*ux_r_l+0.3535533905932737*rhouz_r[0]*ux_r_l-0.7905694150420948*rhouz_r[2]*ux_max_r+0.7905694150420948*rhouz_c[2]*ux_max_r+0.6123724356957945*rhouz_r[1]*ux_max_r+0.6123724356957945*rhouz_c[1]*ux_max_r-0.3535533905932737*rhouz_r[0]*ux_max_r+0.3535533905932737*rhouz_c[0]*ux_max_r+0.7905694150420948*rhouz_c[2]*ux_c_r+0.6123724356957945*rhouz_c[1]*ux_c_r+0.3535533905932737*rhouz_c[0]*ux_c_r; 
  double Ghat_p_perp_l = 0.7905694150420948*p_perp_l[2]*ux_max_l-0.7905694150420948*p_perp_c[2]*ux_max_l+0.6123724356957945*p_perp_l[1]*ux_max_l+0.6123724356957945*p_perp_c[1]*ux_max_l+0.3535533905932737*p_perp_l[0]*ux_max_l-0.3535533905932737*p_perp_c[0]*ux_max_l+0.7905694150420948*p_perp_l[2]*ux_l_r+0.6123724356957945*p_perp_l[1]*ux_l_r+0.3535533905932737*p_perp_l[0]*ux_l_r+0.7905694150420948*p_perp_c[2]*ux_c_l-0.6123724356957945*p_perp_c[1]*ux_c_l+0.3535533905932737*p_perp_c[0]*ux_c_l; 
  double Ghat_p_perp_r = 0.7905694150420948*p_perp_r[2]*ux_r_l-0.6123724356957945*p_perp_r[1]*ux_r_l+0.3535533905932737*p_perp_r[0]*ux_r_l-0.7905694150420948*p_perp_r[2]*ux_max_r+0.7905694150420948*p_perp_c[2]*ux_max_r+0.6123724356957945*p_perp_r[1]*ux_max_r+0.6123724356957945*p_perp_c[1]*ux_max_r-0.3535533905932737*p_perp_r[0]*ux_max_r+0.3535533905932737*p_perp_c[0]*ux_max_r+0.7905694150420948*p_perp_c[2]*ux_c_r+0.6123724356957945*p_perp_c[1]*ux_c_r+0.3535533905932737*p_perp_c[0]*ux_c_r; 
  outrhou0[0] += 0.7071067811865475*Ghat_rhoux_l*dx1-0.7071067811865475*Ghat_rhoux_r*dx1; 
  outrhou0[1] += (-1.224744871391589*Ghat_rhoux_r*dx1)-1.224744871391589*Ghat_rhoux_l*dx1; 
  outrhou0[2] += 1.58113883008419*Ghat_rhoux_l*dx1-1.58113883008419*Ghat_rhoux_r*dx1; 

  outrhou1[0] += 0.7071067811865475*Ghat_rhouy_l*dx1-0.7071067811865475*Ghat_rhouy_r*dx1; 
  outrhou1[1] += (-1.224744871391589*Ghat_rhouy_r*dx1)-1.224744871391589*Ghat_rhouy_l*dx1; 
  outrhou1[2] += 1.58113883008419*Ghat_rhouy_l*dx1-1.58113883008419*Ghat_rhouy_r*dx1; 

  outrhou2[0] += 0.7071067811865475*Ghat_rhouz_l*dx1-0.7071067811865475*Ghat_rhouz_r*dx1; 
  outrhou2[1] += (-1.224744871391589*Ghat_rhouz_r*dx1)-1.224744871391589*Ghat_rhouz_l*dx1; 
  outrhou2[2] += 1.58113883008419*Ghat_rhouz_l*dx1-1.58113883008419*Ghat_rhouz_r*dx1; 

  outp_perp[0] += 0.7071067811865475*Ghat_p_perp_l*dx1-0.7071067811865475*Ghat_p_perp_r*dx1; 
  outp_perp[1] += (-1.224744871391589*Ghat_p_perp_r*dx1)-1.224744871391589*Ghat_p_perp_l*dx1; 
  outp_perp[2] += 1.58113883008419*Ghat_p_perp_l*dx1-1.58113883008419*Ghat_p_perp_r*dx1; 

} 
