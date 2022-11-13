#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p1_surfx1_eval_quad.h> 
GKYL_CU_DH void euler_pkpm_surfx_1x_ser_p1(const double *w, const double *dxv, 
  const double *u_il, const double *u_ic, const double *u_ir,
  const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_il/u_ic/u_ir:  Input bulk velocity (ux,uy,uz) in left/center/right cells.
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

  const double *ux_l = &u_il[0]; 
  const double *ux_c = &u_ic[0]; 
  const double *ux_r = &u_ir[0]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[2]; 
  double *outrhou2 = &out[4]; 
  double *outp_perp = &out[6]; 

  double uxrec_l = 0.408248290463863*ux_l[1]-0.408248290463863*ux_c[1]+0.3535533905932737*ux_l[0]+0.3535533905932737*ux_c[0]; 
  double uxrec_r = (-0.408248290463863*ux_r[1])+0.408248290463863*ux_c[1]+0.3535533905932737*ux_r[0]+0.3535533905932737*ux_c[0]; 

  double Ghat_rhoux_l = 0.0; 
  double Ghat_rhoux_r = 0.0; 
  double Ghat_rhouy_l = 0.0; 
  double Ghat_rhouy_r = 0.0; 
  double Ghat_rhouz_l = 0.0; 
  double Ghat_rhouz_r = 0.0; 
  double Ghat_p_perp_l = 0.0; 
  double Ghat_p_perp_r = 0.0; 
  if (uxrec_l > 0) { 
  Ghat_rhoux_l = (1.224744871391589*rhoux_l[1]+0.7071067811865475*rhoux_l[0])*uxrec_l; 
  Ghat_rhouy_l = (1.224744871391589*rhouy_l[1]+0.7071067811865475*rhouy_l[0])*uxrec_l; 
  Ghat_rhouz_l = (1.224744871391589*rhouz_l[1]+0.7071067811865475*rhouz_l[0])*uxrec_l; 
  Ghat_p_perp_l = (1.224744871391589*p_perp_l[1]+0.7071067811865475*p_perp_l[0])*uxrec_l; 
  } else { 
  Ghat_rhoux_l = (0.7071067811865475*rhoux_c[0]-1.224744871391589*rhoux_c[1])*uxrec_l; 
  Ghat_rhouy_l = (0.7071067811865475*rhouy_c[0]-1.224744871391589*rhouy_c[1])*uxrec_l; 
  Ghat_rhouz_l = (0.7071067811865475*rhouz_c[0]-1.224744871391589*rhouz_c[1])*uxrec_l; 
  Ghat_p_perp_l = (0.7071067811865475*p_perp_c[0]-1.224744871391589*p_perp_c[1])*uxrec_l; 
  } 
  if (uxrec_r > 0) { 
  Ghat_rhoux_r = (1.224744871391589*rhoux_c[1]+0.7071067811865475*rhoux_c[0])*uxrec_r; 
  Ghat_rhouy_r = (1.224744871391589*rhouy_c[1]+0.7071067811865475*rhouy_c[0])*uxrec_r; 
  Ghat_rhouz_r = (1.224744871391589*rhouz_c[1]+0.7071067811865475*rhouz_c[0])*uxrec_r; 
  Ghat_p_perp_r = (1.224744871391589*p_perp_c[1]+0.7071067811865475*p_perp_c[0])*uxrec_r; 
  } else { 
  Ghat_rhoux_r = (0.7071067811865475*rhoux_r[0]-1.224744871391589*rhoux_r[1])*uxrec_r; 
  Ghat_rhouy_r = (0.7071067811865475*rhouy_r[0]-1.224744871391589*rhouy_r[1])*uxrec_r; 
  Ghat_rhouz_r = (0.7071067811865475*rhouz_r[0]-1.224744871391589*rhouz_r[1])*uxrec_r; 
  Ghat_p_perp_r = (0.7071067811865475*p_perp_r[0]-1.224744871391589*p_perp_r[1])*uxrec_r; 
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
