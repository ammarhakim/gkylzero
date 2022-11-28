#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void euler_pkpm_surfx_2x_ser_p1(const double *w, const double *dxv,
  const double *u_il, const double *u_ic, const double *u_ir,
  const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                       Cell-center coordinates.
  // dxv[NDIM]:                     Cell spacing.
  // u_il/u_ic/u_ir:                Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // statevecl/statevecc/statevecr: [rho ux, rho uy, rho uz, E_perp], Fluid input state vector in left/center/right cells.
  // out: Incremented output.

  const double dx1 = 2.0/dxv[0]; 

  const double *rhoux_l = &statevecl[0]; 
  const double *rhouy_l = &statevecl[4]; 
  const double *rhouz_l = &statevecl[8]; 
  const double *p_perp_l = &statevecl[12]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[4]; 
  const double *rhouz_c = &statevecc[8]; 
  const double *p_perp_c = &statevecc[12]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[4]; 
  const double *rhouz_r = &statevecr[8]; 
  const double *p_perp_r = &statevecr[12]; 

  const double *u_l = &u_il[0]; 
  const double *u_c = &u_ic[0]; 
  const double *u_r = &u_ir[0]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[4]; 
  double *outrhou2 = &out[8]; 
  double *outp_perp = &out[12]; 

  double uSurf_l[2] = {0.0}; 
  uSurf_l[0] = 0.408248290463863*u_l[1]-0.408248290463863*u_c[1]+0.3535533905932737*u_l[0]+0.3535533905932737*u_c[0]; 
  uSurf_l[1] = 0.408248290463863*u_l[3]-0.408248290463863*u_c[3]+0.3535533905932737*u_l[2]+0.3535533905932737*u_c[2]; 

  double uSurf_r[2] = {0.0}; 
  uSurf_r[0] = (-0.408248290463863*u_r[1])+0.408248290463863*u_c[1]+0.3535533905932737*u_r[0]+0.3535533905932737*u_c[0]; 
  uSurf_r[1] = (-0.408248290463863*u_r[3])+0.408248290463863*u_c[3]+0.3535533905932737*u_r[2]+0.3535533905932737*u_c[2]; 

  double rhouxUpwindQuad_l[2] = {0.0};
  double rhouxUpwindQuad_r[2] = {0.0};
  double rhouxUpwind_l[2] = {0.0};
  double rhouxUpwind_r[2] = {0.0};
  double Ghat_rhoux_l[2] = {0.0}; 
  double Ghat_rhoux_r[2] = {0.0}; 

  double rhouyUpwindQuad_l[2] = {0.0};
  double rhouyUpwindQuad_r[2] = {0.0};
  double rhouyUpwind_l[2] = {0.0};
  double rhouyUpwind_r[2] = {0.0};
  double Ghat_rhouy_l[2] = {0.0}; 
  double Ghat_rhouy_r[2] = {0.0}; 

  double rhouzUpwindQuad_l[2] = {0.0};
  double rhouzUpwindQuad_r[2] = {0.0};
  double rhouzUpwind_l[2] = {0.0};
  double rhouzUpwind_r[2] = {0.0};
  double Ghat_rhouz_l[2] = {0.0}; 
  double Ghat_rhouz_r[2] = {0.0}; 

  double p_perpUpwindQuad_l[2] = {0.0};
  double p_perpUpwindQuad_r[2] = {0.0};
  double p_perpUpwind_l[2] = {0.0};
  double p_perpUpwind_r[2] = {0.0};
  double Ghat_p_perp_l[2] = {0.0}; 
  double Ghat_p_perp_r[2] = {0.0}; 

  if (0.7071067811865475*uSurf_l[0]-0.7071067811865475*uSurf_l[1] > 0) { 
    rhouxUpwindQuad_l[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(rhoux_l); 
    rhouyUpwindQuad_l[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(rhouy_l); 
    rhouzUpwindQuad_l[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(rhouz_l); 
    p_perpUpwindQuad_l[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(p_perp_l); 
  } else { 
    rhouxUpwindQuad_l[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(rhoux_c); 
    rhouyUpwindQuad_l[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(rhouy_c); 
    rhouzUpwindQuad_l[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(rhouz_c); 
    p_perpUpwindQuad_l[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(p_perp_c); 
  } 
  if (0.7071067811865475*uSurf_r[0]-0.7071067811865475*uSurf_r[1] > 0) { 
    rhouxUpwindQuad_r[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(rhoux_c); 
    rhouyUpwindQuad_r[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(rhouy_c); 
    rhouzUpwindQuad_r[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(rhouz_c); 
    p_perpUpwindQuad_r[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(p_perp_c); 
  } else { 
    rhouxUpwindQuad_r[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(rhoux_r); 
    rhouyUpwindQuad_r[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(rhouy_r); 
    rhouzUpwindQuad_r[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(rhouz_r); 
    p_perpUpwindQuad_r[0] = ser_2x_p1_surfx1_eval_quad_node_0_r(p_perp_r); 
  } 
  if (0.7071067811865475*(uSurf_l[1]+uSurf_l[0]) > 0) { 
    rhouxUpwindQuad_l[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(rhoux_l); 
    rhouyUpwindQuad_l[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(rhouy_l); 
    rhouzUpwindQuad_l[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(rhouz_l); 
    p_perpUpwindQuad_l[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(p_perp_l); 
  } else { 
    rhouxUpwindQuad_l[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(rhoux_c); 
    rhouyUpwindQuad_l[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(rhouy_c); 
    rhouzUpwindQuad_l[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(rhouz_c); 
    p_perpUpwindQuad_l[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(p_perp_c); 
  } 
  if (0.7071067811865475*(uSurf_r[1]+uSurf_r[0]) > 0) { 
    rhouxUpwindQuad_r[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(rhoux_c); 
    rhouyUpwindQuad_r[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(rhouy_c); 
    rhouzUpwindQuad_r[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(rhouz_c); 
    p_perpUpwindQuad_r[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(p_perp_c); 
  } else { 
    rhouxUpwindQuad_r[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(rhoux_r); 
    rhouyUpwindQuad_r[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(rhouy_r); 
    rhouzUpwindQuad_r[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(rhouz_r); 
    p_perpUpwindQuad_r[1] = ser_2x_p1_surfx1_eval_quad_node_1_r(p_perp_r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p1_upwind_quad_to_modal(rhouxUpwindQuad_l, rhouxUpwind_l); 
  ser_2x_p1_upwind_quad_to_modal(rhouxUpwindQuad_r, rhouxUpwind_r); 
  ser_2x_p1_upwind_quad_to_modal(rhouyUpwindQuad_l, rhouyUpwind_l); 
  ser_2x_p1_upwind_quad_to_modal(rhouyUpwindQuad_r, rhouyUpwind_r); 
  ser_2x_p1_upwind_quad_to_modal(rhouzUpwindQuad_l, rhouzUpwind_l); 
  ser_2x_p1_upwind_quad_to_modal(rhouzUpwindQuad_r, rhouzUpwind_r); 
  ser_2x_p1_upwind_quad_to_modal(p_perpUpwindQuad_l, p_perpUpwind_l); 
  ser_2x_p1_upwind_quad_to_modal(p_perpUpwindQuad_r, p_perpUpwind_r); 
  Ghat_rhoux_l[0] = 0.7071067811865475*rhouxUpwind_l[1]*uSurf_l[1]+0.7071067811865475*rhouxUpwind_l[0]*uSurf_l[0]; 
  Ghat_rhoux_l[1] = 0.7071067811865475*rhouxUpwind_l[0]*uSurf_l[1]+0.7071067811865475*uSurf_l[0]*rhouxUpwind_l[1]; 

  Ghat_rhouy_l[0] = 0.7071067811865475*rhouyUpwind_l[1]*uSurf_l[1]+0.7071067811865475*rhouyUpwind_l[0]*uSurf_l[0]; 
  Ghat_rhouy_l[1] = 0.7071067811865475*rhouyUpwind_l[0]*uSurf_l[1]+0.7071067811865475*uSurf_l[0]*rhouyUpwind_l[1]; 

  Ghat_rhouz_l[0] = 0.7071067811865475*rhouzUpwind_l[1]*uSurf_l[1]+0.7071067811865475*rhouzUpwind_l[0]*uSurf_l[0]; 
  Ghat_rhouz_l[1] = 0.7071067811865475*rhouzUpwind_l[0]*uSurf_l[1]+0.7071067811865475*uSurf_l[0]*rhouzUpwind_l[1]; 

  Ghat_p_perp_l[0] = 0.7071067811865475*p_perpUpwind_l[1]*uSurf_l[1]+0.7071067811865475*p_perpUpwind_l[0]*uSurf_l[0]; 
  Ghat_p_perp_l[1] = 0.7071067811865475*p_perpUpwind_l[0]*uSurf_l[1]+0.7071067811865475*uSurf_l[0]*p_perpUpwind_l[1]; 

  Ghat_rhoux_r[0] = 0.7071067811865475*rhouxUpwind_r[1]*uSurf_r[1]+0.7071067811865475*rhouxUpwind_r[0]*uSurf_r[0]; 
  Ghat_rhoux_r[1] = 0.7071067811865475*rhouxUpwind_r[0]*uSurf_r[1]+0.7071067811865475*uSurf_r[0]*rhouxUpwind_r[1]; 

  Ghat_rhouy_r[0] = 0.7071067811865475*rhouyUpwind_r[1]*uSurf_r[1]+0.7071067811865475*rhouyUpwind_r[0]*uSurf_r[0]; 
  Ghat_rhouy_r[1] = 0.7071067811865475*rhouyUpwind_r[0]*uSurf_r[1]+0.7071067811865475*uSurf_r[0]*rhouyUpwind_r[1]; 

  Ghat_rhouz_r[0] = 0.7071067811865475*rhouzUpwind_r[1]*uSurf_r[1]+0.7071067811865475*rhouzUpwind_r[0]*uSurf_r[0]; 
  Ghat_rhouz_r[1] = 0.7071067811865475*rhouzUpwind_r[0]*uSurf_r[1]+0.7071067811865475*uSurf_r[0]*rhouzUpwind_r[1]; 

  Ghat_p_perp_r[0] = 0.7071067811865475*p_perpUpwind_r[1]*uSurf_r[1]+0.7071067811865475*p_perpUpwind_r[0]*uSurf_r[0]; 
  Ghat_p_perp_r[1] = 0.7071067811865475*p_perpUpwind_r[0]*uSurf_r[1]+0.7071067811865475*uSurf_r[0]*p_perpUpwind_r[1]; 

  outrhou0[0] += 0.7071067811865475*Ghat_rhoux_l[0]*dx1-0.7071067811865475*Ghat_rhoux_r[0]*dx1; 
  outrhou0[1] += (-1.224744871391589*Ghat_rhoux_r[0]*dx1)-1.224744871391589*Ghat_rhoux_l[0]*dx1; 
  outrhou0[2] += 0.7071067811865475*Ghat_rhoux_l[1]*dx1-0.7071067811865475*Ghat_rhoux_r[1]*dx1; 
  outrhou0[3] += (-1.224744871391589*Ghat_rhoux_r[1]*dx1)-1.224744871391589*Ghat_rhoux_l[1]*dx1; 

  outrhou1[0] += 0.7071067811865475*Ghat_rhouy_l[0]*dx1-0.7071067811865475*Ghat_rhouy_r[0]*dx1; 
  outrhou1[1] += (-1.224744871391589*Ghat_rhouy_r[0]*dx1)-1.224744871391589*Ghat_rhouy_l[0]*dx1; 
  outrhou1[2] += 0.7071067811865475*Ghat_rhouy_l[1]*dx1-0.7071067811865475*Ghat_rhouy_r[1]*dx1; 
  outrhou1[3] += (-1.224744871391589*Ghat_rhouy_r[1]*dx1)-1.224744871391589*Ghat_rhouy_l[1]*dx1; 

  outrhou2[0] += 0.7071067811865475*Ghat_rhouz_l[0]*dx1-0.7071067811865475*Ghat_rhouz_r[0]*dx1; 
  outrhou2[1] += (-1.224744871391589*Ghat_rhouz_r[0]*dx1)-1.224744871391589*Ghat_rhouz_l[0]*dx1; 
  outrhou2[2] += 0.7071067811865475*Ghat_rhouz_l[1]*dx1-0.7071067811865475*Ghat_rhouz_r[1]*dx1; 
  outrhou2[3] += (-1.224744871391589*Ghat_rhouz_r[1]*dx1)-1.224744871391589*Ghat_rhouz_l[1]*dx1; 

  outp_perp[0] += 0.7071067811865475*Ghat_p_perp_l[0]*dx1-0.7071067811865475*Ghat_p_perp_r[0]*dx1; 
  outp_perp[1] += (-1.224744871391589*Ghat_p_perp_r[0]*dx1)-1.224744871391589*Ghat_p_perp_l[0]*dx1; 
  outp_perp[2] += 0.7071067811865475*Ghat_p_perp_l[1]*dx1-0.7071067811865475*Ghat_p_perp_r[1]*dx1; 
  outp_perp[3] += (-1.224744871391589*Ghat_p_perp_r[1]*dx1)-1.224744871391589*Ghat_p_perp_l[1]*dx1; 

} 
