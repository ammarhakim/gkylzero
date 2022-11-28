#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void euler_pkpm_surfy_2x_ser_p1(const double *w, const double *dxv,
  const double *u_il, const double *u_ic, const double *u_ir,
  const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                       Cell-center coordinates.
  // dxv[NDIM]:                     Cell spacing.
  // u_il/u_ic/u_ir:                Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // statevecl/statevecc/statevecr: [rho ux, rho uy, rho uz, E_perp], Fluid input state vector in left/center/right cells.
  // out: Incremented output.

  const double dx1 = 2.0/dxv[1]; 

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

  const double *u_l = &u_il[4]; 
  const double *u_c = &u_ic[4]; 
  const double *u_r = &u_ir[4]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[4]; 
  double *outrhou2 = &out[8]; 
  double *outp_perp = &out[12]; 

  double jump_rhoux_quad_l[2] = {0.0};
  double jump_rhoux_quad_r[2] = {0.0};
  double jump_rhoux_modal_l[2] = {0.0};
  double jump_rhoux_modal_r[2] = {0.0};
  double jump_rhouy_quad_l[2] = {0.0};
  double jump_rhouy_quad_r[2] = {0.0};
  double jump_rhouy_modal_l[2] = {0.0};
  double jump_rhouy_modal_r[2] = {0.0};
  double jump_rhouz_quad_l[2] = {0.0};
  double jump_rhouz_quad_r[2] = {0.0};
  double jump_rhouz_modal_l[2] = {0.0};
  double jump_rhouz_modal_r[2] = {0.0};
  double jump_p_perp_quad_l[2] = {0.0};
  double jump_p_perp_quad_r[2] = {0.0};
  double jump_p_perp_modal_l[2] = {0.0};;
  double jump_p_perp_modal_r[2] = {0.0};

  double Ghat_rhoux_l[2] = {0.0}; 
  double Ghat_rhoux_r[2] = {0.0}; 
  double Ghat_rhouy_l[2] = {0.0}; 
  double Ghat_rhouy_r[2] = {0.0}; 
  double Ghat_rhouz_l[2] = {0.0}; 
  double Ghat_rhouz_r[2] = {0.0}; 
  double Ghat_p_perp_l[2] = {0.0}; 
  double Ghat_p_perp_r[2] = {0.0}; 

  double u_l_r = 0.0; 
  double u_c_l = 0.0; 
  double u_c_r = 0.0; 
  double u_r_l = 0.0; 
  double uQuad_l = 0.0; 
  double uQuad_r = 0.0; 
  double rhoux_l_r = 0.0; 
  double rhoux_c_l = 0.0; 
  double rhoux_c_r = 0.0; 
  double rhoux_r_l = 0.0; 
  double rhouy_l_r = 0.0; 
  double rhouy_c_l = 0.0; 
  double rhouy_c_r = 0.0; 
  double rhouy_r_l = 0.0; 
  double rhouz_l_r = 0.0; 
  double rhouz_c_l = 0.0; 
  double rhouz_c_r = 0.0; 
  double rhouz_r_l = 0.0; 
  double p_perp_l_r = 0.0; 
  double p_perp_c_l = 0.0; 
  double p_perp_c_r = 0.0; 
  double p_perp_r_l = 0.0; 

  u_l_r = ser_2x_p1_surfx2_eval_quad_node_0_r(u_l); 
  u_c_l = ser_2x_p1_surfx2_eval_quad_node_0_l(u_c); 
  u_c_r = ser_2x_p1_surfx2_eval_quad_node_0_r(u_c); 
  u_r_l = ser_2x_p1_surfx2_eval_quad_node_0_l(u_r); 
  uQuad_l = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uQuad_r = fmax(fabs(u_c_r), fabs(u_r_l)); 
  rhoux_l_r = ser_2x_p1_surfx2_eval_quad_node_0_r(rhoux_l); 
  rhoux_c_l = ser_2x_p1_surfx2_eval_quad_node_0_l(rhoux_c); 
  rhoux_c_r = ser_2x_p1_surfx2_eval_quad_node_0_r(rhoux_c); 
  rhoux_r_l = ser_2x_p1_surfx2_eval_quad_node_0_l(rhoux_r); 
  rhouy_l_r = ser_2x_p1_surfx2_eval_quad_node_0_r(rhouy_l); 
  rhouy_c_l = ser_2x_p1_surfx2_eval_quad_node_0_l(rhouy_c); 
  rhouy_c_r = ser_2x_p1_surfx2_eval_quad_node_0_r(rhouy_c); 
  rhouy_r_l = ser_2x_p1_surfx2_eval_quad_node_0_l(rhouy_r); 
  rhouz_l_r = ser_2x_p1_surfx2_eval_quad_node_0_r(rhouz_l); 
  rhouz_c_l = ser_2x_p1_surfx2_eval_quad_node_0_l(rhouz_c); 
  rhouz_c_r = ser_2x_p1_surfx2_eval_quad_node_0_r(rhouz_c); 
  rhouz_r_l = ser_2x_p1_surfx2_eval_quad_node_0_l(rhouz_r); 
  p_perp_l_r = ser_2x_p1_surfx2_eval_quad_node_0_r(p_perp_l); 
  p_perp_c_l = ser_2x_p1_surfx2_eval_quad_node_0_l(p_perp_c); 
  p_perp_c_r = ser_2x_p1_surfx2_eval_quad_node_0_r(p_perp_c); 
  p_perp_r_l = ser_2x_p1_surfx2_eval_quad_node_0_l(p_perp_r); 
  jump_rhoux_quad_l[0] = 0.5*uQuad_l*(rhoux_c_l - rhoux_l_r); 
  jump_rhoux_quad_r[0] = 0.5*uQuad_r*(rhoux_r_l - rhoux_c_r); 
  jump_rhouy_quad_l[0] = 0.5*uQuad_l*(rhouy_c_l - rhouy_l_r); 
  jump_rhouy_quad_r[0] = 0.5*uQuad_r*(rhouy_r_l - rhouy_c_r); 
  jump_rhouz_quad_l[0] = 0.5*uQuad_l*(rhouz_c_l - rhouz_l_r); 
  jump_rhouz_quad_r[0] = 0.5*uQuad_r*(rhouz_r_l - rhouz_c_r); 
  jump_p_perp_quad_l[0] = 0.5*uQuad_l*(p_perp_c_l - p_perp_l_r); 
  jump_p_perp_quad_r[0] = 0.5*uQuad_r*(p_perp_r_l - p_perp_c_r); 

  u_l_r = ser_2x_p1_surfx2_eval_quad_node_1_r(u_l); 
  u_c_l = ser_2x_p1_surfx2_eval_quad_node_1_l(u_c); 
  u_c_r = ser_2x_p1_surfx2_eval_quad_node_1_r(u_c); 
  u_r_l = ser_2x_p1_surfx2_eval_quad_node_1_l(u_r); 
  uQuad_l = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uQuad_r = fmax(fabs(u_c_r), fabs(u_r_l)); 
  rhoux_l_r = ser_2x_p1_surfx2_eval_quad_node_1_r(rhoux_l); 
  rhoux_c_l = ser_2x_p1_surfx2_eval_quad_node_1_l(rhoux_c); 
  rhoux_c_r = ser_2x_p1_surfx2_eval_quad_node_1_r(rhoux_c); 
  rhoux_r_l = ser_2x_p1_surfx2_eval_quad_node_1_l(rhoux_r); 
  rhouy_l_r = ser_2x_p1_surfx2_eval_quad_node_1_r(rhouy_l); 
  rhouy_c_l = ser_2x_p1_surfx2_eval_quad_node_1_l(rhouy_c); 
  rhouy_c_r = ser_2x_p1_surfx2_eval_quad_node_1_r(rhouy_c); 
  rhouy_r_l = ser_2x_p1_surfx2_eval_quad_node_1_l(rhouy_r); 
  rhouz_l_r = ser_2x_p1_surfx2_eval_quad_node_1_r(rhouz_l); 
  rhouz_c_l = ser_2x_p1_surfx2_eval_quad_node_1_l(rhouz_c); 
  rhouz_c_r = ser_2x_p1_surfx2_eval_quad_node_1_r(rhouz_c); 
  rhouz_r_l = ser_2x_p1_surfx2_eval_quad_node_1_l(rhouz_r); 
  p_perp_l_r = ser_2x_p1_surfx2_eval_quad_node_1_r(p_perp_l); 
  p_perp_c_l = ser_2x_p1_surfx2_eval_quad_node_1_l(p_perp_c); 
  p_perp_c_r = ser_2x_p1_surfx2_eval_quad_node_1_r(p_perp_c); 
  p_perp_r_l = ser_2x_p1_surfx2_eval_quad_node_1_l(p_perp_r); 
  jump_rhoux_quad_l[1] = 0.5*uQuad_l*(rhoux_c_l - rhoux_l_r); 
  jump_rhoux_quad_r[1] = 0.5*uQuad_r*(rhoux_r_l - rhoux_c_r); 
  jump_rhouy_quad_l[1] = 0.5*uQuad_l*(rhouy_c_l - rhouy_l_r); 
  jump_rhouy_quad_r[1] = 0.5*uQuad_r*(rhouy_r_l - rhouy_c_r); 
  jump_rhouz_quad_l[1] = 0.5*uQuad_l*(rhouz_c_l - rhouz_l_r); 
  jump_rhouz_quad_r[1] = 0.5*uQuad_r*(rhouz_r_l - rhouz_c_r); 
  jump_p_perp_quad_l[1] = 0.5*uQuad_l*(p_perp_c_l - p_perp_l_r); 
  jump_p_perp_quad_r[1] = 0.5*uQuad_r*(p_perp_r_l - p_perp_c_r); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p1_upwind_quad_to_modal(jump_rhoux_quad_l, jump_rhoux_modal_l); 
  ser_2x_p1_upwind_quad_to_modal(jump_rhoux_quad_r, jump_rhoux_modal_r); 
  ser_2x_p1_upwind_quad_to_modal(jump_rhouy_quad_l, jump_rhouy_modal_l); 
  ser_2x_p1_upwind_quad_to_modal(jump_rhouy_quad_r, jump_rhouy_modal_r); 
  ser_2x_p1_upwind_quad_to_modal(jump_rhouz_quad_l, jump_rhouz_modal_l); 
  ser_2x_p1_upwind_quad_to_modal(jump_rhouz_quad_r, jump_rhouz_modal_r); 
  ser_2x_p1_upwind_quad_to_modal(jump_p_perp_quad_l, jump_p_perp_modal_l); 
  ser_2x_p1_upwind_quad_to_modal(jump_p_perp_quad_r, jump_p_perp_modal_r); 

  Ghat_rhoux_l[0] = 0.5303300858899105*rhoux_l[3]*u_l[3]+0.3061862178478971*rhoux_l[1]*u_l[3]+0.5303300858899105*rhoux_c[3]*u_c[3]-0.3061862178478971*rhoux_c[1]*u_c[3]+0.3061862178478971*u_l[1]*rhoux_l[3]-0.3061862178478971*u_c[1]*rhoux_c[3]+0.5303300858899105*rhoux_l[2]*u_l[2]+0.3061862178478971*rhoux_l[0]*u_l[2]+0.5303300858899105*rhoux_c[2]*u_c[2]-0.3061862178478971*rhoux_c[0]*u_c[2]+0.3061862178478971*u_l[0]*rhoux_l[2]-0.3061862178478971*u_c[0]*rhoux_c[2]+0.1767766952966368*rhoux_l[1]*u_l[1]+0.1767766952966368*rhoux_c[1]*u_c[1]+0.1767766952966368*rhoux_l[0]*u_l[0]+0.1767766952966368*rhoux_c[0]*u_c[0]-1.0*jump_rhoux_modal_l[0]; 
  Ghat_rhoux_l[1] = 0.5303300858899105*rhoux_l[2]*u_l[3]+0.3061862178478971*rhoux_l[0]*u_l[3]+0.5303300858899105*rhoux_c[2]*u_c[3]-0.3061862178478971*rhoux_c[0]*u_c[3]+0.5303300858899105*u_l[2]*rhoux_l[3]+0.3061862178478971*u_l[0]*rhoux_l[3]+0.5303300858899105*u_c[2]*rhoux_c[3]-0.3061862178478971*u_c[0]*rhoux_c[3]+0.3061862178478971*rhoux_l[1]*u_l[2]-0.3061862178478971*rhoux_c[1]*u_c[2]+0.3061862178478971*u_l[1]*rhoux_l[2]-0.3061862178478971*u_c[1]*rhoux_c[2]+0.1767766952966368*rhoux_l[0]*u_l[1]+0.1767766952966368*rhoux_c[0]*u_c[1]+0.1767766952966368*u_l[0]*rhoux_l[1]+0.1767766952966368*u_c[0]*rhoux_c[1]-1.0*jump_rhoux_modal_l[1]; 

  Ghat_rhoux_r[0] = 0.5303300858899105*rhoux_r[3]*u_r[3]-0.3061862178478971*rhoux_r[1]*u_r[3]+0.5303300858899105*rhoux_c[3]*u_c[3]+0.3061862178478971*rhoux_c[1]*u_c[3]-0.3061862178478971*u_r[1]*rhoux_r[3]+0.3061862178478971*u_c[1]*rhoux_c[3]+0.5303300858899105*rhoux_r[2]*u_r[2]-0.3061862178478971*rhoux_r[0]*u_r[2]+0.5303300858899105*rhoux_c[2]*u_c[2]+0.3061862178478971*rhoux_c[0]*u_c[2]-0.3061862178478971*u_r[0]*rhoux_r[2]+0.3061862178478971*u_c[0]*rhoux_c[2]+0.1767766952966368*rhoux_r[1]*u_r[1]+0.1767766952966368*rhoux_c[1]*u_c[1]+0.1767766952966368*rhoux_r[0]*u_r[0]+0.1767766952966368*rhoux_c[0]*u_c[0]-1.0*jump_rhoux_modal_r[0]; 
  Ghat_rhoux_r[1] = 0.5303300858899105*rhoux_r[2]*u_r[3]-0.3061862178478971*rhoux_r[0]*u_r[3]+0.5303300858899105*rhoux_c[2]*u_c[3]+0.3061862178478971*rhoux_c[0]*u_c[3]+0.5303300858899105*u_r[2]*rhoux_r[3]-0.3061862178478971*u_r[0]*rhoux_r[3]+0.5303300858899105*u_c[2]*rhoux_c[3]+0.3061862178478971*u_c[0]*rhoux_c[3]-0.3061862178478971*rhoux_r[1]*u_r[2]+0.3061862178478971*rhoux_c[1]*u_c[2]-0.3061862178478971*u_r[1]*rhoux_r[2]+0.3061862178478971*u_c[1]*rhoux_c[2]+0.1767766952966368*rhoux_r[0]*u_r[1]+0.1767766952966368*rhoux_c[0]*u_c[1]+0.1767766952966368*u_r[0]*rhoux_r[1]+0.1767766952966368*u_c[0]*rhoux_c[1]-1.0*jump_rhoux_modal_r[1]; 

  Ghat_rhouy_l[0] = 0.5303300858899105*rhouy_l[3]*u_l[3]+0.3061862178478971*rhouy_l[1]*u_l[3]+0.5303300858899105*rhouy_c[3]*u_c[3]-0.3061862178478971*rhouy_c[1]*u_c[3]+0.3061862178478971*u_l[1]*rhouy_l[3]-0.3061862178478971*u_c[1]*rhouy_c[3]+0.5303300858899105*rhouy_l[2]*u_l[2]+0.3061862178478971*rhouy_l[0]*u_l[2]+0.5303300858899105*rhouy_c[2]*u_c[2]-0.3061862178478971*rhouy_c[0]*u_c[2]+0.3061862178478971*u_l[0]*rhouy_l[2]-0.3061862178478971*u_c[0]*rhouy_c[2]+0.1767766952966368*rhouy_l[1]*u_l[1]+0.1767766952966368*rhouy_c[1]*u_c[1]+0.1767766952966368*rhouy_l[0]*u_l[0]+0.1767766952966368*rhouy_c[0]*u_c[0]-1.0*jump_rhouy_modal_l[0]; 
  Ghat_rhouy_l[1] = 0.5303300858899105*rhouy_l[2]*u_l[3]+0.3061862178478971*rhouy_l[0]*u_l[3]+0.5303300858899105*rhouy_c[2]*u_c[3]-0.3061862178478971*rhouy_c[0]*u_c[3]+0.5303300858899105*u_l[2]*rhouy_l[3]+0.3061862178478971*u_l[0]*rhouy_l[3]+0.5303300858899105*u_c[2]*rhouy_c[3]-0.3061862178478971*u_c[0]*rhouy_c[3]+0.3061862178478971*rhouy_l[1]*u_l[2]-0.3061862178478971*rhouy_c[1]*u_c[2]+0.3061862178478971*u_l[1]*rhouy_l[2]-0.3061862178478971*u_c[1]*rhouy_c[2]+0.1767766952966368*rhouy_l[0]*u_l[1]+0.1767766952966368*rhouy_c[0]*u_c[1]+0.1767766952966368*u_l[0]*rhouy_l[1]+0.1767766952966368*u_c[0]*rhouy_c[1]-1.0*jump_rhouy_modal_l[1]; 

  Ghat_rhouy_r[0] = 0.5303300858899105*rhouy_r[3]*u_r[3]-0.3061862178478971*rhouy_r[1]*u_r[3]+0.5303300858899105*rhouy_c[3]*u_c[3]+0.3061862178478971*rhouy_c[1]*u_c[3]-0.3061862178478971*u_r[1]*rhouy_r[3]+0.3061862178478971*u_c[1]*rhouy_c[3]+0.5303300858899105*rhouy_r[2]*u_r[2]-0.3061862178478971*rhouy_r[0]*u_r[2]+0.5303300858899105*rhouy_c[2]*u_c[2]+0.3061862178478971*rhouy_c[0]*u_c[2]-0.3061862178478971*u_r[0]*rhouy_r[2]+0.3061862178478971*u_c[0]*rhouy_c[2]+0.1767766952966368*rhouy_r[1]*u_r[1]+0.1767766952966368*rhouy_c[1]*u_c[1]+0.1767766952966368*rhouy_r[0]*u_r[0]+0.1767766952966368*rhouy_c[0]*u_c[0]-1.0*jump_rhouy_modal_r[0]; 
  Ghat_rhouy_r[1] = 0.5303300858899105*rhouy_r[2]*u_r[3]-0.3061862178478971*rhouy_r[0]*u_r[3]+0.5303300858899105*rhouy_c[2]*u_c[3]+0.3061862178478971*rhouy_c[0]*u_c[3]+0.5303300858899105*u_r[2]*rhouy_r[3]-0.3061862178478971*u_r[0]*rhouy_r[3]+0.5303300858899105*u_c[2]*rhouy_c[3]+0.3061862178478971*u_c[0]*rhouy_c[3]-0.3061862178478971*rhouy_r[1]*u_r[2]+0.3061862178478971*rhouy_c[1]*u_c[2]-0.3061862178478971*u_r[1]*rhouy_r[2]+0.3061862178478971*u_c[1]*rhouy_c[2]+0.1767766952966368*rhouy_r[0]*u_r[1]+0.1767766952966368*rhouy_c[0]*u_c[1]+0.1767766952966368*u_r[0]*rhouy_r[1]+0.1767766952966368*u_c[0]*rhouy_c[1]-1.0*jump_rhouy_modal_r[1]; 

  Ghat_rhouz_l[0] = 0.5303300858899105*rhouz_l[3]*u_l[3]+0.3061862178478971*rhouz_l[1]*u_l[3]+0.5303300858899105*rhouz_c[3]*u_c[3]-0.3061862178478971*rhouz_c[1]*u_c[3]+0.3061862178478971*u_l[1]*rhouz_l[3]-0.3061862178478971*u_c[1]*rhouz_c[3]+0.5303300858899105*rhouz_l[2]*u_l[2]+0.3061862178478971*rhouz_l[0]*u_l[2]+0.5303300858899105*rhouz_c[2]*u_c[2]-0.3061862178478971*rhouz_c[0]*u_c[2]+0.3061862178478971*u_l[0]*rhouz_l[2]-0.3061862178478971*u_c[0]*rhouz_c[2]+0.1767766952966368*rhouz_l[1]*u_l[1]+0.1767766952966368*rhouz_c[1]*u_c[1]+0.1767766952966368*rhouz_l[0]*u_l[0]+0.1767766952966368*rhouz_c[0]*u_c[0]-1.0*jump_rhouz_modal_l[0]; 
  Ghat_rhouz_l[1] = 0.5303300858899105*rhouz_l[2]*u_l[3]+0.3061862178478971*rhouz_l[0]*u_l[3]+0.5303300858899105*rhouz_c[2]*u_c[3]-0.3061862178478971*rhouz_c[0]*u_c[3]+0.5303300858899105*u_l[2]*rhouz_l[3]+0.3061862178478971*u_l[0]*rhouz_l[3]+0.5303300858899105*u_c[2]*rhouz_c[3]-0.3061862178478971*u_c[0]*rhouz_c[3]+0.3061862178478971*rhouz_l[1]*u_l[2]-0.3061862178478971*rhouz_c[1]*u_c[2]+0.3061862178478971*u_l[1]*rhouz_l[2]-0.3061862178478971*u_c[1]*rhouz_c[2]+0.1767766952966368*rhouz_l[0]*u_l[1]+0.1767766952966368*rhouz_c[0]*u_c[1]+0.1767766952966368*u_l[0]*rhouz_l[1]+0.1767766952966368*u_c[0]*rhouz_c[1]-1.0*jump_rhouz_modal_l[1]; 

  Ghat_rhouz_r[0] = 0.5303300858899105*rhouz_r[3]*u_r[3]-0.3061862178478971*rhouz_r[1]*u_r[3]+0.5303300858899105*rhouz_c[3]*u_c[3]+0.3061862178478971*rhouz_c[1]*u_c[3]-0.3061862178478971*u_r[1]*rhouz_r[3]+0.3061862178478971*u_c[1]*rhouz_c[3]+0.5303300858899105*rhouz_r[2]*u_r[2]-0.3061862178478971*rhouz_r[0]*u_r[2]+0.5303300858899105*rhouz_c[2]*u_c[2]+0.3061862178478971*rhouz_c[0]*u_c[2]-0.3061862178478971*u_r[0]*rhouz_r[2]+0.3061862178478971*u_c[0]*rhouz_c[2]+0.1767766952966368*rhouz_r[1]*u_r[1]+0.1767766952966368*rhouz_c[1]*u_c[1]+0.1767766952966368*rhouz_r[0]*u_r[0]+0.1767766952966368*rhouz_c[0]*u_c[0]-1.0*jump_rhouz_modal_r[0]; 
  Ghat_rhouz_r[1] = 0.5303300858899105*rhouz_r[2]*u_r[3]-0.3061862178478971*rhouz_r[0]*u_r[3]+0.5303300858899105*rhouz_c[2]*u_c[3]+0.3061862178478971*rhouz_c[0]*u_c[3]+0.5303300858899105*u_r[2]*rhouz_r[3]-0.3061862178478971*u_r[0]*rhouz_r[3]+0.5303300858899105*u_c[2]*rhouz_c[3]+0.3061862178478971*u_c[0]*rhouz_c[3]-0.3061862178478971*rhouz_r[1]*u_r[2]+0.3061862178478971*rhouz_c[1]*u_c[2]-0.3061862178478971*u_r[1]*rhouz_r[2]+0.3061862178478971*u_c[1]*rhouz_c[2]+0.1767766952966368*rhouz_r[0]*u_r[1]+0.1767766952966368*rhouz_c[0]*u_c[1]+0.1767766952966368*u_r[0]*rhouz_r[1]+0.1767766952966368*u_c[0]*rhouz_c[1]-1.0*jump_rhouz_modal_r[1]; 

  Ghat_p_perp_l[0] = 0.5303300858899105*p_perp_l[3]*u_l[3]+0.3061862178478971*p_perp_l[1]*u_l[3]+0.5303300858899105*p_perp_c[3]*u_c[3]-0.3061862178478971*p_perp_c[1]*u_c[3]+0.3061862178478971*u_l[1]*p_perp_l[3]-0.3061862178478971*u_c[1]*p_perp_c[3]+0.5303300858899105*p_perp_l[2]*u_l[2]+0.3061862178478971*p_perp_l[0]*u_l[2]+0.5303300858899105*p_perp_c[2]*u_c[2]-0.3061862178478971*p_perp_c[0]*u_c[2]+0.3061862178478971*u_l[0]*p_perp_l[2]-0.3061862178478971*u_c[0]*p_perp_c[2]+0.1767766952966368*p_perp_l[1]*u_l[1]+0.1767766952966368*p_perp_c[1]*u_c[1]+0.1767766952966368*p_perp_l[0]*u_l[0]+0.1767766952966368*p_perp_c[0]*u_c[0]-1.0*jump_p_perp_modal_l[0]; 
  Ghat_p_perp_l[1] = 0.5303300858899105*p_perp_l[2]*u_l[3]+0.3061862178478971*p_perp_l[0]*u_l[3]+0.5303300858899105*p_perp_c[2]*u_c[3]-0.3061862178478971*p_perp_c[0]*u_c[3]+0.5303300858899105*u_l[2]*p_perp_l[3]+0.3061862178478971*u_l[0]*p_perp_l[3]+0.5303300858899105*u_c[2]*p_perp_c[3]-0.3061862178478971*u_c[0]*p_perp_c[3]+0.3061862178478971*p_perp_l[1]*u_l[2]-0.3061862178478971*p_perp_c[1]*u_c[2]+0.3061862178478971*u_l[1]*p_perp_l[2]-0.3061862178478971*u_c[1]*p_perp_c[2]+0.1767766952966368*p_perp_l[0]*u_l[1]+0.1767766952966368*p_perp_c[0]*u_c[1]+0.1767766952966368*u_l[0]*p_perp_l[1]+0.1767766952966368*u_c[0]*p_perp_c[1]-1.0*jump_p_perp_modal_l[1]; 

  Ghat_p_perp_r[0] = 0.5303300858899105*p_perp_r[3]*u_r[3]-0.3061862178478971*p_perp_r[1]*u_r[3]+0.5303300858899105*p_perp_c[3]*u_c[3]+0.3061862178478971*p_perp_c[1]*u_c[3]-0.3061862178478971*u_r[1]*p_perp_r[3]+0.3061862178478971*u_c[1]*p_perp_c[3]+0.5303300858899105*p_perp_r[2]*u_r[2]-0.3061862178478971*p_perp_r[0]*u_r[2]+0.5303300858899105*p_perp_c[2]*u_c[2]+0.3061862178478971*p_perp_c[0]*u_c[2]-0.3061862178478971*u_r[0]*p_perp_r[2]+0.3061862178478971*u_c[0]*p_perp_c[2]+0.1767766952966368*p_perp_r[1]*u_r[1]+0.1767766952966368*p_perp_c[1]*u_c[1]+0.1767766952966368*p_perp_r[0]*u_r[0]+0.1767766952966368*p_perp_c[0]*u_c[0]-1.0*jump_p_perp_modal_r[0]; 
  Ghat_p_perp_r[1] = 0.5303300858899105*p_perp_r[2]*u_r[3]-0.3061862178478971*p_perp_r[0]*u_r[3]+0.5303300858899105*p_perp_c[2]*u_c[3]+0.3061862178478971*p_perp_c[0]*u_c[3]+0.5303300858899105*u_r[2]*p_perp_r[3]-0.3061862178478971*u_r[0]*p_perp_r[3]+0.5303300858899105*u_c[2]*p_perp_c[3]+0.3061862178478971*u_c[0]*p_perp_c[3]-0.3061862178478971*p_perp_r[1]*u_r[2]+0.3061862178478971*p_perp_c[1]*u_c[2]-0.3061862178478971*u_r[1]*p_perp_r[2]+0.3061862178478971*u_c[1]*p_perp_c[2]+0.1767766952966368*p_perp_r[0]*u_r[1]+0.1767766952966368*p_perp_c[0]*u_c[1]+0.1767766952966368*u_r[0]*p_perp_r[1]+0.1767766952966368*u_c[0]*p_perp_c[1]-1.0*jump_p_perp_modal_r[1]; 

  outrhou0[0] += 0.7071067811865475*Ghat_rhoux_l[0]*dx1-0.7071067811865475*Ghat_rhoux_r[0]*dx1; 
  outrhou0[1] += 0.7071067811865475*Ghat_rhoux_l[1]*dx1-0.7071067811865475*Ghat_rhoux_r[1]*dx1; 
  outrhou0[2] += (-1.224744871391589*Ghat_rhoux_r[0]*dx1)-1.224744871391589*Ghat_rhoux_l[0]*dx1; 
  outrhou0[3] += (-1.224744871391589*Ghat_rhoux_r[1]*dx1)-1.224744871391589*Ghat_rhoux_l[1]*dx1; 

  outrhou1[0] += 0.7071067811865475*Ghat_rhouy_l[0]*dx1-0.7071067811865475*Ghat_rhouy_r[0]*dx1; 
  outrhou1[1] += 0.7071067811865475*Ghat_rhouy_l[1]*dx1-0.7071067811865475*Ghat_rhouy_r[1]*dx1; 
  outrhou1[2] += (-1.224744871391589*Ghat_rhouy_r[0]*dx1)-1.224744871391589*Ghat_rhouy_l[0]*dx1; 
  outrhou1[3] += (-1.224744871391589*Ghat_rhouy_r[1]*dx1)-1.224744871391589*Ghat_rhouy_l[1]*dx1; 

  outrhou2[0] += 0.7071067811865475*Ghat_rhouz_l[0]*dx1-0.7071067811865475*Ghat_rhouz_r[0]*dx1; 
  outrhou2[1] += 0.7071067811865475*Ghat_rhouz_l[1]*dx1-0.7071067811865475*Ghat_rhouz_r[1]*dx1; 
  outrhou2[2] += (-1.224744871391589*Ghat_rhouz_r[0]*dx1)-1.224744871391589*Ghat_rhouz_l[0]*dx1; 
  outrhou2[3] += (-1.224744871391589*Ghat_rhouz_r[1]*dx1)-1.224744871391589*Ghat_rhouz_l[1]*dx1; 

  outp_perp[0] += 0.7071067811865475*Ghat_p_perp_l[0]*dx1-0.7071067811865475*Ghat_p_perp_r[0]*dx1; 
  outp_perp[1] += 0.7071067811865475*Ghat_p_perp_l[1]*dx1-0.7071067811865475*Ghat_p_perp_r[1]*dx1; 
  outp_perp[2] += (-1.224744871391589*Ghat_p_perp_r[0]*dx1)-1.224744871391589*Ghat_p_perp_l[0]*dx1; 
  outp_perp[3] += (-1.224744871391589*Ghat_p_perp_r[1]*dx1)-1.224744871391589*Ghat_p_perp_l[1]*dx1; 

} 
