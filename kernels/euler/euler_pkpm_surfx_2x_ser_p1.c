#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void euler_pkpm_surfx_2x_ser_p1(const double *w, const double *dxv,
  const double *u_il, const double *u_ic, const double *u_ir,
  const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
  const double *p_perpl, const double *p_perpc, const double *p_perpr, 
  const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                       Cell-center coordinates.
  // dxv[NDIM]:                     Cell spacing.
  // u_il/u_ic/u_ir:                Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // u_perp_il/u_perp_ic/u_perp_ir: Input perpendicular bulk velocity [u_perp_x, u_perp_y, u_perp_z] in left/center/right cells.
  // p_perpl/p_perpc/p_perpr:       Input perpendicular pressure in left/center/right cells.
  // statevecl/statevecc/statevecr: [rho ux, rho uy, rho uz, E_perp], Fluid input state vector in left/center/right cells.
  // out: Incremented output.

  const double dx1 = 2.0/dxv[0]; 

  const double *rhoux_l = &statevecl[0]; 
  const double *rhouy_l = &statevecl[4]; 
  const double *rhouz_l = &statevecl[8]; 
  const double *E_perp_l = &statevecl[12]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[4]; 
  const double *rhouz_c = &statevecc[8]; 
  const double *E_perp_c = &statevecc[12]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[4]; 
  const double *rhouz_r = &statevecr[8]; 
  const double *E_perp_r = &statevecr[12]; 

  const double *u_l = &u_il[0]; 
  const double *u_c = &u_ic[0]; 
  const double *u_r = &u_ir[0]; 

  const double *u_perp_l = &u_perp_il[0]; 
  const double *u_perp_c = &u_perp_ic[0]; 
  const double *u_perp_r = &u_perp_ir[0]; 

  const double *p_perp_l = &p_perpl[0]; 
  const double *p_perp_c = &p_perpc[0]; 
  const double *p_perp_r = &p_perpr[0]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[4]; 
  double *outrhou2 = &out[8]; 
  double *outE_perp = &out[12]; 

  double uQuad_l[2] = {0.0};
  double uQuad_r[2] = {0.0};
  double uMax_l[2] = {0.0};;
  double uMax_r[2] = {0.0};
  double u_perpQuad_l[2] = {0.0};
  double u_perpQuad_r[2] = {0.0};
  double u_perpMax_l[2] = {0.0};;
  double u_perpMax_r[2] = {0.0};

  double Ghat_rhoux_l[2] = {0.0}; 
  double Ghat_rhoux_r[2] = {0.0}; 
  double Ghat_rhouy_l[2] = {0.0}; 
  double Ghat_rhouy_r[2] = {0.0}; 
  double Ghat_rhouz_l[2] = {0.0}; 
  double Ghat_rhouz_r[2] = {0.0}; 
  double Ghat_E_perp_l[2] = {0.0}; 
  double Ghat_E_perp_r[2] = {0.0}; 

  double u_l_r = 0.0; 
  double u_c_l = 0.0; 
  double u_c_r = 0.0; 
  double u_r_l = 0.0; 
  double u_perp_l_r = 0.0; 
  double u_perp_c_l = 0.0; 
  double u_perp_c_r = 0.0; 
  double u_perp_r_l = 0.0; 

  u_l_r = ser_2x_p1_surfx1_eval_quad_node_0_r(u_l); 
  u_c_l = ser_2x_p1_surfx1_eval_quad_node_0_l(u_c); 
  u_c_r = ser_2x_p1_surfx1_eval_quad_node_0_r(u_c); 
  u_r_l = ser_2x_p1_surfx1_eval_quad_node_0_l(u_r); 
  uQuad_l[0] = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uQuad_r[0] = fmax(fabs(u_c_r), fabs(u_r_l)); 

  u_perp_l_r = ser_2x_p1_surfx1_eval_quad_node_0_r(u_perp_l); 
  u_perp_c_l = ser_2x_p1_surfx1_eval_quad_node_0_l(u_perp_c); 
  u_perp_c_r = ser_2x_p1_surfx1_eval_quad_node_0_r(u_perp_c); 
  u_perp_r_l = ser_2x_p1_surfx1_eval_quad_node_0_l(u_perp_r); 
  u_perpQuad_l[0] = fmax(fabs(u_perp_l_r), fabs(u_perp_c_l)); 
  u_perpQuad_r[0] = fmax(fabs(u_perp_c_r), fabs(u_perp_r_l)); 

  u_l_r = ser_2x_p1_surfx1_eval_quad_node_1_r(u_l); 
  u_c_l = ser_2x_p1_surfx1_eval_quad_node_1_l(u_c); 
  u_c_r = ser_2x_p1_surfx1_eval_quad_node_1_r(u_c); 
  u_r_l = ser_2x_p1_surfx1_eval_quad_node_1_l(u_r); 
  uQuad_l[1] = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uQuad_r[1] = fmax(fabs(u_c_r), fabs(u_r_l)); 

  u_perp_l_r = ser_2x_p1_surfx1_eval_quad_node_1_r(u_perp_l); 
  u_perp_c_l = ser_2x_p1_surfx1_eval_quad_node_1_l(u_perp_c); 
  u_perp_c_r = ser_2x_p1_surfx1_eval_quad_node_1_r(u_perp_c); 
  u_perp_r_l = ser_2x_p1_surfx1_eval_quad_node_1_l(u_perp_r); 
  u_perpQuad_l[1] = fmax(fabs(u_perp_l_r), fabs(u_perp_c_l)); 
  u_perpQuad_r[1] = fmax(fabs(u_perp_c_r), fabs(u_perp_r_l)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p1_upwind_quad_to_modal(uQuad_l, uMax_l); 
  ser_2x_p1_upwind_quad_to_modal(uQuad_r, uMax_r); 

  ser_2x_p1_upwind_quad_to_modal(u_perpQuad_l, u_perpMax_l); 
  ser_2x_p1_upwind_quad_to_modal(u_perpQuad_r, u_perpMax_r); 

  Ghat_rhoux_l[0] = 0.5303300858899105*rhoux_l[3]*u_l[3]+0.3061862178478971*rhoux_l[2]*u_l[3]+0.5303300858899105*rhoux_c[3]*u_c[3]-0.3061862178478971*rhoux_c[2]*u_c[3]+0.3061862178478971*u_l[2]*rhoux_l[3]+0.4330127018922193*uMax_l[1]*rhoux_l[3]-0.3061862178478971*u_c[2]*rhoux_c[3]+0.4330127018922193*uMax_l[1]*rhoux_c[3]+0.1767766952966368*rhoux_l[2]*u_l[2]+0.1767766952966368*rhoux_c[2]*u_c[2]+0.25*uMax_l[1]*rhoux_l[2]-0.25*uMax_l[1]*rhoux_c[2]+0.5303300858899105*rhoux_l[1]*u_l[1]+0.3061862178478971*rhoux_l[0]*u_l[1]+0.5303300858899105*rhoux_c[1]*u_c[1]-0.3061862178478971*rhoux_c[0]*u_c[1]+0.3061862178478971*u_l[0]*rhoux_l[1]+0.4330127018922193*uMax_l[0]*rhoux_l[1]-0.3061862178478971*u_c[0]*rhoux_c[1]+0.4330127018922193*uMax_l[0]*rhoux_c[1]+0.1767766952966368*rhoux_l[0]*u_l[0]+0.1767766952966368*rhoux_c[0]*u_c[0]+0.25*rhoux_l[0]*uMax_l[0]-0.25*rhoux_c[0]*uMax_l[0]; 
  Ghat_rhoux_l[1] = 0.5303300858899105*rhoux_l[1]*u_l[3]+0.3061862178478971*rhoux_l[0]*u_l[3]+0.5303300858899105*rhoux_c[1]*u_c[3]-0.3061862178478971*rhoux_c[0]*u_c[3]+0.5303300858899105*u_l[1]*rhoux_l[3]+0.3061862178478971*u_l[0]*rhoux_l[3]+0.4330127018922193*uMax_l[0]*rhoux_l[3]+0.5303300858899105*u_c[1]*rhoux_c[3]-0.3061862178478971*u_c[0]*rhoux_c[3]+0.4330127018922193*uMax_l[0]*rhoux_c[3]+0.3061862178478971*rhoux_l[1]*u_l[2]+0.1767766952966368*rhoux_l[0]*u_l[2]-0.3061862178478971*rhoux_c[1]*u_c[2]+0.1767766952966368*rhoux_c[0]*u_c[2]+0.3061862178478971*u_l[1]*rhoux_l[2]+0.1767766952966368*u_l[0]*rhoux_l[2]+0.25*uMax_l[0]*rhoux_l[2]-0.3061862178478971*u_c[1]*rhoux_c[2]+0.1767766952966368*u_c[0]*rhoux_c[2]-0.25*uMax_l[0]*rhoux_c[2]+0.4330127018922193*rhoux_l[1]*uMax_l[1]+0.4330127018922193*rhoux_c[1]*uMax_l[1]+0.25*rhoux_l[0]*uMax_l[1]-0.25*rhoux_c[0]*uMax_l[1]; 

  Ghat_rhoux_r[0] = 0.5303300858899105*rhoux_r[3]*u_r[3]-0.3061862178478971*rhoux_r[2]*u_r[3]+0.5303300858899105*rhoux_c[3]*u_c[3]+0.3061862178478971*rhoux_c[2]*u_c[3]-0.3061862178478971*u_r[2]*rhoux_r[3]+0.4330127018922193*uMax_r[1]*rhoux_r[3]+0.3061862178478971*u_c[2]*rhoux_c[3]+0.4330127018922193*uMax_r[1]*rhoux_c[3]+0.1767766952966368*rhoux_r[2]*u_r[2]+0.1767766952966368*rhoux_c[2]*u_c[2]-0.25*uMax_r[1]*rhoux_r[2]+0.25*uMax_r[1]*rhoux_c[2]+0.5303300858899105*rhoux_r[1]*u_r[1]-0.3061862178478971*rhoux_r[0]*u_r[1]+0.5303300858899105*rhoux_c[1]*u_c[1]+0.3061862178478971*rhoux_c[0]*u_c[1]-0.3061862178478971*u_r[0]*rhoux_r[1]+0.4330127018922193*uMax_r[0]*rhoux_r[1]+0.3061862178478971*u_c[0]*rhoux_c[1]+0.4330127018922193*uMax_r[0]*rhoux_c[1]+0.1767766952966368*rhoux_r[0]*u_r[0]+0.1767766952966368*rhoux_c[0]*u_c[0]-0.25*rhoux_r[0]*uMax_r[0]+0.25*rhoux_c[0]*uMax_r[0]; 
  Ghat_rhoux_r[1] = 0.5303300858899105*rhoux_r[1]*u_r[3]-0.3061862178478971*rhoux_r[0]*u_r[3]+0.5303300858899105*rhoux_c[1]*u_c[3]+0.3061862178478971*rhoux_c[0]*u_c[3]+0.5303300858899105*u_r[1]*rhoux_r[3]-0.3061862178478971*u_r[0]*rhoux_r[3]+0.4330127018922193*uMax_r[0]*rhoux_r[3]+0.5303300858899105*u_c[1]*rhoux_c[3]+0.3061862178478971*u_c[0]*rhoux_c[3]+0.4330127018922193*uMax_r[0]*rhoux_c[3]-0.3061862178478971*rhoux_r[1]*u_r[2]+0.1767766952966368*rhoux_r[0]*u_r[2]+0.3061862178478971*rhoux_c[1]*u_c[2]+0.1767766952966368*rhoux_c[0]*u_c[2]-0.3061862178478971*u_r[1]*rhoux_r[2]+0.1767766952966368*u_r[0]*rhoux_r[2]-0.25*uMax_r[0]*rhoux_r[2]+0.3061862178478971*u_c[1]*rhoux_c[2]+0.1767766952966368*u_c[0]*rhoux_c[2]+0.25*uMax_r[0]*rhoux_c[2]+0.4330127018922193*rhoux_r[1]*uMax_r[1]+0.4330127018922193*rhoux_c[1]*uMax_r[1]-0.25*rhoux_r[0]*uMax_r[1]+0.25*rhoux_c[0]*uMax_r[1]; 

  Ghat_rhouy_l[0] = 0.5303300858899105*rhouy_l[3]*u_l[3]+0.3061862178478971*rhouy_l[2]*u_l[3]+0.5303300858899105*rhouy_c[3]*u_c[3]-0.3061862178478971*rhouy_c[2]*u_c[3]+0.3061862178478971*u_l[2]*rhouy_l[3]+0.4330127018922193*uMax_l[1]*rhouy_l[3]-0.3061862178478971*u_c[2]*rhouy_c[3]+0.4330127018922193*uMax_l[1]*rhouy_c[3]+0.1767766952966368*rhouy_l[2]*u_l[2]+0.1767766952966368*rhouy_c[2]*u_c[2]+0.25*uMax_l[1]*rhouy_l[2]-0.25*uMax_l[1]*rhouy_c[2]+0.5303300858899105*rhouy_l[1]*u_l[1]+0.3061862178478971*rhouy_l[0]*u_l[1]+0.5303300858899105*rhouy_c[1]*u_c[1]-0.3061862178478971*rhouy_c[0]*u_c[1]+0.3061862178478971*u_l[0]*rhouy_l[1]+0.4330127018922193*uMax_l[0]*rhouy_l[1]-0.3061862178478971*u_c[0]*rhouy_c[1]+0.4330127018922193*uMax_l[0]*rhouy_c[1]+0.1767766952966368*rhouy_l[0]*u_l[0]+0.1767766952966368*rhouy_c[0]*u_c[0]+0.25*rhouy_l[0]*uMax_l[0]-0.25*rhouy_c[0]*uMax_l[0]; 
  Ghat_rhouy_l[1] = 0.5303300858899105*rhouy_l[1]*u_l[3]+0.3061862178478971*rhouy_l[0]*u_l[3]+0.5303300858899105*rhouy_c[1]*u_c[3]-0.3061862178478971*rhouy_c[0]*u_c[3]+0.5303300858899105*u_l[1]*rhouy_l[3]+0.3061862178478971*u_l[0]*rhouy_l[3]+0.4330127018922193*uMax_l[0]*rhouy_l[3]+0.5303300858899105*u_c[1]*rhouy_c[3]-0.3061862178478971*u_c[0]*rhouy_c[3]+0.4330127018922193*uMax_l[0]*rhouy_c[3]+0.3061862178478971*rhouy_l[1]*u_l[2]+0.1767766952966368*rhouy_l[0]*u_l[2]-0.3061862178478971*rhouy_c[1]*u_c[2]+0.1767766952966368*rhouy_c[0]*u_c[2]+0.3061862178478971*u_l[1]*rhouy_l[2]+0.1767766952966368*u_l[0]*rhouy_l[2]+0.25*uMax_l[0]*rhouy_l[2]-0.3061862178478971*u_c[1]*rhouy_c[2]+0.1767766952966368*u_c[0]*rhouy_c[2]-0.25*uMax_l[0]*rhouy_c[2]+0.4330127018922193*rhouy_l[1]*uMax_l[1]+0.4330127018922193*rhouy_c[1]*uMax_l[1]+0.25*rhouy_l[0]*uMax_l[1]-0.25*rhouy_c[0]*uMax_l[1]; 

  Ghat_rhouy_r[0] = 0.5303300858899105*rhouy_r[3]*u_r[3]-0.3061862178478971*rhouy_r[2]*u_r[3]+0.5303300858899105*rhouy_c[3]*u_c[3]+0.3061862178478971*rhouy_c[2]*u_c[3]-0.3061862178478971*u_r[2]*rhouy_r[3]+0.4330127018922193*uMax_r[1]*rhouy_r[3]+0.3061862178478971*u_c[2]*rhouy_c[3]+0.4330127018922193*uMax_r[1]*rhouy_c[3]+0.1767766952966368*rhouy_r[2]*u_r[2]+0.1767766952966368*rhouy_c[2]*u_c[2]-0.25*uMax_r[1]*rhouy_r[2]+0.25*uMax_r[1]*rhouy_c[2]+0.5303300858899105*rhouy_r[1]*u_r[1]-0.3061862178478971*rhouy_r[0]*u_r[1]+0.5303300858899105*rhouy_c[1]*u_c[1]+0.3061862178478971*rhouy_c[0]*u_c[1]-0.3061862178478971*u_r[0]*rhouy_r[1]+0.4330127018922193*uMax_r[0]*rhouy_r[1]+0.3061862178478971*u_c[0]*rhouy_c[1]+0.4330127018922193*uMax_r[0]*rhouy_c[1]+0.1767766952966368*rhouy_r[0]*u_r[0]+0.1767766952966368*rhouy_c[0]*u_c[0]-0.25*rhouy_r[0]*uMax_r[0]+0.25*rhouy_c[0]*uMax_r[0]; 
  Ghat_rhouy_r[1] = 0.5303300858899105*rhouy_r[1]*u_r[3]-0.3061862178478971*rhouy_r[0]*u_r[3]+0.5303300858899105*rhouy_c[1]*u_c[3]+0.3061862178478971*rhouy_c[0]*u_c[3]+0.5303300858899105*u_r[1]*rhouy_r[3]-0.3061862178478971*u_r[0]*rhouy_r[3]+0.4330127018922193*uMax_r[0]*rhouy_r[3]+0.5303300858899105*u_c[1]*rhouy_c[3]+0.3061862178478971*u_c[0]*rhouy_c[3]+0.4330127018922193*uMax_r[0]*rhouy_c[3]-0.3061862178478971*rhouy_r[1]*u_r[2]+0.1767766952966368*rhouy_r[0]*u_r[2]+0.3061862178478971*rhouy_c[1]*u_c[2]+0.1767766952966368*rhouy_c[0]*u_c[2]-0.3061862178478971*u_r[1]*rhouy_r[2]+0.1767766952966368*u_r[0]*rhouy_r[2]-0.25*uMax_r[0]*rhouy_r[2]+0.3061862178478971*u_c[1]*rhouy_c[2]+0.1767766952966368*u_c[0]*rhouy_c[2]+0.25*uMax_r[0]*rhouy_c[2]+0.4330127018922193*rhouy_r[1]*uMax_r[1]+0.4330127018922193*rhouy_c[1]*uMax_r[1]-0.25*rhouy_r[0]*uMax_r[1]+0.25*rhouy_c[0]*uMax_r[1]; 

  Ghat_rhouz_l[0] = 0.5303300858899105*rhouz_l[3]*u_l[3]+0.3061862178478971*rhouz_l[2]*u_l[3]+0.5303300858899105*rhouz_c[3]*u_c[3]-0.3061862178478971*rhouz_c[2]*u_c[3]+0.3061862178478971*u_l[2]*rhouz_l[3]+0.4330127018922193*uMax_l[1]*rhouz_l[3]-0.3061862178478971*u_c[2]*rhouz_c[3]+0.4330127018922193*uMax_l[1]*rhouz_c[3]+0.1767766952966368*rhouz_l[2]*u_l[2]+0.1767766952966368*rhouz_c[2]*u_c[2]+0.25*uMax_l[1]*rhouz_l[2]-0.25*uMax_l[1]*rhouz_c[2]+0.5303300858899105*rhouz_l[1]*u_l[1]+0.3061862178478971*rhouz_l[0]*u_l[1]+0.5303300858899105*rhouz_c[1]*u_c[1]-0.3061862178478971*rhouz_c[0]*u_c[1]+0.3061862178478971*u_l[0]*rhouz_l[1]+0.4330127018922193*uMax_l[0]*rhouz_l[1]-0.3061862178478971*u_c[0]*rhouz_c[1]+0.4330127018922193*uMax_l[0]*rhouz_c[1]+0.1767766952966368*rhouz_l[0]*u_l[0]+0.1767766952966368*rhouz_c[0]*u_c[0]+0.25*rhouz_l[0]*uMax_l[0]-0.25*rhouz_c[0]*uMax_l[0]; 
  Ghat_rhouz_l[1] = 0.5303300858899105*rhouz_l[1]*u_l[3]+0.3061862178478971*rhouz_l[0]*u_l[3]+0.5303300858899105*rhouz_c[1]*u_c[3]-0.3061862178478971*rhouz_c[0]*u_c[3]+0.5303300858899105*u_l[1]*rhouz_l[3]+0.3061862178478971*u_l[0]*rhouz_l[3]+0.4330127018922193*uMax_l[0]*rhouz_l[3]+0.5303300858899105*u_c[1]*rhouz_c[3]-0.3061862178478971*u_c[0]*rhouz_c[3]+0.4330127018922193*uMax_l[0]*rhouz_c[3]+0.3061862178478971*rhouz_l[1]*u_l[2]+0.1767766952966368*rhouz_l[0]*u_l[2]-0.3061862178478971*rhouz_c[1]*u_c[2]+0.1767766952966368*rhouz_c[0]*u_c[2]+0.3061862178478971*u_l[1]*rhouz_l[2]+0.1767766952966368*u_l[0]*rhouz_l[2]+0.25*uMax_l[0]*rhouz_l[2]-0.3061862178478971*u_c[1]*rhouz_c[2]+0.1767766952966368*u_c[0]*rhouz_c[2]-0.25*uMax_l[0]*rhouz_c[2]+0.4330127018922193*rhouz_l[1]*uMax_l[1]+0.4330127018922193*rhouz_c[1]*uMax_l[1]+0.25*rhouz_l[0]*uMax_l[1]-0.25*rhouz_c[0]*uMax_l[1]; 

  Ghat_rhouz_r[0] = 0.5303300858899105*rhouz_r[3]*u_r[3]-0.3061862178478971*rhouz_r[2]*u_r[3]+0.5303300858899105*rhouz_c[3]*u_c[3]+0.3061862178478971*rhouz_c[2]*u_c[3]-0.3061862178478971*u_r[2]*rhouz_r[3]+0.4330127018922193*uMax_r[1]*rhouz_r[3]+0.3061862178478971*u_c[2]*rhouz_c[3]+0.4330127018922193*uMax_r[1]*rhouz_c[3]+0.1767766952966368*rhouz_r[2]*u_r[2]+0.1767766952966368*rhouz_c[2]*u_c[2]-0.25*uMax_r[1]*rhouz_r[2]+0.25*uMax_r[1]*rhouz_c[2]+0.5303300858899105*rhouz_r[1]*u_r[1]-0.3061862178478971*rhouz_r[0]*u_r[1]+0.5303300858899105*rhouz_c[1]*u_c[1]+0.3061862178478971*rhouz_c[0]*u_c[1]-0.3061862178478971*u_r[0]*rhouz_r[1]+0.4330127018922193*uMax_r[0]*rhouz_r[1]+0.3061862178478971*u_c[0]*rhouz_c[1]+0.4330127018922193*uMax_r[0]*rhouz_c[1]+0.1767766952966368*rhouz_r[0]*u_r[0]+0.1767766952966368*rhouz_c[0]*u_c[0]-0.25*rhouz_r[0]*uMax_r[0]+0.25*rhouz_c[0]*uMax_r[0]; 
  Ghat_rhouz_r[1] = 0.5303300858899105*rhouz_r[1]*u_r[3]-0.3061862178478971*rhouz_r[0]*u_r[3]+0.5303300858899105*rhouz_c[1]*u_c[3]+0.3061862178478971*rhouz_c[0]*u_c[3]+0.5303300858899105*u_r[1]*rhouz_r[3]-0.3061862178478971*u_r[0]*rhouz_r[3]+0.4330127018922193*uMax_r[0]*rhouz_r[3]+0.5303300858899105*u_c[1]*rhouz_c[3]+0.3061862178478971*u_c[0]*rhouz_c[3]+0.4330127018922193*uMax_r[0]*rhouz_c[3]-0.3061862178478971*rhouz_r[1]*u_r[2]+0.1767766952966368*rhouz_r[0]*u_r[2]+0.3061862178478971*rhouz_c[1]*u_c[2]+0.1767766952966368*rhouz_c[0]*u_c[2]-0.3061862178478971*u_r[1]*rhouz_r[2]+0.1767766952966368*u_r[0]*rhouz_r[2]-0.25*uMax_r[0]*rhouz_r[2]+0.3061862178478971*u_c[1]*rhouz_c[2]+0.1767766952966368*u_c[0]*rhouz_c[2]+0.25*uMax_r[0]*rhouz_c[2]+0.4330127018922193*rhouz_r[1]*uMax_r[1]+0.4330127018922193*rhouz_c[1]*uMax_r[1]-0.25*rhouz_r[0]*uMax_r[1]+0.25*rhouz_c[0]*uMax_r[1]; 

  Ghat_E_perp_l[0] = 0.5303300858899105*p_perp_l[3]*u_perp_l[3]+0.3061862178478971*p_perp_l[2]*u_perp_l[3]+0.5303300858899105*p_perp_c[3]*u_perp_c[3]-0.3061862178478971*p_perp_c[2]*u_perp_c[3]+0.5303300858899105*E_perp_l[3]*u_l[3]+0.3061862178478971*E_perp_l[2]*u_l[3]+0.5303300858899105*E_perp_c[3]*u_c[3]-0.3061862178478971*E_perp_c[2]*u_c[3]+0.3061862178478971*u_perp_l[2]*p_perp_l[3]+0.4330127018922193*u_perpMax_l[1]*p_perp_l[3]-0.3061862178478971*u_perp_c[2]*p_perp_c[3]+0.4330127018922193*u_perpMax_l[1]*p_perp_c[3]+0.3061862178478971*u_l[2]*E_perp_l[3]+0.4330127018922193*uMax_l[1]*E_perp_l[3]-0.3061862178478971*u_c[2]*E_perp_c[3]+0.4330127018922193*uMax_l[1]*E_perp_c[3]+0.1767766952966368*p_perp_l[2]*u_perp_l[2]+0.1767766952966368*p_perp_c[2]*u_perp_c[2]+0.1767766952966368*E_perp_l[2]*u_l[2]+0.1767766952966368*E_perp_c[2]*u_c[2]+0.25*u_perpMax_l[1]*p_perp_l[2]-0.25*u_perpMax_l[1]*p_perp_c[2]+0.25*uMax_l[1]*E_perp_l[2]-0.25*uMax_l[1]*E_perp_c[2]+0.5303300858899105*p_perp_l[1]*u_perp_l[1]+0.3061862178478971*p_perp_l[0]*u_perp_l[1]+0.5303300858899105*p_perp_c[1]*u_perp_c[1]-0.3061862178478971*p_perp_c[0]*u_perp_c[1]+0.5303300858899105*E_perp_l[1]*u_l[1]+0.3061862178478971*E_perp_l[0]*u_l[1]+0.5303300858899105*E_perp_c[1]*u_c[1]-0.3061862178478971*E_perp_c[0]*u_c[1]+0.3061862178478971*u_perp_l[0]*p_perp_l[1]+0.4330127018922193*u_perpMax_l[0]*p_perp_l[1]-0.3061862178478971*u_perp_c[0]*p_perp_c[1]+0.4330127018922193*u_perpMax_l[0]*p_perp_c[1]+0.3061862178478971*u_l[0]*E_perp_l[1]+0.4330127018922193*uMax_l[0]*E_perp_l[1]-0.3061862178478971*u_c[0]*E_perp_c[1]+0.4330127018922193*uMax_l[0]*E_perp_c[1]+0.1767766952966368*p_perp_l[0]*u_perp_l[0]+0.1767766952966368*p_perp_c[0]*u_perp_c[0]+0.25*p_perp_l[0]*u_perpMax_l[0]-0.25*p_perp_c[0]*u_perpMax_l[0]+0.1767766952966368*E_perp_l[0]*u_l[0]+0.1767766952966368*E_perp_c[0]*u_c[0]+0.25*E_perp_l[0]*uMax_l[0]-0.25*E_perp_c[0]*uMax_l[0]; 
  Ghat_E_perp_l[1] = 0.5303300858899105*p_perp_l[1]*u_perp_l[3]+0.3061862178478971*p_perp_l[0]*u_perp_l[3]+0.5303300858899105*p_perp_c[1]*u_perp_c[3]-0.3061862178478971*p_perp_c[0]*u_perp_c[3]+0.5303300858899105*E_perp_l[1]*u_l[3]+0.3061862178478971*E_perp_l[0]*u_l[3]+0.5303300858899105*E_perp_c[1]*u_c[3]-0.3061862178478971*E_perp_c[0]*u_c[3]+0.5303300858899105*u_perp_l[1]*p_perp_l[3]+0.3061862178478971*u_perp_l[0]*p_perp_l[3]+0.4330127018922193*u_perpMax_l[0]*p_perp_l[3]+0.5303300858899105*u_perp_c[1]*p_perp_c[3]-0.3061862178478971*u_perp_c[0]*p_perp_c[3]+0.4330127018922193*u_perpMax_l[0]*p_perp_c[3]+0.5303300858899105*u_l[1]*E_perp_l[3]+0.3061862178478971*u_l[0]*E_perp_l[3]+0.4330127018922193*uMax_l[0]*E_perp_l[3]+0.5303300858899105*u_c[1]*E_perp_c[3]-0.3061862178478971*u_c[0]*E_perp_c[3]+0.4330127018922193*uMax_l[0]*E_perp_c[3]+0.3061862178478971*p_perp_l[1]*u_perp_l[2]+0.1767766952966368*p_perp_l[0]*u_perp_l[2]-0.3061862178478971*p_perp_c[1]*u_perp_c[2]+0.1767766952966368*p_perp_c[0]*u_perp_c[2]+0.3061862178478971*E_perp_l[1]*u_l[2]+0.1767766952966368*E_perp_l[0]*u_l[2]-0.3061862178478971*E_perp_c[1]*u_c[2]+0.1767766952966368*E_perp_c[0]*u_c[2]+0.3061862178478971*u_perp_l[1]*p_perp_l[2]+0.1767766952966368*u_perp_l[0]*p_perp_l[2]+0.25*u_perpMax_l[0]*p_perp_l[2]-0.3061862178478971*u_perp_c[1]*p_perp_c[2]+0.1767766952966368*u_perp_c[0]*p_perp_c[2]-0.25*u_perpMax_l[0]*p_perp_c[2]+0.3061862178478971*u_l[1]*E_perp_l[2]+0.1767766952966368*u_l[0]*E_perp_l[2]+0.25*uMax_l[0]*E_perp_l[2]-0.3061862178478971*u_c[1]*E_perp_c[2]+0.1767766952966368*u_c[0]*E_perp_c[2]-0.25*uMax_l[0]*E_perp_c[2]+0.4330127018922193*p_perp_l[1]*u_perpMax_l[1]+0.4330127018922193*p_perp_c[1]*u_perpMax_l[1]+0.25*p_perp_l[0]*u_perpMax_l[1]-0.25*p_perp_c[0]*u_perpMax_l[1]+0.4330127018922193*E_perp_l[1]*uMax_l[1]+0.4330127018922193*E_perp_c[1]*uMax_l[1]+0.25*E_perp_l[0]*uMax_l[1]-0.25*E_perp_c[0]*uMax_l[1]; 

  Ghat_E_perp_r[0] = 0.5303300858899105*E_perp_r[3]*u_r[3]-0.3061862178478971*E_perp_r[2]*u_r[3]+0.5303300858899105*p_perp_r[3]*u_perp_r[3]-0.3061862178478971*p_perp_r[2]*u_perp_r[3]+0.5303300858899105*p_perp_c[3]*u_perp_c[3]+0.3061862178478971*p_perp_c[2]*u_perp_c[3]+0.5303300858899105*E_perp_c[3]*u_c[3]+0.3061862178478971*E_perp_c[2]*u_c[3]-0.3061862178478971*u_perp_r[2]*p_perp_r[3]+0.4330127018922193*u_perpMax_r[1]*p_perp_r[3]+0.3061862178478971*u_perp_c[2]*p_perp_c[3]+0.4330127018922193*u_perpMax_r[1]*p_perp_c[3]-0.3061862178478971*u_r[2]*E_perp_r[3]+0.4330127018922193*uMax_r[1]*E_perp_r[3]+0.3061862178478971*u_c[2]*E_perp_c[3]+0.4330127018922193*uMax_r[1]*E_perp_c[3]+0.1767766952966368*E_perp_r[2]*u_r[2]+0.1767766952966368*p_perp_r[2]*u_perp_r[2]+0.1767766952966368*p_perp_c[2]*u_perp_c[2]+0.1767766952966368*E_perp_c[2]*u_c[2]-0.25*u_perpMax_r[1]*p_perp_r[2]+0.25*u_perpMax_r[1]*p_perp_c[2]-0.25*uMax_r[1]*E_perp_r[2]+0.25*uMax_r[1]*E_perp_c[2]+0.5303300858899105*E_perp_r[1]*u_r[1]-0.3061862178478971*E_perp_r[0]*u_r[1]+0.5303300858899105*p_perp_r[1]*u_perp_r[1]-0.3061862178478971*p_perp_r[0]*u_perp_r[1]+0.5303300858899105*p_perp_c[1]*u_perp_c[1]+0.3061862178478971*p_perp_c[0]*u_perp_c[1]+0.5303300858899105*E_perp_c[1]*u_c[1]+0.3061862178478971*E_perp_c[0]*u_c[1]-0.3061862178478971*u_perp_r[0]*p_perp_r[1]+0.4330127018922193*u_perpMax_r[0]*p_perp_r[1]+0.3061862178478971*u_perp_c[0]*p_perp_c[1]+0.4330127018922193*u_perpMax_r[0]*p_perp_c[1]-0.3061862178478971*u_r[0]*E_perp_r[1]+0.4330127018922193*uMax_r[0]*E_perp_r[1]+0.3061862178478971*u_c[0]*E_perp_c[1]+0.4330127018922193*uMax_r[0]*E_perp_c[1]+0.1767766952966368*E_perp_r[0]*u_r[0]+0.1767766952966368*p_perp_r[0]*u_perp_r[0]+0.1767766952966368*p_perp_c[0]*u_perp_c[0]-0.25*p_perp_r[0]*u_perpMax_r[0]+0.25*p_perp_c[0]*u_perpMax_r[0]+0.1767766952966368*E_perp_c[0]*u_c[0]-0.25*E_perp_r[0]*uMax_r[0]+0.25*E_perp_c[0]*uMax_r[0]; 
  Ghat_E_perp_r[1] = 0.5303300858899105*E_perp_r[1]*u_r[3]-0.3061862178478971*E_perp_r[0]*u_r[3]+0.5303300858899105*p_perp_r[1]*u_perp_r[3]-0.3061862178478971*p_perp_r[0]*u_perp_r[3]+0.5303300858899105*p_perp_c[1]*u_perp_c[3]+0.3061862178478971*p_perp_c[0]*u_perp_c[3]+0.5303300858899105*E_perp_c[1]*u_c[3]+0.3061862178478971*E_perp_c[0]*u_c[3]+0.5303300858899105*u_perp_r[1]*p_perp_r[3]-0.3061862178478971*u_perp_r[0]*p_perp_r[3]+0.4330127018922193*u_perpMax_r[0]*p_perp_r[3]+0.5303300858899105*u_perp_c[1]*p_perp_c[3]+0.3061862178478971*u_perp_c[0]*p_perp_c[3]+0.4330127018922193*u_perpMax_r[0]*p_perp_c[3]+0.5303300858899105*u_r[1]*E_perp_r[3]-0.3061862178478971*u_r[0]*E_perp_r[3]+0.4330127018922193*uMax_r[0]*E_perp_r[3]+0.5303300858899105*u_c[1]*E_perp_c[3]+0.3061862178478971*u_c[0]*E_perp_c[3]+0.4330127018922193*uMax_r[0]*E_perp_c[3]-0.3061862178478971*E_perp_r[1]*u_r[2]+0.1767766952966368*E_perp_r[0]*u_r[2]-0.3061862178478971*p_perp_r[1]*u_perp_r[2]+0.1767766952966368*p_perp_r[0]*u_perp_r[2]+0.3061862178478971*p_perp_c[1]*u_perp_c[2]+0.1767766952966368*p_perp_c[0]*u_perp_c[2]+0.3061862178478971*E_perp_c[1]*u_c[2]+0.1767766952966368*E_perp_c[0]*u_c[2]-0.3061862178478971*u_perp_r[1]*p_perp_r[2]+0.1767766952966368*u_perp_r[0]*p_perp_r[2]-0.25*u_perpMax_r[0]*p_perp_r[2]+0.3061862178478971*u_perp_c[1]*p_perp_c[2]+0.1767766952966368*u_perp_c[0]*p_perp_c[2]+0.25*u_perpMax_r[0]*p_perp_c[2]-0.3061862178478971*u_r[1]*E_perp_r[2]+0.1767766952966368*u_r[0]*E_perp_r[2]-0.25*uMax_r[0]*E_perp_r[2]+0.3061862178478971*u_c[1]*E_perp_c[2]+0.1767766952966368*u_c[0]*E_perp_c[2]+0.25*uMax_r[0]*E_perp_c[2]+0.4330127018922193*p_perp_r[1]*u_perpMax_r[1]+0.4330127018922193*p_perp_c[1]*u_perpMax_r[1]-0.25*p_perp_r[0]*u_perpMax_r[1]+0.25*p_perp_c[0]*u_perpMax_r[1]+0.4330127018922193*E_perp_r[1]*uMax_r[1]+0.4330127018922193*E_perp_c[1]*uMax_r[1]-0.25*E_perp_r[0]*uMax_r[1]+0.25*E_perp_c[0]*uMax_r[1]; 

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

  outE_perp[0] += 0.7071067811865475*Ghat_E_perp_l[0]*dx1-0.7071067811865475*Ghat_E_perp_r[0]*dx1; 
  outE_perp[1] += (-1.224744871391589*Ghat_E_perp_r[0]*dx1)-1.224744871391589*Ghat_E_perp_l[0]*dx1; 
  outE_perp[2] += 0.7071067811865475*Ghat_E_perp_l[1]*dx1-0.7071067811865475*Ghat_E_perp_r[1]*dx1; 
  outE_perp[3] += (-1.224744871391589*Ghat_E_perp_r[1]*dx1)-1.224744871391589*Ghat_E_perp_l[1]*dx1; 

} 
