#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void euler_pkpm_surfy_2x_ser_p1(const double *w, const double *dxv,
  const double *u_il, const double *u_ic, const double *u_ir,
  const double *T_ijl, const double *T_ijc, const double *T_ijr,
  const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                       Cell-center coordinates.
  // dxv[NDIM]:                     Cell spacing.
  // u_il/u_ic/u_ir:                Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // T_ijl/T_ijc/T_ijr:             Input Temperature tensor/mass (for penalization) in left/center/right cells.
  // statevecl/statevecc/statevecr: [rho ux, rho uy, rho uz, E_perp], Fluid input state vector in left/center/right cells.
  // out: Incremented output.

  const double dx1 = 2.0/dxv[1]; 

  const double *rhoux_l = &statevecl[0]; 
  const double *rhouy_l = &statevecl[4]; 
  const double *rhouz_l = &statevecl[8]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[4]; 
  const double *rhouz_c = &statevecc[8]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[4]; 
  const double *rhouz_r = &statevecr[8]; 

  const double *u_l = &u_il[4]; 
  const double *u_c = &u_ic[4]; 
  const double *u_r = &u_ir[4]; 

  // Get thermal velocity in direction of update for penalization vth^2 = 3.0*T_ii/m. 
  const double *vth_sql = &T_ijl[12]; 
  const double *vth_sqc = &T_ijc[12]; 
  const double *vth_sqr = &T_ijr[12]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[4]; 
  double *outrhou2 = &out[8]; 

  double lax_rhoux_quad_l[2] = {0.0};
  double lax_rhoux_quad_r[2] = {0.0};
  double lax_rhoux_modal_l[2] = {0.0};
  double lax_rhoux_modal_r[2] = {0.0};
  double lax_rhouy_quad_l[2] = {0.0};
  double lax_rhouy_quad_r[2] = {0.0};
  double lax_rhouy_modal_l[2] = {0.0};
  double lax_rhouy_modal_r[2] = {0.0};
  double lax_rhouz_quad_l[2] = {0.0};
  double lax_rhouz_quad_r[2] = {0.0};
  double lax_rhouz_modal_l[2] = {0.0};
  double lax_rhouz_modal_r[2] = {0.0};

  double u_l_r = 0.0; 
  double u_c_l = 0.0; 
  double u_c_r = 0.0; 
  double u_r_l = 0.0; 
  double uQuad_l = 0.0; 
  double uQuad_r = 0.0; 
  double vth_sq_l_r = 0.0; 
  double vth_sq_c_l = 0.0; 
  double vth_sq_c_r = 0.0; 
  double vth_sq_r_l = 0.0; 
  double vthQuad_l = 0.0; 
  double vthQuad_r = 0.0; 
  double max_speed_l = 0.0; 
  double max_speed_r = 0.0; 
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

  u_l_r = ser_2x_p1_surfx2_eval_quad_node_0_r(u_l); 
  u_c_l = ser_2x_p1_surfx2_eval_quad_node_0_l(u_c); 
  u_c_r = ser_2x_p1_surfx2_eval_quad_node_0_r(u_c); 
  u_r_l = ser_2x_p1_surfx2_eval_quad_node_0_l(u_r); 
  uQuad_l = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uQuad_r = fmax(fabs(u_c_r), fabs(u_r_l)); 
  vth_sq_l_r = ser_2x_p1_surfx2_eval_quad_node_0_r(vth_sql); 
  vth_sq_c_l = ser_2x_p1_surfx2_eval_quad_node_0_l(vth_sqc); 
  vth_sq_c_r = ser_2x_p1_surfx2_eval_quad_node_0_r(vth_sqc); 
  vth_sq_r_l = ser_2x_p1_surfx2_eval_quad_node_0_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = uQuad_l + vthQuad_l; 
  max_speed_r = uQuad_r + vthQuad_r; 
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
  lax_rhoux_quad_l[0] = 0.5*(rhoux_l_r*u_l_r + rhoux_c_l*u_c_l) - 0.5*max_speed_l*(rhoux_c_l - rhoux_l_r); 
  lax_rhoux_quad_r[0] = 0.5*(rhoux_c_r*u_c_r + rhoux_r_l*u_r_l) - 0.5*max_speed_r*(rhoux_r_l - rhoux_c_r); 
  lax_rhouy_quad_l[0] = 0.5*(rhouy_l_r*u_l_r + rhouy_c_l*u_c_l) - 0.5*max_speed_l*(rhouy_c_l - rhouy_l_r); 
  lax_rhouy_quad_r[0] = 0.5*(rhouy_c_r*u_c_r + rhouy_r_l*u_r_l) - 0.5*max_speed_r*(rhouy_r_l - rhouy_c_r); 
  lax_rhouz_quad_l[0] = 0.5*(rhouz_l_r*u_l_r + rhouz_c_l*u_c_l) - 0.5*max_speed_l*(rhouz_c_l - rhouz_l_r); 
  lax_rhouz_quad_r[0] = 0.5*(rhouz_c_r*u_c_r + rhouz_r_l*u_r_l) - 0.5*max_speed_r*(rhouz_r_l - rhouz_c_r); 

  u_l_r = ser_2x_p1_surfx2_eval_quad_node_1_r(u_l); 
  u_c_l = ser_2x_p1_surfx2_eval_quad_node_1_l(u_c); 
  u_c_r = ser_2x_p1_surfx2_eval_quad_node_1_r(u_c); 
  u_r_l = ser_2x_p1_surfx2_eval_quad_node_1_l(u_r); 
  uQuad_l = fmax(fabs(u_l_r), fabs(u_c_l)); 
  uQuad_r = fmax(fabs(u_c_r), fabs(u_r_l)); 
  vth_sq_l_r = ser_2x_p1_surfx2_eval_quad_node_1_r(vth_sql); 
  vth_sq_c_l = ser_2x_p1_surfx2_eval_quad_node_1_l(vth_sqc); 
  vth_sq_c_r = ser_2x_p1_surfx2_eval_quad_node_1_r(vth_sqc); 
  vth_sq_r_l = ser_2x_p1_surfx2_eval_quad_node_1_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = uQuad_l + vthQuad_l; 
  max_speed_r = uQuad_r + vthQuad_r; 
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
  lax_rhoux_quad_l[1] = 0.5*(rhoux_l_r*u_l_r + rhoux_c_l*u_c_l) - 0.5*max_speed_l*(rhoux_c_l - rhoux_l_r); 
  lax_rhoux_quad_r[1] = 0.5*(rhoux_c_r*u_c_r + rhoux_r_l*u_r_l) - 0.5*max_speed_r*(rhoux_r_l - rhoux_c_r); 
  lax_rhouy_quad_l[1] = 0.5*(rhouy_l_r*u_l_r + rhouy_c_l*u_c_l) - 0.5*max_speed_l*(rhouy_c_l - rhouy_l_r); 
  lax_rhouy_quad_r[1] = 0.5*(rhouy_c_r*u_c_r + rhouy_r_l*u_r_l) - 0.5*max_speed_r*(rhouy_r_l - rhouy_c_r); 
  lax_rhouz_quad_l[1] = 0.5*(rhouz_l_r*u_l_r + rhouz_c_l*u_c_l) - 0.5*max_speed_l*(rhouz_c_l - rhouz_l_r); 
  lax_rhouz_quad_r[1] = 0.5*(rhouz_c_r*u_c_r + rhouz_r_l*u_r_l) - 0.5*max_speed_r*(rhouz_r_l - rhouz_c_r); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p1_upwind_quad_to_modal(lax_rhoux_quad_l, lax_rhoux_modal_l); 
  ser_2x_p1_upwind_quad_to_modal(lax_rhoux_quad_r, lax_rhoux_modal_r); 
  ser_2x_p1_upwind_quad_to_modal(lax_rhouy_quad_l, lax_rhouy_modal_l); 
  ser_2x_p1_upwind_quad_to_modal(lax_rhouy_quad_r, lax_rhouy_modal_r); 
  ser_2x_p1_upwind_quad_to_modal(lax_rhouz_quad_l, lax_rhouz_modal_l); 
  ser_2x_p1_upwind_quad_to_modal(lax_rhouz_quad_r, lax_rhouz_modal_r); 

  outrhou0[0] += 0.7071067811865475*lax_rhoux_modal_l[0]*dx1-0.7071067811865475*lax_rhoux_modal_r[0]*dx1; 
  outrhou0[1] += 0.7071067811865475*lax_rhoux_modal_l[1]*dx1-0.7071067811865475*lax_rhoux_modal_r[1]*dx1; 
  outrhou0[2] += (-1.224744871391589*lax_rhoux_modal_r[0]*dx1)-1.224744871391589*lax_rhoux_modal_l[0]*dx1; 
  outrhou0[3] += (-1.224744871391589*lax_rhoux_modal_r[1]*dx1)-1.224744871391589*lax_rhoux_modal_l[1]*dx1; 

  outrhou1[0] += 0.7071067811865475*lax_rhouy_modal_l[0]*dx1-0.7071067811865475*lax_rhouy_modal_r[0]*dx1; 
  outrhou1[1] += 0.7071067811865475*lax_rhouy_modal_l[1]*dx1-0.7071067811865475*lax_rhouy_modal_r[1]*dx1; 
  outrhou1[2] += (-1.224744871391589*lax_rhouy_modal_r[0]*dx1)-1.224744871391589*lax_rhouy_modal_l[0]*dx1; 
  outrhou1[3] += (-1.224744871391589*lax_rhouy_modal_r[1]*dx1)-1.224744871391589*lax_rhouy_modal_l[1]*dx1; 

  outrhou2[0] += 0.7071067811865475*lax_rhouz_modal_l[0]*dx1-0.7071067811865475*lax_rhouz_modal_r[0]*dx1; 
  outrhou2[1] += 0.7071067811865475*lax_rhouz_modal_l[1]*dx1-0.7071067811865475*lax_rhouz_modal_r[1]*dx1; 
  outrhou2[2] += (-1.224744871391589*lax_rhouz_modal_r[0]*dx1)-1.224744871391589*lax_rhouz_modal_l[0]*dx1; 
  outrhou2[3] += (-1.224744871391589*lax_rhouz_modal_r[1]*dx1)-1.224744871391589*lax_rhouz_modal_l[1]*dx1; 

} 
