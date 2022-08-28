#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p2_surfx1_eval_quad.h> 
GKYL_CU_DH void euler_surfx_1x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gas_gamma: Adiabatic index.
  // ul/uc/ur: [ux, uy, uz] Fluid flow in left/center/right cells.
  // pl/pc/pr: Fluid pressure in left/center/right cells.
  // statevecl/statevecc/statevecr: [rho, rho ux, rho uy, rho uz, energy], Fluid input state vector in left/center/right cells.
  // out: Incremented output.
  const double dx1 = 2.0/dxv[0]; 
  const double *rho_l = &statevecl[0]; 
  const double *rhou0_l = &statevecl[3]; 
  const double *rhou1_l = &statevecl[6]; 
  const double *rhou2_l = &statevecl[9]; 
  const double *energy_l = &statevecl[12]; 

  const double *rho_c = &statevecc[0]; 
  const double *rhou0_c = &statevecc[3]; 
  const double *rhou1_c = &statevecc[6]; 
  const double *rhou2_c = &statevecc[9]; 
  const double *energy_c = &statevecc[12]; 

  const double *rho_r = &statevecr[0]; 
  const double *rhou0_r = &statevecr[3]; 
  const double *rhou1_r = &statevecr[6]; 
  const double *rhou2_r = &statevecr[9]; 
  const double *energy_r = &statevecr[12]; 

  const double *ul_0 = &ul[0]; 
  const double *uc_0 = &uc[0]; 
  const double *ur_0 = &ur[0]; 

  const double *ul_1 = &ul[3]; 
  const double *uc_1 = &uc[3]; 
  const double *ur_1 = &ur[3]; 

  const double *ul_2 = &ul[6]; 
  const double *uc_2 = &uc[6]; 
  const double *ur_2 = &ur[6]; 

  double *outrho = &out[0]; 
  double *outrhou0 = &out[3]; 
  double *outrhou1 = &out[6]; 
  double *outrhou2 = &out[9]; 
  double *outenergy = &out[12]; 

  double u_l_r = ser_1x_p2_surfx1_eval_quad_node_0_r(ul_0); 
  double u_c_l = ser_1x_p2_surfx1_eval_quad_node_0_l(uc_0); 
  double u_c_r = ser_1x_p2_surfx1_eval_quad_node_0_r(uc_0); 
  double u_r_l = ser_1x_p2_surfx1_eval_quad_node_0_l(ur_0); 
  double u_max_l = fmax(fabs(u_l_r), fabs(u_c_l)); 
  double u_max_r = fmax(fabs(u_c_r), fabs(u_r_l)); 

  double Ghat_rho_l = 0.7905694150420948*rho_l[2]*u_max_l-0.7905694150420948*rho_c[2]*u_max_l+0.6123724356957945*rho_l[1]*u_max_l+0.6123724356957945*rho_c[1]*u_max_l+0.3535533905932737*rho_l[0]*u_max_l-0.3535533905932737*rho_c[0]*u_max_l+1.25*rho_l[2]*ul_0[2]+0.9682458365518543*rho_l[1]*ul_0[2]+0.5590169943749475*rho_l[0]*ul_0[2]+1.25*rho_c[2]*uc_0[2]-0.9682458365518543*rho_c[1]*uc_0[2]+0.5590169943749475*rho_c[0]*uc_0[2]+0.9682458365518543*ul_0[1]*rho_l[2]+0.5590169943749475*ul_0[0]*rho_l[2]-0.9682458365518543*uc_0[1]*rho_c[2]+0.5590169943749475*uc_0[0]*rho_c[2]+0.75*rho_l[1]*ul_0[1]+0.4330127018922193*rho_l[0]*ul_0[1]+0.75*rho_c[1]*uc_0[1]-0.4330127018922193*rho_c[0]*uc_0[1]+0.4330127018922193*ul_0[0]*rho_l[1]-0.4330127018922193*uc_0[0]*rho_c[1]+0.25*rho_l[0]*ul_0[0]+0.25*rho_c[0]*uc_0[0]; 
  double Ghat_rho_r = (-0.7905694150420948*rho_r[2]*u_max_r)+0.7905694150420948*rho_c[2]*u_max_r+0.6123724356957945*rho_r[1]*u_max_r+0.6123724356957945*rho_c[1]*u_max_r-0.3535533905932737*rho_r[0]*u_max_r+0.3535533905932737*rho_c[0]*u_max_r+1.25*rho_r[2]*ur_0[2]-0.9682458365518543*rho_r[1]*ur_0[2]+0.5590169943749475*rho_r[0]*ur_0[2]+1.25*rho_c[2]*uc_0[2]+0.9682458365518543*rho_c[1]*uc_0[2]+0.5590169943749475*rho_c[0]*uc_0[2]-0.9682458365518543*ur_0[1]*rho_r[2]+0.5590169943749475*ur_0[0]*rho_r[2]+0.9682458365518543*uc_0[1]*rho_c[2]+0.5590169943749475*uc_0[0]*rho_c[2]+0.75*rho_r[1]*ur_0[1]-0.4330127018922193*rho_r[0]*ur_0[1]+0.75*rho_c[1]*uc_0[1]+0.4330127018922193*rho_c[0]*uc_0[1]-0.4330127018922193*ur_0[0]*rho_r[1]+0.4330127018922193*uc_0[0]*rho_c[1]+0.25*rho_r[0]*ur_0[0]+0.25*rho_c[0]*uc_0[0]; 
  double Ghat_energy_l = (-0.7905694150420948*pl[2]*u_max_l)+0.7905694150420948*pc[2]*u_max_l-0.7905694150420948*energy_l[2]*u_max_l+0.7905694150420948*energy_c[2]*u_max_l-0.6123724356957945*pl[1]*u_max_l-0.6123724356957945*pc[1]*u_max_l-0.6123724356957945*energy_l[1]*u_max_l-0.6123724356957945*energy_c[1]*u_max_l-0.3535533905932737*pl[0]*u_max_l+0.3535533905932737*pc[0]*u_max_l-0.3535533905932737*energy_l[0]*u_max_l+0.3535533905932737*energy_c[0]*u_max_l+1.25*pl[2]*ul_0[2]+1.25*energy_l[2]*ul_0[2]+0.9682458365518543*pl[1]*ul_0[2]+0.9682458365518543*energy_l[1]*ul_0[2]+0.5590169943749475*pl[0]*ul_0[2]+0.5590169943749475*energy_l[0]*ul_0[2]+1.25*pc[2]*uc_0[2]+1.25*energy_c[2]*uc_0[2]-0.9682458365518543*pc[1]*uc_0[2]-0.9682458365518543*energy_c[1]*uc_0[2]+0.5590169943749475*pc[0]*uc_0[2]+0.5590169943749475*energy_c[0]*uc_0[2]+0.9682458365518543*ul_0[1]*pl[2]+0.5590169943749475*ul_0[0]*pl[2]-0.9682458365518543*uc_0[1]*pc[2]+0.5590169943749475*uc_0[0]*pc[2]+0.9682458365518543*ul_0[1]*energy_l[2]+0.5590169943749475*ul_0[0]*energy_l[2]-0.9682458365518543*uc_0[1]*energy_c[2]+0.5590169943749475*uc_0[0]*energy_c[2]+0.75*pl[1]*ul_0[1]+0.75*energy_l[1]*ul_0[1]+0.4330127018922193*pl[0]*ul_0[1]+0.4330127018922193*energy_l[0]*ul_0[1]+0.75*pc[1]*uc_0[1]+0.75*energy_c[1]*uc_0[1]-0.4330127018922193*pc[0]*uc_0[1]-0.4330127018922193*energy_c[0]*uc_0[1]+0.4330127018922193*ul_0[0]*pl[1]-0.4330127018922193*uc_0[0]*pc[1]+0.4330127018922193*ul_0[0]*energy_l[1]-0.4330127018922193*uc_0[0]*energy_c[1]+0.25*pl[0]*ul_0[0]+0.25*energy_l[0]*ul_0[0]+0.25*pc[0]*uc_0[0]+0.25*energy_c[0]*uc_0[0]; 
  double Ghat_energy_r = 0.7905694150420948*pr[2]*u_max_r-0.7905694150420948*pc[2]*u_max_r+0.7905694150420948*energy_r[2]*u_max_r-0.7905694150420948*energy_c[2]*u_max_r-0.6123724356957945*pr[1]*u_max_r-0.6123724356957945*pc[1]*u_max_r-0.6123724356957945*energy_r[1]*u_max_r-0.6123724356957945*energy_c[1]*u_max_r+0.3535533905932737*pr[0]*u_max_r-0.3535533905932737*pc[0]*u_max_r+0.3535533905932737*energy_r[0]*u_max_r-0.3535533905932737*energy_c[0]*u_max_r+1.25*pr[2]*ur_0[2]+1.25*energy_r[2]*ur_0[2]-0.9682458365518543*pr[1]*ur_0[2]-0.9682458365518543*energy_r[1]*ur_0[2]+0.5590169943749475*pr[0]*ur_0[2]+0.5590169943749475*energy_r[0]*ur_0[2]+1.25*pc[2]*uc_0[2]+1.25*energy_c[2]*uc_0[2]+0.9682458365518543*pc[1]*uc_0[2]+0.9682458365518543*energy_c[1]*uc_0[2]+0.5590169943749475*pc[0]*uc_0[2]+0.5590169943749475*energy_c[0]*uc_0[2]-0.9682458365518543*ur_0[1]*pr[2]+0.5590169943749475*ur_0[0]*pr[2]+0.9682458365518543*uc_0[1]*pc[2]+0.5590169943749475*uc_0[0]*pc[2]-0.9682458365518543*ur_0[1]*energy_r[2]+0.5590169943749475*ur_0[0]*energy_r[2]+0.9682458365518543*uc_0[1]*energy_c[2]+0.5590169943749475*uc_0[0]*energy_c[2]+0.75*pr[1]*ur_0[1]+0.75*energy_r[1]*ur_0[1]-0.4330127018922193*pr[0]*ur_0[1]-0.4330127018922193*energy_r[0]*ur_0[1]+0.75*pc[1]*uc_0[1]+0.75*energy_c[1]*uc_0[1]+0.4330127018922193*pc[0]*uc_0[1]+0.4330127018922193*energy_c[0]*uc_0[1]-0.4330127018922193*ur_0[0]*pr[1]+0.4330127018922193*uc_0[0]*pc[1]-0.4330127018922193*ur_0[0]*energy_r[1]+0.4330127018922193*uc_0[0]*energy_c[1]+0.25*pr[0]*ur_0[0]+0.25*energy_r[0]*ur_0[0]+0.25*pc[0]*uc_0[0]+0.25*energy_c[0]*uc_0[0]; 

  double urec_0_l = 0.3458741190809163*ul_0[2]+0.3458741190809163*uc_0[2]+0.4975526040028326*ul_0[1]-0.4975526040028326*uc_0[1]+0.3535533905932737*ul_0[0]+0.3535533905932737*uc_0[0]; 
  double urec_0_r = 0.3458741190809163*ur_0[2]+0.3458741190809163*uc_0[2]-0.4975526040028326*ur_0[1]+0.4975526040028326*uc_0[1]+0.3535533905932737*ur_0[0]+0.3535533905932737*uc_0[0]; 
  double urec_1_l = 0.3458741190809163*ul_1[2]+0.3458741190809163*uc_1[2]+0.4975526040028326*ul_1[1]-0.4975526040028326*uc_1[1]+0.3535533905932737*ul_1[0]+0.3535533905932737*uc_1[0]; 
  double urec_1_r = 0.3458741190809163*ur_1[2]+0.3458741190809163*uc_1[2]-0.4975526040028326*ur_1[1]+0.4975526040028326*uc_1[1]+0.3535533905932737*ur_1[0]+0.3535533905932737*uc_1[0]; 
  double urec_2_l = 0.3458741190809163*ul_2[2]+0.3458741190809163*uc_2[2]+0.4975526040028326*ul_2[1]-0.4975526040028326*uc_2[1]+0.3535533905932737*ul_2[0]+0.3535533905932737*uc_2[0]; 
  double urec_2_r = 0.3458741190809163*ur_2[2]+0.3458741190809163*uc_2[2]-0.4975526040028326*ur_2[1]+0.4975526040028326*uc_2[1]+0.3535533905932737*ur_2[0]+0.3535533905932737*uc_2[0]; 
  double prec_l = 0.3458741190809163*pl[2]+0.3458741190809163*pc[2]+0.4975526040028326*pl[1]-0.4975526040028326*pc[1]+0.3535533905932737*pl[0]+0.3535533905932737*pc[0]; 
  double prec_r = 0.3458741190809163*pr[2]+0.3458741190809163*pc[2]-0.4975526040028326*pr[1]+0.4975526040028326*pc[1]+0.3535533905932737*pr[0]+0.3535533905932737*pc[0]; 

  outrho[0] += 0.7071067811865475*Ghat_rho_l*dx1-0.7071067811865475*Ghat_rho_r*dx1; 
  outrho[1] += (-1.224744871391589*Ghat_rho_r*dx1)-1.224744871391589*Ghat_rho_l*dx1; 
  outrho[2] += 1.58113883008419*Ghat_rho_l*dx1-1.58113883008419*Ghat_rho_r*dx1; 

  outrhou0[0] += (-0.7071067811865475*Ghat_rho_r*dx1*urec_0_r)+0.7071067811865475*Ghat_rho_l*dx1*urec_0_l-0.7071067811865475*dx1*prec_r+0.7071067811865475*dx1*prec_l; 
  outrhou0[1] += (-1.224744871391589*Ghat_rho_r*dx1*urec_0_r)-1.224744871391589*Ghat_rho_l*dx1*urec_0_l-1.224744871391589*dx1*prec_r-1.224744871391589*dx1*prec_l; 
  outrhou0[2] += (-1.58113883008419*Ghat_rho_r*dx1*urec_0_r)+1.58113883008419*Ghat_rho_l*dx1*urec_0_l-1.58113883008419*dx1*prec_r+1.58113883008419*dx1*prec_l; 

  outrhou1[0] += 0.7071067811865475*Ghat_rho_l*dx1*urec_1_l-0.7071067811865475*Ghat_rho_r*dx1*urec_1_r; 
  outrhou1[1] += (-1.224744871391589*Ghat_rho_r*dx1*urec_1_r)-1.224744871391589*Ghat_rho_l*dx1*urec_1_l; 
  outrhou1[2] += 1.58113883008419*Ghat_rho_l*dx1*urec_1_l-1.58113883008419*Ghat_rho_r*dx1*urec_1_r; 

  outrhou2[0] += 0.7071067811865475*Ghat_rho_l*dx1*urec_2_l-0.7071067811865475*Ghat_rho_r*dx1*urec_2_r; 
  outrhou2[1] += (-1.224744871391589*Ghat_rho_r*dx1*urec_2_r)-1.224744871391589*Ghat_rho_l*dx1*urec_2_l; 
  outrhou2[2] += 1.58113883008419*Ghat_rho_l*dx1*urec_2_l-1.58113883008419*Ghat_rho_r*dx1*urec_2_r; 

  outenergy[0] += 0.7071067811865475*Ghat_energy_l*dx1-0.7071067811865475*Ghat_energy_r*dx1; 
  outenergy[1] += (-1.224744871391589*Ghat_energy_r*dx1)-1.224744871391589*Ghat_energy_l*dx1; 
  outenergy[2] += 1.58113883008419*Ghat_energy_l*dx1-1.58113883008419*Ghat_energy_r*dx1; 

} 
