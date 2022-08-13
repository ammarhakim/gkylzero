#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p1_surfx1_eval_quad.h> 
GKYL_CU_DH void euler_surfx_1x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
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
  const double *rhou0_l = &statevecl[2]; 
  const double *rhou1_l = &statevecl[4]; 
  const double *rhou2_l = &statevecl[6]; 
  const double *energy_l = &statevecl[8]; 

  const double *rho_c = &statevecc[0]; 
  const double *rhou0_c = &statevecc[2]; 
  const double *rhou1_c = &statevecc[4]; 
  const double *rhou2_c = &statevecc[6]; 
  const double *energy_c = &statevecc[8]; 

  const double *rho_r = &statevecr[0]; 
  const double *rhou0_r = &statevecr[2]; 
  const double *rhou1_r = &statevecr[4]; 
  const double *rhou2_r = &statevecr[6]; 
  const double *energy_r = &statevecr[8]; 

  const double *ul_0 = &ul[0]; 
  const double *uc_0 = &uc[0]; 
  const double *ur_0 = &ur[0]; 

  const double *ul_1 = &ul[2]; 
  const double *uc_1 = &uc[2]; 
  const double *ur_1 = &ur[2]; 

  const double *ul_2 = &ul[4]; 
  const double *uc_2 = &uc[4]; 
  const double *ur_2 = &ur[4]; 

  double *outrho = &out[0]; 
  double *outrhou0 = &out[2]; 
  double *outrhou1 = &out[4]; 
  double *outrhou2 = &out[6]; 
  double *outenergy = &out[8]; 

  double u_l_r = ser_1x_p1_surfx1_eval_quad_node_0_r(ul_0); 
  double u_c_l = ser_1x_p1_surfx1_eval_quad_node_0_l(uc_0); 
  double u_c_r = ser_1x_p1_surfx1_eval_quad_node_0_r(uc_0); 
  double u_r_l = ser_1x_p1_surfx1_eval_quad_node_0_l(ur_0); 
  double u_max_l = fabs(fmax(u_l_r, u_c_l)); 
  double u_max_r= fabs(fmax(u_c_r, u_r_l)); 

  double Ghat_rho_l = 0.6123724356957945*rho_l[1]*u_max_l+0.6123724356957945*rho_c[1]*u_max_l+0.3535533905932737*rho_l[0]*u_max_l-0.3535533905932737*rho_c[0]*u_max_l+0.75*rho_l[1]*ul_0[1]+0.4330127018922193*rho_l[0]*ul_0[1]+0.75*rho_c[1]*uc_0[1]-0.4330127018922193*rho_c[0]*uc_0[1]+0.4330127018922193*ul_0[0]*rho_l[1]-0.4330127018922193*uc_0[0]*rho_c[1]+0.25*rho_l[0]*ul_0[0]+0.25*rho_c[0]*uc_0[0]; 
  double Ghat_rho_r = 0.6123724356957945*rho_r[1]*u_max_r+0.6123724356957945*rho_c[1]*u_max_r-0.3535533905932737*rho_r[0]*u_max_r+0.3535533905932737*rho_c[0]*u_max_r+0.75*rho_r[1]*ur_0[1]-0.4330127018922193*rho_r[0]*ur_0[1]+0.75*rho_c[1]*uc_0[1]+0.4330127018922193*rho_c[0]*uc_0[1]-0.4330127018922193*ur_0[0]*rho_r[1]+0.4330127018922193*uc_0[0]*rho_c[1]+0.25*rho_r[0]*ur_0[0]+0.25*rho_c[0]*uc_0[0]; 
  double Ghat_energy_l = (-0.6123724356957945*pl[1]*u_max_l)-0.6123724356957945*pc[1]*u_max_l-0.6123724356957945*energy_l[1]*u_max_l-0.6123724356957945*energy_c[1]*u_max_l-0.3535533905932737*pl[0]*u_max_l+0.3535533905932737*pc[0]*u_max_l-0.3535533905932737*energy_l[0]*u_max_l+0.3535533905932737*energy_c[0]*u_max_l+0.75*pl[1]*ul_0[1]+0.75*energy_l[1]*ul_0[1]+0.4330127018922193*pl[0]*ul_0[1]+0.4330127018922193*energy_l[0]*ul_0[1]+0.75*pc[1]*uc_0[1]+0.75*energy_c[1]*uc_0[1]-0.4330127018922193*pc[0]*uc_0[1]-0.4330127018922193*energy_c[0]*uc_0[1]+0.4330127018922193*ul_0[0]*pl[1]-0.4330127018922193*uc_0[0]*pc[1]+0.4330127018922193*ul_0[0]*energy_l[1]-0.4330127018922193*uc_0[0]*energy_c[1]+0.25*pl[0]*ul_0[0]+0.25*energy_l[0]*ul_0[0]+0.25*pc[0]*uc_0[0]+0.25*energy_c[0]*uc_0[0]; 
  double Ghat_energy_r = (-0.6123724356957945*pr[1]*u_max_r)-0.6123724356957945*pc[1]*u_max_r-0.6123724356957945*energy_r[1]*u_max_r-0.6123724356957945*energy_c[1]*u_max_r+0.3535533905932737*pr[0]*u_max_r-0.3535533905932737*pc[0]*u_max_r+0.3535533905932737*energy_r[0]*u_max_r-0.3535533905932737*energy_c[0]*u_max_r+0.75*pr[1]*ur_0[1]+0.75*energy_r[1]*ur_0[1]-0.4330127018922193*pr[0]*ur_0[1]-0.4330127018922193*energy_r[0]*ur_0[1]+0.75*pc[1]*uc_0[1]+0.75*energy_c[1]*uc_0[1]+0.4330127018922193*pc[0]*uc_0[1]+0.4330127018922193*energy_c[0]*uc_0[1]-0.4330127018922193*ur_0[0]*pr[1]+0.4330127018922193*uc_0[0]*pc[1]-0.4330127018922193*ur_0[0]*energy_r[1]+0.4330127018922193*uc_0[0]*energy_c[1]+0.25*pr[0]*ur_0[0]+0.25*energy_r[0]*ur_0[0]+0.25*pc[0]*uc_0[0]+0.25*energy_c[0]*uc_0[0]; 

  outrho[0] += 0.7071067811865475*Ghat_rho_l*dx1-0.7071067811865475*Ghat_rho_r*dx1; 
  outrho[1] += (-1.224744871391589*Ghat_rho_r*dx1)-1.224744871391589*Ghat_rho_l*dx1; 

  outrhou0[0] += 0.8660254037844386*ur_0[1]*Ghat_rho_r*dx1-0.8660254037844386*uc_0[1]*Ghat_rho_r*dx1-0.5*ur_0[0]*Ghat_rho_r*dx1-0.5*uc_0[0]*Ghat_rho_r*dx1+0.4330127018922193*ul_0[1]*Ghat_rho_l*dx1-0.4330127018922193*uc_0[1]*Ghat_rho_l*dx1+0.25*ul_0[0]*Ghat_rho_l*dx1+0.25*uc_0[0]*Ghat_rho_l*dx1+0.4330127018922193*pr[1]*dx1+0.4330127018922193*pl[1]*dx1-0.8660254037844386*pc[1]*dx1-0.25*pr[0]*dx1+0.25*pl[0]*dx1; 
  outrhou0[1] += 1.5*ur_0[1]*Ghat_rho_r*dx1-1.5*uc_0[1]*Ghat_rho_r*dx1-0.8660254037844386*ur_0[0]*Ghat_rho_r*dx1-0.8660254037844386*uc_0[0]*Ghat_rho_r*dx1-0.75*ul_0[1]*Ghat_rho_l*dx1+0.75*uc_0[1]*Ghat_rho_l*dx1-0.4330127018922193*ul_0[0]*Ghat_rho_l*dx1-0.4330127018922193*uc_0[0]*Ghat_rho_l*dx1+0.75*pr[1]*dx1-0.75*pl[1]*dx1-0.4330127018922193*pr[0]*dx1-0.4330127018922193*pl[0]*dx1-0.8660254037844386*pc[0]*dx1; 

  outrhou1[0] += 0.8660254037844386*ur_1[1]*Ghat_rho_r*dx1-0.8660254037844386*uc_1[1]*Ghat_rho_r*dx1-0.5*ur_1[0]*Ghat_rho_r*dx1-0.5*uc_1[0]*Ghat_rho_r*dx1+0.4330127018922193*ul_1[1]*Ghat_rho_l*dx1-0.4330127018922193*uc_1[1]*Ghat_rho_l*dx1+0.25*ul_1[0]*Ghat_rho_l*dx1+0.25*uc_1[0]*Ghat_rho_l*dx1; 
  outrhou1[1] += 1.5*ur_1[1]*Ghat_rho_r*dx1-1.5*uc_1[1]*Ghat_rho_r*dx1-0.8660254037844386*ur_1[0]*Ghat_rho_r*dx1-0.8660254037844386*uc_1[0]*Ghat_rho_r*dx1-0.75*ul_1[1]*Ghat_rho_l*dx1+0.75*uc_1[1]*Ghat_rho_l*dx1-0.4330127018922193*ul_1[0]*Ghat_rho_l*dx1-0.4330127018922193*uc_1[0]*Ghat_rho_l*dx1; 

  outrhou2[0] += 0.8660254037844386*ur_2[1]*Ghat_rho_r*dx1-0.8660254037844386*uc_2[1]*Ghat_rho_r*dx1-0.5*ur_2[0]*Ghat_rho_r*dx1-0.5*uc_2[0]*Ghat_rho_r*dx1+0.4330127018922193*ul_2[1]*Ghat_rho_l*dx1-0.4330127018922193*uc_2[1]*Ghat_rho_l*dx1+0.25*ul_2[0]*Ghat_rho_l*dx1+0.25*uc_2[0]*Ghat_rho_l*dx1; 
  outrhou2[1] += 1.5*ur_2[1]*Ghat_rho_r*dx1-1.5*uc_2[1]*Ghat_rho_r*dx1-0.8660254037844386*ur_2[0]*Ghat_rho_r*dx1-0.8660254037844386*uc_2[0]*Ghat_rho_r*dx1-0.75*ul_2[1]*Ghat_rho_l*dx1+0.75*uc_2[1]*Ghat_rho_l*dx1-0.4330127018922193*ul_2[0]*Ghat_rho_l*dx1-0.4330127018922193*uc_2[0]*Ghat_rho_l*dx1; 

  outenergy[0] += 0.7071067811865475*Ghat_energy_l*dx1-0.7071067811865475*Ghat_energy_r*dx1; 
  outenergy[1] += (-1.224744871391589*Ghat_energy_r*dx1)-1.224744871391589*Ghat_energy_l*dx1; 

} 
