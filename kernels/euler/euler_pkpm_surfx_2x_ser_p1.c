#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void euler_pkpm_surfx_2x_ser_p1(const double *w, const double *dxv, 
  const double *u_il, const double *u_ic, const double *u_ir,
  const double *p_ijl, const double *p_ijc, const double *p_ijr,
  const double *vlasov_pkpm_surf_momsl, const double *vlasov_pkpm_surf_momsr, 
  const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_il/u_ic/u_ir: [ux, uy, uz] Fluid flow in left/center/right cells.
  // p_ijl/p_ijc/p_ijr: Pressure tensor in left/center/right cells.
  // vlasov_pkpm_surf_momsl/vlasov_pkpm_surf_momsr: Mass flux and heat flux at left edge and right edge (computed externally) .
  // statevecl/statevecc/statevecr: [rho ux, rho uy, rho uz, energy], Fluid input state vector in left/center/right cells.
  // out: Incremented output.

  const double dx1 = 2.0/dxv[0]; 
  // Only need to fetch input energy variable, other fluxes are computed from input flow and pressure tensor.
  const double *energy_l = &statevecl[12]; 
  const double *energy_c = &statevecc[12]; 
  const double *energy_r = &statevecr[12]; 

  const double *ux_l = &u_il[0]; 
  const double *ux_c = &u_ic[0]; 
  const double *ux_r = &u_ir[0]; 

  const double *uy_l = &u_il[4]; 
  const double *uy_c = &u_ic[4]; 
  const double *uy_r = &u_ir[4]; 

  const double *uz_l = &u_il[8]; 
  const double *uz_c = &u_ic[8]; 
  const double *uz_r = &u_ir[8]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxx_r = &p_ijr[0]; 

  const double *Pxy_l = &p_ijl[4]; 
  const double *Pxy_c = &p_ijc[4]; 
  const double *Pxy_r = &p_ijr[4]; 

  const double *Pxz_l = &p_ijl[8]; 
  const double *Pxz_c = &p_ijc[8]; 
  const double *Pxz_r = &p_ijr[8]; 

  const double *Pyy_l = &p_ijl[12]; 
  const double *Pyy_c = &p_ijc[12]; 
  const double *Pyy_r = &p_ijr[12]; 

  const double *Pyz_l = &p_ijl[16]; 
  const double *Pyz_c = &p_ijc[16]; 
  const double *Pyz_r = &p_ijr[16]; 

  const double *Pzz_l = &p_ijl[20]; 
  const double *Pzz_c = &p_ijc[20]; 
  const double *Pzz_r = &p_ijr[20]; 

  const double *rho_flux_l = &vlasov_pkpm_surf_momsl[0]; 
  const double *heat_flux_l = &vlasov_pkpm_surf_momsl[2]; 
  const double *rho_flux_r = &vlasov_pkpm_surf_momsr[0]; 
  const double *heat_flux_r = &vlasov_pkpm_surf_momsr[2]; 
  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[4]; 
  double *outrhou2 = &out[8]; 
  double *outenergy = &out[12]; 

  double uxQuad_l[2] = {0.0};
  double uxQuad_r[2] = {0.0};
  double uxMax_l[2] = {0.0};;
  double uxMax_r[2] = {0.0};
  double uyQuad_l[2] = {0.0};
  double uyQuad_r[2] = {0.0};
  double uyMax_l[2] = {0.0};;
  double uyMax_r[2] = {0.0};
  double uzQuad_l[2] = {0.0};
  double uzQuad_r[2] = {0.0};
  double uzMax_l[2] = {0.0};;
  double uzMax_r[2] = {0.0};

  double Ghat_rhoux_l[2] = {0.0}; 
  double Ghat_rhoux_r[2] = {0.0}; 
  double Ghat_rhouy_l[2] = {0.0}; 
  double Ghat_rhouy_r[2] = {0.0}; 
  double Ghat_rhouz_l[2] = {0.0}; 
  double Ghat_rhouz_r[2] = {0.0}; 
  double Ghat_energy_l[2] = {0.0}; 
  double Ghat_energy_r[2] = {0.0}; 

  double ux_l_r = 0.0; 
  double ux_c_l = 0.0; 
  double ux_c_r = 0.0; 
  double ux_r_l = 0.0; 
  double uy_l_r = 0.0; 
  double uy_c_l = 0.0; 
  double uy_c_r = 0.0; 
  double uy_r_l = 0.0; 
  double uz_l_r = 0.0; 
  double uz_c_l = 0.0; 
  double uz_c_r = 0.0; 
  double uz_r_l = 0.0; 

  ux_l_r = ser_2x_p1_surfx1_eval_quad_node_0_r(ux_l); 
  ux_c_l = ser_2x_p1_surfx1_eval_quad_node_0_l(ux_c); 
  ux_c_r = ser_2x_p1_surfx1_eval_quad_node_0_r(ux_c); 
  ux_r_l = ser_2x_p1_surfx1_eval_quad_node_0_l(ux_r); 
  uxQuad_l[0] = fmax(fabs(ux_l_r), fabs(ux_c_l)); 
  uxQuad_r[0] = fmax(fabs(ux_c_r), fabs(ux_r_l)); 

  uy_l_r = ser_2x_p1_surfx1_eval_quad_node_0_r(uy_l); 
  uy_c_l = ser_2x_p1_surfx1_eval_quad_node_0_l(uy_c); 
  uy_c_r = ser_2x_p1_surfx1_eval_quad_node_0_r(uy_c); 
  uy_r_l = ser_2x_p1_surfx1_eval_quad_node_0_l(uy_r); 
  uyQuad_l[0] = fmax(fabs(uy_l_r), fabs(uy_c_l)); 
  uyQuad_r[0] = fmax(fabs(uy_c_r), fabs(uy_r_l)); 

  uz_l_r = ser_2x_p1_surfx1_eval_quad_node_0_r(uz_l); 
  uz_c_l = ser_2x_p1_surfx1_eval_quad_node_0_l(uz_c); 
  uz_c_r = ser_2x_p1_surfx1_eval_quad_node_0_r(uz_c); 
  uz_r_l = ser_2x_p1_surfx1_eval_quad_node_0_l(uz_r); 
  uzQuad_l[0] = fmax(fabs(uz_l_r), fabs(uz_c_l)); 
  uzQuad_r[0] = fmax(fabs(uz_c_r), fabs(uz_r_l)); 

  ux_l_r = ser_2x_p1_surfx1_eval_quad_node_1_r(ux_l); 
  ux_c_l = ser_2x_p1_surfx1_eval_quad_node_1_l(ux_c); 
  ux_c_r = ser_2x_p1_surfx1_eval_quad_node_1_r(ux_c); 
  ux_r_l = ser_2x_p1_surfx1_eval_quad_node_1_l(ux_r); 
  uxQuad_l[1] = fmax(fabs(ux_l_r), fabs(ux_c_l)); 
  uxQuad_r[1] = fmax(fabs(ux_c_r), fabs(ux_r_l)); 

  uy_l_r = ser_2x_p1_surfx1_eval_quad_node_1_r(uy_l); 
  uy_c_l = ser_2x_p1_surfx1_eval_quad_node_1_l(uy_c); 
  uy_c_r = ser_2x_p1_surfx1_eval_quad_node_1_r(uy_c); 
  uy_r_l = ser_2x_p1_surfx1_eval_quad_node_1_l(uy_r); 
  uyQuad_l[1] = fmax(fabs(uy_l_r), fabs(uy_c_l)); 
  uyQuad_r[1] = fmax(fabs(uy_c_r), fabs(uy_r_l)); 

  uz_l_r = ser_2x_p1_surfx1_eval_quad_node_1_r(uz_l); 
  uz_c_l = ser_2x_p1_surfx1_eval_quad_node_1_l(uz_c); 
  uz_c_r = ser_2x_p1_surfx1_eval_quad_node_1_r(uz_c); 
  uz_r_l = ser_2x_p1_surfx1_eval_quad_node_1_l(uz_r); 
  uzQuad_l[1] = fmax(fabs(uz_l_r), fabs(uz_c_l)); 
  uzQuad_r[1] = fmax(fabs(uz_c_r), fabs(uz_r_l)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p1_upwind_quad_to_modal(uxQuad_l, uxMax_l); 
  ser_2x_p1_upwind_quad_to_modal(uxQuad_r, uxMax_r); 

  ser_2x_p1_upwind_quad_to_modal(uyQuad_l, uyMax_l); 
  ser_2x_p1_upwind_quad_to_modal(uyQuad_r, uyMax_r); 

  ser_2x_p1_upwind_quad_to_modal(uzQuad_l, uzMax_l); 
  ser_2x_p1_upwind_quad_to_modal(uzQuad_r, uzMax_r); 

  Ghat_rhoux_l[0] = 0.4330127018922193*rho_flux_l[1]*ux_l[3]-0.4330127018922193*rho_flux_l[1]*ux_c[3]+0.25*rho_flux_l[1]*ux_l[2]+0.25*rho_flux_l[1]*ux_c[2]+0.4330127018922193*rho_flux_l[0]*ux_l[1]-0.4330127018922193*rho_flux_l[0]*ux_c[1]+0.6123724356957944*Pxx_l[1]-0.6123724356957944*Pxx_c[1]+0.25*rho_flux_l[0]*ux_l[0]+0.25*rho_flux_l[0]*ux_c[0]+0.3535533905932737*Pxx_l[0]+0.3535533905932737*Pxx_c[0]; 
  Ghat_rhoux_l[1] = 0.4330127018922193*rho_flux_l[0]*ux_l[3]-0.4330127018922193*rho_flux_l[0]*ux_c[3]+0.6123724356957944*Pxx_l[3]-0.6123724356957944*Pxx_c[3]+0.25*rho_flux_l[0]*ux_l[2]+0.25*rho_flux_l[0]*ux_c[2]+0.3535533905932737*Pxx_l[2]+0.3535533905932737*Pxx_c[2]+0.4330127018922193*rho_flux_l[1]*ux_l[1]-0.4330127018922193*rho_flux_l[1]*ux_c[1]+0.25*ux_l[0]*rho_flux_l[1]+0.25*ux_c[0]*rho_flux_l[1]; 

  Ghat_rhouy_l[0] = 0.4330127018922193*rho_flux_l[1]*uy_l[3]-0.4330127018922193*rho_flux_l[1]*uy_c[3]+0.25*rho_flux_l[1]*uy_l[2]+0.25*rho_flux_l[1]*uy_c[2]+0.4330127018922193*rho_flux_l[0]*uy_l[1]-0.4330127018922193*rho_flux_l[0]*uy_c[1]+0.6123724356957944*Pxy_l[1]-0.6123724356957944*Pxy_c[1]+0.25*rho_flux_l[0]*uy_l[0]+0.25*rho_flux_l[0]*uy_c[0]+0.3535533905932737*Pxy_l[0]+0.3535533905932737*Pxy_c[0]; 
  Ghat_rhouy_l[1] = 0.4330127018922193*rho_flux_l[0]*uy_l[3]-0.4330127018922193*rho_flux_l[0]*uy_c[3]+0.6123724356957944*Pxy_l[3]-0.6123724356957944*Pxy_c[3]+0.25*rho_flux_l[0]*uy_l[2]+0.25*rho_flux_l[0]*uy_c[2]+0.3535533905932737*Pxy_l[2]+0.3535533905932737*Pxy_c[2]+0.4330127018922193*rho_flux_l[1]*uy_l[1]-0.4330127018922193*rho_flux_l[1]*uy_c[1]+0.25*uy_l[0]*rho_flux_l[1]+0.25*uy_c[0]*rho_flux_l[1]; 

  Ghat_rhouz_l[0] = 0.4330127018922193*rho_flux_l[1]*uz_l[3]-0.4330127018922193*rho_flux_l[1]*uz_c[3]+0.25*rho_flux_l[1]*uz_l[2]+0.25*rho_flux_l[1]*uz_c[2]+0.4330127018922193*rho_flux_l[0]*uz_l[1]-0.4330127018922193*rho_flux_l[0]*uz_c[1]+0.6123724356957944*Pxz_l[1]-0.6123724356957944*Pxz_c[1]+0.25*rho_flux_l[0]*uz_l[0]+0.25*rho_flux_l[0]*uz_c[0]+0.3535533905932737*Pxz_l[0]+0.3535533905932737*Pxz_c[0]; 
  Ghat_rhouz_l[1] = 0.4330127018922193*rho_flux_l[0]*uz_l[3]-0.4330127018922193*rho_flux_l[0]*uz_c[3]+0.6123724356957944*Pxz_l[3]-0.6123724356957944*Pxz_c[3]+0.25*rho_flux_l[0]*uz_l[2]+0.25*rho_flux_l[0]*uz_c[2]+0.3535533905932737*Pxz_l[2]+0.3535533905932737*Pxz_c[2]+0.4330127018922193*rho_flux_l[1]*uz_l[1]-0.4330127018922193*rho_flux_l[1]*uz_c[1]+0.25*uz_l[0]*rho_flux_l[1]+0.25*uz_c[0]*rho_flux_l[1]; 

  Ghat_energy_l[0] = 0.5303300858899105*Pxz_l[3]*uz_l[3]+0.3061862178478971*Pxz_l[2]*uz_l[3]+0.5303300858899105*Pxz_c[3]*uz_c[3]-0.3061862178478971*Pxz_c[2]*uz_c[3]+0.5303300858899105*Pxy_l[3]*uy_l[3]+0.3061862178478971*Pxy_l[2]*uy_l[3]+0.5303300858899105*Pxy_c[3]*uy_c[3]-0.3061862178478971*Pxy_c[2]*uy_c[3]+0.5303300858899105*energy_l[3]*ux_l[3]+0.5303300858899105*Pxx_l[3]*ux_l[3]+0.3061862178478971*energy_l[2]*ux_l[3]+0.3061862178478971*Pxx_l[2]*ux_l[3]+0.5303300858899105*energy_c[3]*ux_c[3]+0.5303300858899105*Pxx_c[3]*ux_c[3]-0.3061862178478971*energy_c[2]*ux_c[3]-0.3061862178478971*Pxx_c[2]*ux_c[3]+0.3061862178478971*ux_l[2]*energy_l[3]-0.4330127018922193*uxMax_l[1]*energy_l[3]-0.3061862178478971*ux_c[2]*energy_c[3]-0.4330127018922193*uxMax_l[1]*energy_c[3]+0.3061862178478971*uz_l[2]*Pxz_l[3]-0.4330127018922193*uzMax_l[1]*Pxz_l[3]-0.3061862178478971*uz_c[2]*Pxz_c[3]-0.4330127018922193*uzMax_l[1]*Pxz_c[3]+0.3061862178478971*uy_l[2]*Pxy_l[3]-0.4330127018922193*uyMax_l[1]*Pxy_l[3]-0.3061862178478971*uy_c[2]*Pxy_c[3]-0.4330127018922193*uyMax_l[1]*Pxy_c[3]+0.3061862178478971*ux_l[2]*Pxx_l[3]-0.4330127018922193*uxMax_l[1]*Pxx_l[3]-0.3061862178478971*ux_c[2]*Pxx_c[3]-0.4330127018922193*uxMax_l[1]*Pxx_c[3]+0.1767766952966368*Pxz_l[2]*uz_l[2]+0.1767766952966368*Pxz_c[2]*uz_c[2]+0.1767766952966368*Pxy_l[2]*uy_l[2]+0.1767766952966368*Pxy_c[2]*uy_c[2]+0.1767766952966368*energy_l[2]*ux_l[2]+0.1767766952966368*Pxx_l[2]*ux_l[2]+0.1767766952966368*energy_c[2]*ux_c[2]+0.1767766952966368*Pxx_c[2]*ux_c[2]-0.25*uxMax_l[1]*energy_l[2]+0.25*uxMax_l[1]*energy_c[2]-0.25*uzMax_l[1]*Pxz_l[2]+0.25*uzMax_l[1]*Pxz_c[2]-0.25*uyMax_l[1]*Pxy_l[2]+0.25*uyMax_l[1]*Pxy_c[2]-0.25*uxMax_l[1]*Pxx_l[2]+0.25*uxMax_l[1]*Pxx_c[2]+0.5303300858899105*Pxz_l[1]*uz_l[1]+0.3061862178478971*Pxz_l[0]*uz_l[1]+0.5303300858899105*Pxz_c[1]*uz_c[1]-0.3061862178478971*Pxz_c[0]*uz_c[1]+0.5303300858899105*Pxy_l[1]*uy_l[1]+0.3061862178478971*Pxy_l[0]*uy_l[1]+0.5303300858899105*Pxy_c[1]*uy_c[1]-0.3061862178478971*Pxy_c[0]*uy_c[1]+0.5303300858899105*energy_l[1]*ux_l[1]+0.5303300858899105*Pxx_l[1]*ux_l[1]+0.3061862178478971*energy_l[0]*ux_l[1]+0.3061862178478971*Pxx_l[0]*ux_l[1]+0.5303300858899105*energy_c[1]*ux_c[1]+0.5303300858899105*Pxx_c[1]*ux_c[1]-0.3061862178478971*energy_c[0]*ux_c[1]-0.3061862178478971*Pxx_c[0]*ux_c[1]+0.3061862178478971*ux_l[0]*energy_l[1]-0.4330127018922193*uxMax_l[0]*energy_l[1]-0.3061862178478971*ux_c[0]*energy_c[1]-0.4330127018922193*uxMax_l[0]*energy_c[1]+0.3061862178478971*uz_l[0]*Pxz_l[1]-0.4330127018922193*uzMax_l[0]*Pxz_l[1]-0.3061862178478971*uz_c[0]*Pxz_c[1]-0.4330127018922193*uzMax_l[0]*Pxz_c[1]+0.3061862178478971*uy_l[0]*Pxy_l[1]-0.4330127018922193*uyMax_l[0]*Pxy_l[1]-0.3061862178478971*uy_c[0]*Pxy_c[1]-0.4330127018922193*uyMax_l[0]*Pxy_c[1]+0.3061862178478971*ux_l[0]*Pxx_l[1]-0.4330127018922193*uxMax_l[0]*Pxx_l[1]-0.3061862178478971*ux_c[0]*Pxx_c[1]-0.4330127018922193*uxMax_l[0]*Pxx_c[1]+0.1767766952966368*Pxz_l[0]*uz_l[0]+0.1767766952966368*Pxz_c[0]*uz_c[0]-0.25*Pxz_l[0]*uzMax_l[0]+0.25*Pxz_c[0]*uzMax_l[0]+0.1767766952966368*Pxy_l[0]*uy_l[0]+0.1767766952966368*Pxy_c[0]*uy_c[0]-0.25*Pxy_l[0]*uyMax_l[0]+0.25*Pxy_c[0]*uyMax_l[0]+0.1767766952966368*energy_l[0]*ux_l[0]+0.1767766952966368*Pxx_l[0]*ux_l[0]+0.1767766952966368*energy_c[0]*ux_c[0]+0.1767766952966368*Pxx_c[0]*ux_c[0]-0.25*energy_l[0]*uxMax_l[0]+0.25*energy_c[0]*uxMax_l[0]-0.25*Pxx_l[0]*uxMax_l[0]+0.25*Pxx_c[0]*uxMax_l[0]; 
  Ghat_energy_l[1] = 0.5303300858899105*Pxz_l[1]*uz_l[3]+0.3061862178478971*Pxz_l[0]*uz_l[3]+0.5303300858899105*Pxz_c[1]*uz_c[3]-0.3061862178478971*Pxz_c[0]*uz_c[3]+0.5303300858899105*Pxy_l[1]*uy_l[3]+0.3061862178478971*Pxy_l[0]*uy_l[3]+0.5303300858899105*Pxy_c[1]*uy_c[3]-0.3061862178478971*Pxy_c[0]*uy_c[3]+0.5303300858899105*energy_l[1]*ux_l[3]+0.5303300858899105*Pxx_l[1]*ux_l[3]+0.3061862178478971*energy_l[0]*ux_l[3]+0.3061862178478971*Pxx_l[0]*ux_l[3]+0.5303300858899105*energy_c[1]*ux_c[3]+0.5303300858899105*Pxx_c[1]*ux_c[3]-0.3061862178478971*energy_c[0]*ux_c[3]-0.3061862178478971*Pxx_c[0]*ux_c[3]+0.5303300858899105*ux_l[1]*energy_l[3]+0.3061862178478971*ux_l[0]*energy_l[3]-0.4330127018922193*uxMax_l[0]*energy_l[3]+0.5303300858899105*ux_c[1]*energy_c[3]-0.3061862178478971*ux_c[0]*energy_c[3]-0.4330127018922193*uxMax_l[0]*energy_c[3]+0.5303300858899105*uz_l[1]*Pxz_l[3]+0.3061862178478971*uz_l[0]*Pxz_l[3]-0.4330127018922193*uzMax_l[0]*Pxz_l[3]+0.5303300858899105*uz_c[1]*Pxz_c[3]-0.3061862178478971*uz_c[0]*Pxz_c[3]-0.4330127018922193*uzMax_l[0]*Pxz_c[3]+0.5303300858899105*uy_l[1]*Pxy_l[3]+0.3061862178478971*uy_l[0]*Pxy_l[3]-0.4330127018922193*uyMax_l[0]*Pxy_l[3]+0.5303300858899105*uy_c[1]*Pxy_c[3]-0.3061862178478971*uy_c[0]*Pxy_c[3]-0.4330127018922193*uyMax_l[0]*Pxy_c[3]+0.5303300858899105*ux_l[1]*Pxx_l[3]+0.3061862178478971*ux_l[0]*Pxx_l[3]-0.4330127018922193*uxMax_l[0]*Pxx_l[3]+0.5303300858899105*ux_c[1]*Pxx_c[3]-0.3061862178478971*ux_c[0]*Pxx_c[3]-0.4330127018922193*uxMax_l[0]*Pxx_c[3]+0.3061862178478971*Pxz_l[1]*uz_l[2]+0.1767766952966368*Pxz_l[0]*uz_l[2]-0.3061862178478971*Pxz_c[1]*uz_c[2]+0.1767766952966368*Pxz_c[0]*uz_c[2]+0.3061862178478971*Pxy_l[1]*uy_l[2]+0.1767766952966368*Pxy_l[0]*uy_l[2]-0.3061862178478971*Pxy_c[1]*uy_c[2]+0.1767766952966368*Pxy_c[0]*uy_c[2]+0.3061862178478971*energy_l[1]*ux_l[2]+0.3061862178478971*Pxx_l[1]*ux_l[2]+0.1767766952966368*energy_l[0]*ux_l[2]+0.1767766952966368*Pxx_l[0]*ux_l[2]-0.3061862178478971*energy_c[1]*ux_c[2]-0.3061862178478971*Pxx_c[1]*ux_c[2]+0.1767766952966368*energy_c[0]*ux_c[2]+0.1767766952966368*Pxx_c[0]*ux_c[2]+0.3061862178478971*ux_l[1]*energy_l[2]+0.1767766952966368*ux_l[0]*energy_l[2]-0.25*uxMax_l[0]*energy_l[2]-0.3061862178478971*ux_c[1]*energy_c[2]+0.1767766952966368*ux_c[0]*energy_c[2]+0.25*uxMax_l[0]*energy_c[2]+0.3061862178478971*uz_l[1]*Pxz_l[2]+0.1767766952966368*uz_l[0]*Pxz_l[2]-0.25*uzMax_l[0]*Pxz_l[2]-0.3061862178478971*uz_c[1]*Pxz_c[2]+0.1767766952966368*uz_c[0]*Pxz_c[2]+0.25*uzMax_l[0]*Pxz_c[2]+0.3061862178478971*uy_l[1]*Pxy_l[2]+0.1767766952966368*uy_l[0]*Pxy_l[2]-0.25*uyMax_l[0]*Pxy_l[2]-0.3061862178478971*uy_c[1]*Pxy_c[2]+0.1767766952966368*uy_c[0]*Pxy_c[2]+0.25*uyMax_l[0]*Pxy_c[2]+0.3061862178478971*ux_l[1]*Pxx_l[2]+0.1767766952966368*ux_l[0]*Pxx_l[2]-0.25*uxMax_l[0]*Pxx_l[2]-0.3061862178478971*ux_c[1]*Pxx_c[2]+0.1767766952966368*ux_c[0]*Pxx_c[2]+0.25*uxMax_l[0]*Pxx_c[2]-0.4330127018922193*Pxz_l[1]*uzMax_l[1]-0.4330127018922193*Pxz_c[1]*uzMax_l[1]-0.25*Pxz_l[0]*uzMax_l[1]+0.25*Pxz_c[0]*uzMax_l[1]-0.4330127018922193*Pxy_l[1]*uyMax_l[1]-0.4330127018922193*Pxy_c[1]*uyMax_l[1]-0.25*Pxy_l[0]*uyMax_l[1]+0.25*Pxy_c[0]*uyMax_l[1]-0.4330127018922193*energy_l[1]*uxMax_l[1]-0.4330127018922193*energy_c[1]*uxMax_l[1]-0.4330127018922193*Pxx_l[1]*uxMax_l[1]-0.4330127018922193*Pxx_c[1]*uxMax_l[1]-0.25*energy_l[0]*uxMax_l[1]+0.25*energy_c[0]*uxMax_l[1]-0.25*Pxx_l[0]*uxMax_l[1]+0.25*Pxx_c[0]*uxMax_l[1]; 

  Ghat_rhoux_r[0] = (-0.4330127018922193*rho_flux_r[1]*ux_r[3])+0.4330127018922193*rho_flux_r[1]*ux_c[3]+0.25*rho_flux_r[1]*ux_r[2]+0.25*rho_flux_r[1]*ux_c[2]-0.4330127018922193*rho_flux_r[0]*ux_r[1]+0.4330127018922193*rho_flux_r[0]*ux_c[1]-0.6123724356957944*Pxx_r[1]+0.6123724356957944*Pxx_c[1]+0.25*rho_flux_r[0]*ux_r[0]+0.25*rho_flux_r[0]*ux_c[0]+0.3535533905932737*Pxx_r[0]+0.3535533905932737*Pxx_c[0]; 
  Ghat_rhoux_r[1] = (-0.4330127018922193*rho_flux_r[0]*ux_r[3])+0.4330127018922193*rho_flux_r[0]*ux_c[3]-0.6123724356957944*Pxx_r[3]+0.6123724356957944*Pxx_c[3]+0.25*rho_flux_r[0]*ux_r[2]+0.25*rho_flux_r[0]*ux_c[2]+0.3535533905932737*Pxx_r[2]+0.3535533905932737*Pxx_c[2]-0.4330127018922193*rho_flux_r[1]*ux_r[1]+0.4330127018922193*rho_flux_r[1]*ux_c[1]+0.25*ux_r[0]*rho_flux_r[1]+0.25*ux_c[0]*rho_flux_r[1]; 

  Ghat_rhouy_r[0] = (-0.4330127018922193*rho_flux_r[1]*uy_r[3])+0.4330127018922193*rho_flux_r[1]*uy_c[3]+0.25*rho_flux_r[1]*uy_r[2]+0.25*rho_flux_r[1]*uy_c[2]-0.4330127018922193*rho_flux_r[0]*uy_r[1]+0.4330127018922193*rho_flux_r[0]*uy_c[1]-0.6123724356957944*Pxy_r[1]+0.6123724356957944*Pxy_c[1]+0.25*rho_flux_r[0]*uy_r[0]+0.25*rho_flux_r[0]*uy_c[0]+0.3535533905932737*Pxy_r[0]+0.3535533905932737*Pxy_c[0]; 
  Ghat_rhouy_r[1] = (-0.4330127018922193*rho_flux_r[0]*uy_r[3])+0.4330127018922193*rho_flux_r[0]*uy_c[3]-0.6123724356957944*Pxy_r[3]+0.6123724356957944*Pxy_c[3]+0.25*rho_flux_r[0]*uy_r[2]+0.25*rho_flux_r[0]*uy_c[2]+0.3535533905932737*Pxy_r[2]+0.3535533905932737*Pxy_c[2]-0.4330127018922193*rho_flux_r[1]*uy_r[1]+0.4330127018922193*rho_flux_r[1]*uy_c[1]+0.25*uy_r[0]*rho_flux_r[1]+0.25*uy_c[0]*rho_flux_r[1]; 

  Ghat_rhouz_r[0] = (-0.4330127018922193*rho_flux_r[1]*uz_r[3])+0.4330127018922193*rho_flux_r[1]*uz_c[3]+0.25*rho_flux_r[1]*uz_r[2]+0.25*rho_flux_r[1]*uz_c[2]-0.4330127018922193*rho_flux_r[0]*uz_r[1]+0.4330127018922193*rho_flux_r[0]*uz_c[1]-0.6123724356957944*Pxz_r[1]+0.6123724356957944*Pxz_c[1]+0.25*rho_flux_r[0]*uz_r[0]+0.25*rho_flux_r[0]*uz_c[0]+0.3535533905932737*Pxz_r[0]+0.3535533905932737*Pxz_c[0]; 
  Ghat_rhouz_r[1] = (-0.4330127018922193*rho_flux_r[0]*uz_r[3])+0.4330127018922193*rho_flux_r[0]*uz_c[3]-0.6123724356957944*Pxz_r[3]+0.6123724356957944*Pxz_c[3]+0.25*rho_flux_r[0]*uz_r[2]+0.25*rho_flux_r[0]*uz_c[2]+0.3535533905932737*Pxz_r[2]+0.3535533905932737*Pxz_c[2]-0.4330127018922193*rho_flux_r[1]*uz_r[1]+0.4330127018922193*rho_flux_r[1]*uz_c[1]+0.25*uz_r[0]*rho_flux_r[1]+0.25*uz_c[0]*rho_flux_r[1]; 

  Ghat_energy_r[0] = 0.5303300858899105*Pxz_r[3]*uz_r[3]-0.3061862178478971*Pxz_r[2]*uz_r[3]+0.5303300858899105*Pxz_c[3]*uz_c[3]+0.3061862178478971*Pxz_c[2]*uz_c[3]+0.5303300858899105*Pxy_r[3]*uy_r[3]-0.3061862178478971*Pxy_r[2]*uy_r[3]+0.5303300858899105*Pxy_c[3]*uy_c[3]+0.3061862178478971*Pxy_c[2]*uy_c[3]+0.5303300858899105*energy_r[3]*ux_r[3]+0.5303300858899105*Pxx_r[3]*ux_r[3]-0.3061862178478971*energy_r[2]*ux_r[3]-0.3061862178478971*Pxx_r[2]*ux_r[3]+0.5303300858899105*energy_c[3]*ux_c[3]+0.5303300858899105*Pxx_c[3]*ux_c[3]+0.3061862178478971*energy_c[2]*ux_c[3]+0.3061862178478971*Pxx_c[2]*ux_c[3]-0.3061862178478971*ux_r[2]*energy_r[3]-0.4330127018922193*uxMax_r[1]*energy_r[3]+0.3061862178478971*ux_c[2]*energy_c[3]-0.4330127018922193*uxMax_r[1]*energy_c[3]-0.3061862178478971*uz_r[2]*Pxz_r[3]-0.4330127018922193*uzMax_r[1]*Pxz_r[3]+0.3061862178478971*uz_c[2]*Pxz_c[3]-0.4330127018922193*uzMax_r[1]*Pxz_c[3]-0.3061862178478971*uy_r[2]*Pxy_r[3]-0.4330127018922193*uyMax_r[1]*Pxy_r[3]+0.3061862178478971*uy_c[2]*Pxy_c[3]-0.4330127018922193*uyMax_r[1]*Pxy_c[3]-0.3061862178478971*ux_r[2]*Pxx_r[3]-0.4330127018922193*uxMax_r[1]*Pxx_r[3]+0.3061862178478971*ux_c[2]*Pxx_c[3]-0.4330127018922193*uxMax_r[1]*Pxx_c[3]+0.1767766952966368*Pxz_r[2]*uz_r[2]+0.1767766952966368*Pxz_c[2]*uz_c[2]+0.1767766952966368*Pxy_r[2]*uy_r[2]+0.1767766952966368*Pxy_c[2]*uy_c[2]+0.1767766952966368*energy_r[2]*ux_r[2]+0.1767766952966368*Pxx_r[2]*ux_r[2]+0.1767766952966368*energy_c[2]*ux_c[2]+0.1767766952966368*Pxx_c[2]*ux_c[2]+0.25*uxMax_r[1]*energy_r[2]-0.25*uxMax_r[1]*energy_c[2]+0.25*uzMax_r[1]*Pxz_r[2]-0.25*uzMax_r[1]*Pxz_c[2]+0.25*uyMax_r[1]*Pxy_r[2]-0.25*uyMax_r[1]*Pxy_c[2]+0.25*uxMax_r[1]*Pxx_r[2]-0.25*uxMax_r[1]*Pxx_c[2]+0.5303300858899105*Pxz_r[1]*uz_r[1]-0.3061862178478971*Pxz_r[0]*uz_r[1]+0.5303300858899105*Pxz_c[1]*uz_c[1]+0.3061862178478971*Pxz_c[0]*uz_c[1]+0.5303300858899105*Pxy_r[1]*uy_r[1]-0.3061862178478971*Pxy_r[0]*uy_r[1]+0.5303300858899105*Pxy_c[1]*uy_c[1]+0.3061862178478971*Pxy_c[0]*uy_c[1]+0.5303300858899105*energy_r[1]*ux_r[1]+0.5303300858899105*Pxx_r[1]*ux_r[1]-0.3061862178478971*energy_r[0]*ux_r[1]-0.3061862178478971*Pxx_r[0]*ux_r[1]+0.5303300858899105*energy_c[1]*ux_c[1]+0.5303300858899105*Pxx_c[1]*ux_c[1]+0.3061862178478971*energy_c[0]*ux_c[1]+0.3061862178478971*Pxx_c[0]*ux_c[1]-0.3061862178478971*ux_r[0]*energy_r[1]-0.4330127018922193*uxMax_r[0]*energy_r[1]+0.3061862178478971*ux_c[0]*energy_c[1]-0.4330127018922193*uxMax_r[0]*energy_c[1]-0.3061862178478971*uz_r[0]*Pxz_r[1]-0.4330127018922193*uzMax_r[0]*Pxz_r[1]+0.3061862178478971*uz_c[0]*Pxz_c[1]-0.4330127018922193*uzMax_r[0]*Pxz_c[1]-0.3061862178478971*uy_r[0]*Pxy_r[1]-0.4330127018922193*uyMax_r[0]*Pxy_r[1]+0.3061862178478971*uy_c[0]*Pxy_c[1]-0.4330127018922193*uyMax_r[0]*Pxy_c[1]-0.3061862178478971*ux_r[0]*Pxx_r[1]-0.4330127018922193*uxMax_r[0]*Pxx_r[1]+0.3061862178478971*ux_c[0]*Pxx_c[1]-0.4330127018922193*uxMax_r[0]*Pxx_c[1]+0.1767766952966368*Pxz_r[0]*uz_r[0]+0.1767766952966368*Pxz_c[0]*uz_c[0]+0.25*Pxz_r[0]*uzMax_r[0]-0.25*Pxz_c[0]*uzMax_r[0]+0.1767766952966368*Pxy_r[0]*uy_r[0]+0.1767766952966368*Pxy_c[0]*uy_c[0]+0.25*Pxy_r[0]*uyMax_r[0]-0.25*Pxy_c[0]*uyMax_r[0]+0.1767766952966368*energy_r[0]*ux_r[0]+0.1767766952966368*Pxx_r[0]*ux_r[0]+0.1767766952966368*energy_c[0]*ux_c[0]+0.1767766952966368*Pxx_c[0]*ux_c[0]+0.25*energy_r[0]*uxMax_r[0]-0.25*energy_c[0]*uxMax_r[0]+0.25*Pxx_r[0]*uxMax_r[0]-0.25*Pxx_c[0]*uxMax_r[0]; 
  Ghat_energy_r[1] = 0.5303300858899105*Pxz_r[1]*uz_r[3]-0.3061862178478971*Pxz_r[0]*uz_r[3]+0.5303300858899105*Pxz_c[1]*uz_c[3]+0.3061862178478971*Pxz_c[0]*uz_c[3]+0.5303300858899105*Pxy_r[1]*uy_r[3]-0.3061862178478971*Pxy_r[0]*uy_r[3]+0.5303300858899105*Pxy_c[1]*uy_c[3]+0.3061862178478971*Pxy_c[0]*uy_c[3]+0.5303300858899105*energy_r[1]*ux_r[3]+0.5303300858899105*Pxx_r[1]*ux_r[3]-0.3061862178478971*energy_r[0]*ux_r[3]-0.3061862178478971*Pxx_r[0]*ux_r[3]+0.5303300858899105*energy_c[1]*ux_c[3]+0.5303300858899105*Pxx_c[1]*ux_c[3]+0.3061862178478971*energy_c[0]*ux_c[3]+0.3061862178478971*Pxx_c[0]*ux_c[3]+0.5303300858899105*ux_r[1]*energy_r[3]-0.3061862178478971*ux_r[0]*energy_r[3]-0.4330127018922193*uxMax_r[0]*energy_r[3]+0.5303300858899105*ux_c[1]*energy_c[3]+0.3061862178478971*ux_c[0]*energy_c[3]-0.4330127018922193*uxMax_r[0]*energy_c[3]+0.5303300858899105*uz_r[1]*Pxz_r[3]-0.3061862178478971*uz_r[0]*Pxz_r[3]-0.4330127018922193*uzMax_r[0]*Pxz_r[3]+0.5303300858899105*uz_c[1]*Pxz_c[3]+0.3061862178478971*uz_c[0]*Pxz_c[3]-0.4330127018922193*uzMax_r[0]*Pxz_c[3]+0.5303300858899105*uy_r[1]*Pxy_r[3]-0.3061862178478971*uy_r[0]*Pxy_r[3]-0.4330127018922193*uyMax_r[0]*Pxy_r[3]+0.5303300858899105*uy_c[1]*Pxy_c[3]+0.3061862178478971*uy_c[0]*Pxy_c[3]-0.4330127018922193*uyMax_r[0]*Pxy_c[3]+0.5303300858899105*ux_r[1]*Pxx_r[3]-0.3061862178478971*ux_r[0]*Pxx_r[3]-0.4330127018922193*uxMax_r[0]*Pxx_r[3]+0.5303300858899105*ux_c[1]*Pxx_c[3]+0.3061862178478971*ux_c[0]*Pxx_c[3]-0.4330127018922193*uxMax_r[0]*Pxx_c[3]-0.3061862178478971*Pxz_r[1]*uz_r[2]+0.1767766952966368*Pxz_r[0]*uz_r[2]+0.3061862178478971*Pxz_c[1]*uz_c[2]+0.1767766952966368*Pxz_c[0]*uz_c[2]-0.3061862178478971*Pxy_r[1]*uy_r[2]+0.1767766952966368*Pxy_r[0]*uy_r[2]+0.3061862178478971*Pxy_c[1]*uy_c[2]+0.1767766952966368*Pxy_c[0]*uy_c[2]-0.3061862178478971*energy_r[1]*ux_r[2]-0.3061862178478971*Pxx_r[1]*ux_r[2]+0.1767766952966368*energy_r[0]*ux_r[2]+0.1767766952966368*Pxx_r[0]*ux_r[2]+0.3061862178478971*energy_c[1]*ux_c[2]+0.3061862178478971*Pxx_c[1]*ux_c[2]+0.1767766952966368*energy_c[0]*ux_c[2]+0.1767766952966368*Pxx_c[0]*ux_c[2]-0.3061862178478971*ux_r[1]*energy_r[2]+0.1767766952966368*ux_r[0]*energy_r[2]+0.25*uxMax_r[0]*energy_r[2]+0.3061862178478971*ux_c[1]*energy_c[2]+0.1767766952966368*ux_c[0]*energy_c[2]-0.25*uxMax_r[0]*energy_c[2]-0.3061862178478971*uz_r[1]*Pxz_r[2]+0.1767766952966368*uz_r[0]*Pxz_r[2]+0.25*uzMax_r[0]*Pxz_r[2]+0.3061862178478971*uz_c[1]*Pxz_c[2]+0.1767766952966368*uz_c[0]*Pxz_c[2]-0.25*uzMax_r[0]*Pxz_c[2]-0.3061862178478971*uy_r[1]*Pxy_r[2]+0.1767766952966368*uy_r[0]*Pxy_r[2]+0.25*uyMax_r[0]*Pxy_r[2]+0.3061862178478971*uy_c[1]*Pxy_c[2]+0.1767766952966368*uy_c[0]*Pxy_c[2]-0.25*uyMax_r[0]*Pxy_c[2]-0.3061862178478971*ux_r[1]*Pxx_r[2]+0.1767766952966368*ux_r[0]*Pxx_r[2]+0.25*uxMax_r[0]*Pxx_r[2]+0.3061862178478971*ux_c[1]*Pxx_c[2]+0.1767766952966368*ux_c[0]*Pxx_c[2]-0.25*uxMax_r[0]*Pxx_c[2]-0.4330127018922193*Pxz_r[1]*uzMax_r[1]-0.4330127018922193*Pxz_c[1]*uzMax_r[1]+0.25*Pxz_r[0]*uzMax_r[1]-0.25*Pxz_c[0]*uzMax_r[1]-0.4330127018922193*Pxy_r[1]*uyMax_r[1]-0.4330127018922193*Pxy_c[1]*uyMax_r[1]+0.25*Pxy_r[0]*uyMax_r[1]-0.25*Pxy_c[0]*uyMax_r[1]-0.4330127018922193*energy_r[1]*uxMax_r[1]-0.4330127018922193*energy_c[1]*uxMax_r[1]-0.4330127018922193*Pxx_r[1]*uxMax_r[1]-0.4330127018922193*Pxx_c[1]*uxMax_r[1]+0.25*energy_r[0]*uxMax_r[1]-0.25*energy_c[0]*uxMax_r[1]+0.25*Pxx_r[0]*uxMax_r[1]-0.25*Pxx_c[0]*uxMax_r[1]; 

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

  outenergy[0] += (-0.7071067811865475*heat_flux_r[0]*dx1)+0.7071067811865475*heat_flux_l[0]*dx1-0.7071067811865475*Ghat_energy_r[0]*dx1+0.7071067811865475*Ghat_energy_l[0]*dx1; 
  outenergy[1] += (-1.224744871391589*heat_flux_r[0]*dx1)-1.224744871391589*heat_flux_l[0]*dx1-1.224744871391589*Ghat_energy_r[0]*dx1-1.224744871391589*Ghat_energy_l[0]*dx1; 
  outenergy[2] += (-0.7071067811865475*heat_flux_r[1]*dx1)+0.7071067811865475*heat_flux_l[1]*dx1-0.7071067811865475*Ghat_energy_r[1]*dx1+0.7071067811865475*Ghat_energy_l[1]*dx1; 
  outenergy[3] += (-1.224744871391589*heat_flux_r[1]*dx1)-1.224744871391589*heat_flux_l[1]*dx1-1.224744871391589*Ghat_energy_r[1]*dx1-1.224744871391589*Ghat_energy_l[1]*dx1; 

} 
