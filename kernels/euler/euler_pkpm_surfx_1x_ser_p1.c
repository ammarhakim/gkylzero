#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p1_surfx1_eval_quad.h> 
GKYL_CU_DH void euler_pkpm_surfx_1x_ser_p1(const double *w, const double *dxv, 
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
  const double *energy_l = &statevecl[6]; 
  const double *energy_c = &statevecc[6]; 
  const double *energy_r = &statevecr[6]; 

  const double *ux_l = &u_il[0]; 
  const double *ux_c = &u_ic[0]; 
  const double *ux_r = &u_ir[0]; 

  const double *uy_l = &u_il[2]; 
  const double *uy_c = &u_ic[2]; 
  const double *uy_r = &u_ir[2]; 

  const double *uz_l = &u_il[4]; 
  const double *uz_c = &u_ic[4]; 
  const double *uz_r = &u_ir[4]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxx_r = &p_ijr[0]; 

  const double *Pxy_l = &p_ijl[2]; 
  const double *Pxy_c = &p_ijc[2]; 
  const double *Pxy_r = &p_ijr[2]; 

  const double *Pxz_l = &p_ijl[4]; 
  const double *Pxz_c = &p_ijc[4]; 
  const double *Pxz_r = &p_ijr[4]; 

  const double *rho_flux_l = &vlasov_pkpm_surf_momsl[0]; 
  const double *heat_flux_l = &vlasov_pkpm_surf_momsl[1]; 
  const double *rho_flux_r = &vlasov_pkpm_surf_momsr[0]; 
  const double *heat_flux_r = &vlasov_pkpm_surf_momsr[1]; 
  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[2]; 
  double *outrhou2 = &out[4]; 
  double *outenergy = &out[6]; 

  double ux_l_r = ser_1x_p1_surfx1_eval_quad_node_0_r(ux_l); 
  double ux_c_l = ser_1x_p1_surfx1_eval_quad_node_0_l(ux_c); 
  double ux_c_r = ser_1x_p1_surfx1_eval_quad_node_0_r(ux_c); 
  double ux_r_l = ser_1x_p1_surfx1_eval_quad_node_0_l(ux_r); 

  double uy_l_r = ser_1x_p1_surfx1_eval_quad_node_0_r(uy_l); 
  double uy_c_l = ser_1x_p1_surfx1_eval_quad_node_0_l(uy_c); 
  double uy_c_r = ser_1x_p1_surfx1_eval_quad_node_0_r(uy_c); 
  double uy_r_l = ser_1x_p1_surfx1_eval_quad_node_0_l(uy_r); 

  double uz_l_r = ser_1x_p1_surfx1_eval_quad_node_0_r(uz_l); 
  double uz_c_l = ser_1x_p1_surfx1_eval_quad_node_0_l(uz_c); 
  double uz_c_r = ser_1x_p1_surfx1_eval_quad_node_0_r(uz_c); 
  double uz_r_l = ser_1x_p1_surfx1_eval_quad_node_0_l(uz_r); 

  double ux_max_l = fmax(fabs(ux_l_r), fabs(ux_c_l)); 
  double ux_max_r = fmax(fabs(ux_c_r), fabs(ux_r_l)); 
  double uy_max_l = fmax(fabs(uy_l_r), fabs(uy_c_l)); 
  double uy_max_r = fmax(fabs(uy_c_r), fabs(uy_r_l)); 
  double uz_max_l = fmax(fabs(uz_l_r), fabs(uz_c_l)); 
  double uz_max_r = fmax(fabs(uz_c_r), fabs(uz_r_l)); 

  double Ghat_rhoux_l = 0.6123724356957945*rho_flux_l[0]*ux_l[1]-0.6123724356957945*rho_flux_l[0]*ux_c[1]+0.6123724356957945*Pxx_l[1]-0.6123724356957945*Pxx_c[1]+0.3535533905932737*rho_flux_l[0]*ux_l[0]+0.3535533905932737*rho_flux_l[0]*ux_c[0]+0.3535533905932737*Pxx_l[0]+0.3535533905932737*Pxx_c[0]; 
  double Ghat_rhoux_r = (-0.6123724356957945*rho_flux_r[0]*ux_r[1])+0.6123724356957945*rho_flux_r[0]*ux_c[1]-0.6123724356957945*Pxx_r[1]+0.6123724356957945*Pxx_c[1]+0.3535533905932737*rho_flux_r[0]*ux_r[0]+0.3535533905932737*rho_flux_r[0]*ux_c[0]+0.3535533905932737*Pxx_r[0]+0.3535533905932737*Pxx_c[0]; 
  double Ghat_rhouy_l = 0.6123724356957945*rho_flux_l[0]*uy_l[1]-0.6123724356957945*rho_flux_l[0]*uy_c[1]+0.6123724356957945*Pxy_l[1]-0.6123724356957945*Pxy_c[1]+0.3535533905932737*rho_flux_l[0]*uy_l[0]+0.3535533905932737*rho_flux_l[0]*uy_c[0]+0.3535533905932737*Pxy_l[0]+0.3535533905932737*Pxy_c[0]; 
  double Ghat_rhouy_r = (-0.6123724356957945*rho_flux_r[0]*uy_r[1])+0.6123724356957945*rho_flux_r[0]*uy_c[1]-0.6123724356957945*Pxy_r[1]+0.6123724356957945*Pxy_c[1]+0.3535533905932737*rho_flux_r[0]*uy_r[0]+0.3535533905932737*rho_flux_r[0]*uy_c[0]+0.3535533905932737*Pxy_r[0]+0.3535533905932737*Pxy_c[0]; 
  double Ghat_rhouz_l = 0.6123724356957945*rho_flux_l[0]*uz_l[1]-0.6123724356957945*rho_flux_l[0]*uz_c[1]+0.6123724356957945*Pxz_l[1]-0.6123724356957945*Pxz_c[1]+0.3535533905932737*rho_flux_l[0]*uz_l[0]+0.3535533905932737*rho_flux_l[0]*uz_c[0]+0.3535533905932737*Pxz_l[0]+0.3535533905932737*Pxz_c[0]; 
  double Ghat_rhouz_r = (-0.6123724356957945*rho_flux_r[0]*uz_r[1])+0.6123724356957945*rho_flux_r[0]*uz_c[1]-0.6123724356957945*Pxz_r[1]+0.6123724356957945*Pxz_c[1]+0.3535533905932737*rho_flux_r[0]*uz_r[0]+0.3535533905932737*rho_flux_r[0]*uz_c[0]+0.3535533905932737*Pxz_r[0]+0.3535533905932737*Pxz_c[0]; 
  double Ghat_energy_l = (-0.6123724356957945*Pxz_l[1]*uz_max_l)-0.6123724356957945*Pxz_c[1]*uz_max_l-0.3535533905932737*Pxz_l[0]*uz_max_l+0.3535533905932737*Pxz_c[0]*uz_max_l-0.6123724356957945*Pxy_l[1]*uy_max_l-0.6123724356957945*Pxy_c[1]*uy_max_l-0.3535533905932737*Pxy_l[0]*uy_max_l+0.3535533905932737*Pxy_c[0]*uy_max_l-0.6123724356957945*energy_l[1]*ux_max_l-0.6123724356957945*energy_c[1]*ux_max_l-0.6123724356957945*Pxx_l[1]*ux_max_l-0.6123724356957945*Pxx_c[1]*ux_max_l-0.3535533905932737*energy_l[0]*ux_max_l+0.3535533905932737*energy_c[0]*ux_max_l-0.3535533905932737*Pxx_l[0]*ux_max_l+0.3535533905932737*Pxx_c[0]*ux_max_l+0.75*Pxz_l[1]*uz_l[1]+0.4330127018922193*Pxz_l[0]*uz_l[1]+0.75*Pxz_c[1]*uz_c[1]-0.4330127018922193*Pxz_c[0]*uz_c[1]+0.75*Pxy_l[1]*uy_l[1]+0.4330127018922193*Pxy_l[0]*uy_l[1]+0.75*Pxy_c[1]*uy_c[1]-0.4330127018922193*Pxy_c[0]*uy_c[1]+0.75*energy_l[1]*ux_l[1]+0.75*Pxx_l[1]*ux_l[1]+0.4330127018922193*energy_l[0]*ux_l[1]+0.4330127018922193*Pxx_l[0]*ux_l[1]+0.75*energy_c[1]*ux_c[1]+0.75*Pxx_c[1]*ux_c[1]-0.4330127018922193*energy_c[0]*ux_c[1]-0.4330127018922193*Pxx_c[0]*ux_c[1]+0.4330127018922193*ux_l[0]*energy_l[1]-0.4330127018922193*ux_c[0]*energy_c[1]+0.4330127018922193*uz_l[0]*Pxz_l[1]-0.4330127018922193*uz_c[0]*Pxz_c[1]+0.4330127018922193*uy_l[0]*Pxy_l[1]-0.4330127018922193*uy_c[0]*Pxy_c[1]+0.4330127018922193*ux_l[0]*Pxx_l[1]-0.4330127018922193*ux_c[0]*Pxx_c[1]+0.25*Pxz_l[0]*uz_l[0]+0.25*Pxz_c[0]*uz_c[0]+0.25*Pxy_l[0]*uy_l[0]+0.25*Pxy_c[0]*uy_c[0]+0.25*energy_l[0]*ux_l[0]+0.25*Pxx_l[0]*ux_l[0]+0.25*energy_c[0]*ux_c[0]+0.25*Pxx_c[0]*ux_c[0]; 
  double Ghat_energy_r = (-0.6123724356957945*Pxz_r[1]*uz_max_r)-0.6123724356957945*Pxz_c[1]*uz_max_r+0.3535533905932737*Pxz_r[0]*uz_max_r-0.3535533905932737*Pxz_c[0]*uz_max_r-0.6123724356957945*Pxy_r[1]*uy_max_r-0.6123724356957945*Pxy_c[1]*uy_max_r+0.3535533905932737*Pxy_r[0]*uy_max_r-0.3535533905932737*Pxy_c[0]*uy_max_r-0.6123724356957945*energy_r[1]*ux_max_r-0.6123724356957945*energy_c[1]*ux_max_r-0.6123724356957945*Pxx_r[1]*ux_max_r-0.6123724356957945*Pxx_c[1]*ux_max_r+0.3535533905932737*energy_r[0]*ux_max_r-0.3535533905932737*energy_c[0]*ux_max_r+0.3535533905932737*Pxx_r[0]*ux_max_r-0.3535533905932737*Pxx_c[0]*ux_max_r+0.75*Pxz_r[1]*uz_r[1]-0.4330127018922193*Pxz_r[0]*uz_r[1]+0.75*Pxz_c[1]*uz_c[1]+0.4330127018922193*Pxz_c[0]*uz_c[1]+0.75*Pxy_r[1]*uy_r[1]-0.4330127018922193*Pxy_r[0]*uy_r[1]+0.75*Pxy_c[1]*uy_c[1]+0.4330127018922193*Pxy_c[0]*uy_c[1]+0.75*energy_r[1]*ux_r[1]+0.75*Pxx_r[1]*ux_r[1]-0.4330127018922193*energy_r[0]*ux_r[1]-0.4330127018922193*Pxx_r[0]*ux_r[1]+0.75*energy_c[1]*ux_c[1]+0.75*Pxx_c[1]*ux_c[1]+0.4330127018922193*energy_c[0]*ux_c[1]+0.4330127018922193*Pxx_c[0]*ux_c[1]-0.4330127018922193*ux_r[0]*energy_r[1]+0.4330127018922193*ux_c[0]*energy_c[1]-0.4330127018922193*uz_r[0]*Pxz_r[1]+0.4330127018922193*uz_c[0]*Pxz_c[1]-0.4330127018922193*uy_r[0]*Pxy_r[1]+0.4330127018922193*uy_c[0]*Pxy_c[1]-0.4330127018922193*ux_r[0]*Pxx_r[1]+0.4330127018922193*ux_c[0]*Pxx_c[1]+0.25*Pxz_r[0]*uz_r[0]+0.25*Pxz_c[0]*uz_c[0]+0.25*Pxy_r[0]*uy_r[0]+0.25*Pxy_c[0]*uy_c[0]+0.25*energy_r[0]*ux_r[0]+0.25*Pxx_r[0]*ux_r[0]+0.25*energy_c[0]*ux_c[0]+0.25*Pxx_c[0]*ux_c[0]; 

  outrhou0[0] += 0.7071067811865475*Ghat_rhoux_l*dx1-0.7071067811865475*Ghat_rhoux_r*dx1; 
  outrhou0[1] += (-1.224744871391589*Ghat_rhoux_r*dx1)-1.224744871391589*Ghat_rhoux_l*dx1; 

  outrhou1[0] += 0.7071067811865475*Ghat_rhouy_l*dx1-0.7071067811865475*Ghat_rhouy_r*dx1; 
  outrhou1[1] += (-1.224744871391589*Ghat_rhouy_r*dx1)-1.224744871391589*Ghat_rhouy_l*dx1; 

  outrhou2[0] += 0.7071067811865475*Ghat_rhouz_l*dx1-0.7071067811865475*Ghat_rhouz_r*dx1; 
  outrhou2[1] += (-1.224744871391589*Ghat_rhouz_r*dx1)-1.224744871391589*Ghat_rhouz_l*dx1; 

  outenergy[0] += (-0.7071067811865475*Ghat_energy_r*dx1)+0.7071067811865475*Ghat_energy_l*dx1-0.7071067811865475*heat_flux_r[0]*dx1+0.7071067811865475*heat_flux_l[0]*dx1; 
  outenergy[1] += (-1.224744871391589*Ghat_energy_r*dx1)-1.224744871391589*Ghat_energy_l*dx1-1.224744871391589*heat_flux_r[0]*dx1-1.224744871391589*heat_flux_l[0]*dx1; 

} 
