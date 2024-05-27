#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_io_1x_tensor_p2(const double *vlasov_pkpm_moms, 
  const double* pkpm_u, const double* p_ij, const double* pkpm_prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io) 
{ 
  // vlasov_pkpm_moms: Input [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model. 
  // pkpm_u:           Input [ux, uy, uz]. 
  // p_ij:             Input p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij. 
  // pkpm_prim:        Input [1/rho div(p_par b), T_perp/m, m/T_perp, 3*Txx/m, 3*Tyy/m, 3*Tzz/m]. 
  // pkpm_accel:       Input volume expansion of pkpm acceleration variables [T_perp/m*div(b), bb:grad(u), p_force, p_perp_source]. 
  // fluid_io:         Output fluid conserved variables. 
  // pkpm_vars_io:     Output pkpm variables (primitive and acceleration). 

  const double *rho = &vlasov_pkpm_moms[0]; 

  const double *ux = &pkpm_u[0]; 
  const double *uy = &pkpm_u[2]; 
  const double *uz = &pkpm_u[4]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[3]; 
  const double *Pxz = &p_ij[6]; 
  const double *Pyy = &p_ij[9]; 
  const double *Pyz = &p_ij[12]; 
  const double *Pzz = &p_ij[15]; 

  const double *pkpm_div_ppar = &pkpm_prim[0]; 
  const double *T_perp_over_m = &pkpm_prim[3]; 
  const double *T_perp_over_m_inv = &pkpm_prim[6]; 

  const double *p_perp_div_b = &pkpm_accel[0]; 
  const double *bb_grad_u = &pkpm_accel[3]; 

  double *fluid_io_rho = &fluid_io[0]; 
  double *fluid_io_rhoux = &fluid_io[3]; 
  double *fluid_io_rhouy = &fluid_io[6]; 
  double *fluid_io_rhouz = &fluid_io[9]; 
  double *fluid_io_Sxx = &fluid_io[12]; 
  double *fluid_io_Sxy = &fluid_io[15]; 
  double *fluid_io_Sxz = &fluid_io[18]; 
  double *fluid_io_Syy = &fluid_io[21]; 
  double *fluid_io_Syz = &fluid_io[24]; 
  double *fluid_io_Szz = &fluid_io[27]; 

  double *pkpm_vars_io_T_perp_over_m = &pkpm_vars_io[0]; 
  double *pkpm_vars_io_T_perp_over_m_inv = &pkpm_vars_io[3]; 
  double *pkpm_vars_io_pkpm_div_ppar = &pkpm_vars_io[6]; 
  double *pkpm_vars_io_p_perp_div_b = &pkpm_vars_io[9]; 
  double *pkpm_vars_io_bb_grad_u = &pkpm_vars_io[12]; 

  double rhoux[3] = {0.0}; 
  rhoux[0] = 0.7071067811865475*rho[1]*ux[1]+0.7071067811865475*rho[0]*ux[0]; 
  rhoux[1] = 0.6324555320336759*ux[1]*rho[2]+0.7071067811865475*rho[0]*ux[1]+0.7071067811865475*ux[0]*rho[1]; 
  rhoux[2] = 0.7071067811865475*ux[0]*rho[2]+0.6324555320336759*rho[1]*ux[1]; 
 
  double rhouy[3] = {0.0}; 
  rhouy[0] = 0.7071067811865475*rho[1]*uy[1]+0.7071067811865475*rho[0]*uy[0]; 
  rhouy[1] = 0.6324555320336759*uy[1]*rho[2]+0.7071067811865475*rho[0]*uy[1]+0.7071067811865475*uy[0]*rho[1]; 
  rhouy[2] = 0.7071067811865475*uy[0]*rho[2]+0.6324555320336759*rho[1]*uy[1]; 
 
  double rhouz[3] = {0.0}; 
  rhouz[0] = 0.7071067811865475*rho[1]*uz[1]+0.7071067811865475*rho[0]*uz[0]; 
  rhouz[1] = 0.6324555320336759*uz[1]*rho[2]+0.7071067811865475*rho[0]*uz[1]+0.7071067811865475*uz[0]*rho[1]; 
  rhouz[2] = 0.7071067811865475*uz[0]*rho[2]+0.6324555320336759*rho[1]*uz[1]; 

  double rhouxux[3] = {0.0}; 
  rhouxux[0] = 0.7071067811865475*rhoux[1]*ux[1]+0.7071067811865475*rhoux[0]*ux[0]; 
  rhouxux[1] = 0.6324555320336759*ux[1]*rhoux[2]+0.7071067811865475*rhoux[0]*ux[1]+0.7071067811865475*ux[0]*rhoux[1]; 
  rhouxux[2] = 0.7071067811865475*ux[0]*rhoux[2]+0.6324555320336759*rhoux[1]*ux[1]; 
 
  double rhouxuy[3] = {0.0}; 
  rhouxuy[0] = 0.7071067811865475*rhoux[1]*uy[1]+0.7071067811865475*rhoux[0]*uy[0]; 
  rhouxuy[1] = 0.6324555320336759*uy[1]*rhoux[2]+0.7071067811865475*rhoux[0]*uy[1]+0.7071067811865475*uy[0]*rhoux[1]; 
  rhouxuy[2] = 0.7071067811865475*uy[0]*rhoux[2]+0.6324555320336759*rhoux[1]*uy[1]; 
 
  double rhouxuz[3] = {0.0}; 
  rhouxuz[0] = 0.7071067811865475*rhoux[1]*uz[1]+0.7071067811865475*rhoux[0]*uz[0]; 
  rhouxuz[1] = 0.6324555320336759*uz[1]*rhoux[2]+0.7071067811865475*rhoux[0]*uz[1]+0.7071067811865475*uz[0]*rhoux[1]; 
  rhouxuz[2] = 0.7071067811865475*uz[0]*rhoux[2]+0.6324555320336759*rhoux[1]*uz[1]; 
 
  double rhouyuy[3] = {0.0}; 
  rhouyuy[0] = 0.7071067811865475*rhouy[1]*uy[1]+0.7071067811865475*rhouy[0]*uy[0]; 
  rhouyuy[1] = 0.6324555320336759*uy[1]*rhouy[2]+0.7071067811865475*rhouy[0]*uy[1]+0.7071067811865475*uy[0]*rhouy[1]; 
  rhouyuy[2] = 0.7071067811865475*uy[0]*rhouy[2]+0.6324555320336759*rhouy[1]*uy[1]; 
 
  double rhouyuz[3] = {0.0}; 
  rhouyuz[0] = 0.7071067811865475*rhouy[1]*uz[1]+0.7071067811865475*rhouy[0]*uz[0]; 
  rhouyuz[1] = 0.6324555320336759*uz[1]*rhouy[2]+0.7071067811865475*rhouy[0]*uz[1]+0.7071067811865475*uz[0]*rhouy[1]; 
  rhouyuz[2] = 0.7071067811865475*uz[0]*rhouy[2]+0.6324555320336759*rhouy[1]*uz[1]; 
 
  double rhouzuz[3] = {0.0}; 
  rhouzuz[0] = 0.7071067811865475*rhouz[1]*uz[1]+0.7071067811865475*rhouz[0]*uz[0]; 
  rhouzuz[1] = 0.6324555320336759*uz[1]*rhouz[2]+0.7071067811865475*rhouz[0]*uz[1]+0.7071067811865475*uz[0]*rhouz[1]; 
  rhouzuz[2] = 0.7071067811865475*uz[0]*rhouz[2]+0.6324555320336759*rhouz[1]*uz[1]; 
 
  fluid_io_rho[0] = rho[0]; 
  fluid_io_rhoux[0] = rhoux[0]; 
  fluid_io_rhouy[0] = rhouy[0]; 
  fluid_io_rhouz[0] = rhouz[0]; 
  fluid_io_Sxx[0] = Pxx[0] + rhouxux[0]; 
  fluid_io_Sxy[0] = Pxy[0] + rhouxuy[0]; 
  fluid_io_Sxz[0] = Pxz[0] + rhouxuz[0]; 
  fluid_io_Syy[0] = Pyy[0] + rhouyuy[0]; 
  fluid_io_Syz[0] = Pyz[0] + rhouyuz[0]; 
  fluid_io_Szz[0] = Pzz[0] + rhouzuz[0]; 
 
  pkpm_vars_io_T_perp_over_m[0] = T_perp_over_m[0]; 
  pkpm_vars_io_T_perp_over_m_inv[0] = T_perp_over_m_inv[0]; 
  pkpm_vars_io_pkpm_div_ppar[0] = pkpm_div_ppar[0]; 
  pkpm_vars_io_p_perp_div_b[0] = p_perp_div_b[0]; 
  pkpm_vars_io_bb_grad_u[0] = bb_grad_u[0]; 
 
  fluid_io_rho[1] = rho[1]; 
  fluid_io_rhoux[1] = rhoux[1]; 
  fluid_io_rhouy[1] = rhouy[1]; 
  fluid_io_rhouz[1] = rhouz[1]; 
  fluid_io_Sxx[1] = Pxx[1] + rhouxux[1]; 
  fluid_io_Sxy[1] = Pxy[1] + rhouxuy[1]; 
  fluid_io_Sxz[1] = Pxz[1] + rhouxuz[1]; 
  fluid_io_Syy[1] = Pyy[1] + rhouyuy[1]; 
  fluid_io_Syz[1] = Pyz[1] + rhouyuz[1]; 
  fluid_io_Szz[1] = Pzz[1] + rhouzuz[1]; 
 
  pkpm_vars_io_T_perp_over_m[1] = T_perp_over_m[1]; 
  pkpm_vars_io_T_perp_over_m_inv[1] = T_perp_over_m_inv[1]; 
  pkpm_vars_io_pkpm_div_ppar[1] = pkpm_div_ppar[1]; 
  pkpm_vars_io_p_perp_div_b[1] = p_perp_div_b[1]; 
  pkpm_vars_io_bb_grad_u[1] = bb_grad_u[1]; 
 
  fluid_io_rho[2] = rho[2]; 
  fluid_io_rhoux[2] = rhoux[2]; 
  fluid_io_rhouy[2] = rhouy[2]; 
  fluid_io_rhouz[2] = rhouz[2]; 
  fluid_io_Sxx[2] = Pxx[2] + rhouxux[2]; 
  fluid_io_Sxy[2] = Pxy[2] + rhouxuy[2]; 
  fluid_io_Sxz[2] = Pxz[2] + rhouxuz[2]; 
  fluid_io_Syy[2] = Pyy[2] + rhouyuy[2]; 
  fluid_io_Syz[2] = Pyz[2] + rhouyuz[2]; 
  fluid_io_Szz[2] = Pzz[2] + rhouzuz[2]; 
 
  pkpm_vars_io_T_perp_over_m[2] = T_perp_over_m[2]; 
  pkpm_vars_io_T_perp_over_m_inv[2] = T_perp_over_m_inv[2]; 
  pkpm_vars_io_pkpm_div_ppar[2] = pkpm_div_ppar[2]; 
  pkpm_vars_io_p_perp_div_b[2] = p_perp_div_b[2]; 
  pkpm_vars_io_bb_grad_u[2] = bb_grad_u[2]; 
 
} 
