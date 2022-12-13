#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void euler_pkpm_prim_vars_1x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, 
  double* u_i, double* p_ij, double* T_ij) 
{ 
  // bvar:             Magnetic field unit vector and tensor.
  // vlasov_pkpm_moms: [rho, p_parallel, q_parallel], Moments computed from kinetic equation in pkpm model.
  // statevec:         [rho ux, rho uy, rho uz, p_perp], Fluid input state vector.
  // u_i:              Output flow velocity [ux, uy, uz].
  // p_ij:             Output pressure tensor, p_ij = (p_parallel - p_perp)bb + p_perp I.
  // T_ij:             Output Temperature tensor for penalization T_ij = 3.0*p_ij/rho.

  const double *rhoux = &statevec[0]; 
  const double *rhouy = &statevec[2]; 
  const double *rhouz = &statevec[4]; 
  const double *p_perp = &statevec[6]; 
  // Parallel pressure is first component of pkpm moment array and unit tensor are last six components of bvar array.
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_parallel = &vlasov_pkpm_moms[2]; 
  const double *bxbx = &bvar[6]; 
  const double *bxby = &bvar[8]; 
  const double *bxbz = &bvar[10]; 
  const double *byby = &bvar[12]; 
  const double *bybz = &bvar[14]; 
  const double *bzbz = &bvar[16]; 

  double *ux = &u_i[0]; 
  double *uy = &u_i[2]; 
  double *uz = &u_i[4]; 

  double *Pxx = &p_ij[0]; 
  double *Pxy = &p_ij[2]; 
  double *Pxz = &p_ij[4]; 
  double *Pyy = &p_ij[6]; 
  double *Pyz = &p_ij[8]; 
  double *Pzz = &p_ij[10]; 

  double *Txx = &T_ij[0]; 
  double *Tyy = &T_ij[6]; 
  double *Tzz = &T_ij[10]; 

  double rho_inv[2] = {0.0}; 
  ser_1x_p1_inv(rho, rho_inv); 
  ux[0] = 0.7071067811865475*rho_inv[1]*rhoux[1]+0.7071067811865475*rho_inv[0]*rhoux[0]; 
  ux[1] = 0.7071067811865475*rho_inv[0]*rhoux[1]+0.7071067811865475*rhoux[0]*rho_inv[1]; 

  uy[0] = 0.7071067811865475*rho_inv[1]*rhouy[1]+0.7071067811865475*rho_inv[0]*rhouy[0]; 
  uy[1] = 0.7071067811865475*rho_inv[0]*rhouy[1]+0.7071067811865475*rhouy[0]*rho_inv[1]; 

  uz[0] = 0.7071067811865475*rho_inv[1]*rhouz[1]+0.7071067811865475*rho_inv[0]*rhouz[0]; 
  uz[1] = 0.7071067811865475*rho_inv[0]*rhouz[1]+0.7071067811865475*rhouz[0]*rho_inv[1]; 

  Pxx[0] = (-0.7071067811865475*bxbx[1]*p_perp[1])+0.7071067811865475*bxbx[1]*p_parallel[1]-0.7071067811865475*bxbx[0]*p_perp[0]+p_perp[0]+0.7071067811865475*bxbx[0]*p_parallel[0]; 
  Pxx[1] = (-0.7071067811865475*bxbx[0]*p_perp[1])+p_perp[1]+0.7071067811865475*bxbx[0]*p_parallel[1]-0.7071067811865475*p_perp[0]*bxbx[1]+0.7071067811865475*p_parallel[0]*bxbx[1]; 

  Pxy[0] = (-0.7071067811865475*bxby[1]*p_perp[1])+0.7071067811865475*bxby[1]*p_parallel[1]-0.7071067811865475*bxby[0]*p_perp[0]+0.7071067811865475*bxby[0]*p_parallel[0]; 
  Pxy[1] = (-0.7071067811865475*bxby[0]*p_perp[1])+0.7071067811865475*bxby[0]*p_parallel[1]-0.7071067811865475*p_perp[0]*bxby[1]+0.7071067811865475*p_parallel[0]*bxby[1]; 

  Pxz[0] = (-0.7071067811865475*bxbz[1]*p_perp[1])+0.7071067811865475*bxbz[1]*p_parallel[1]-0.7071067811865475*bxbz[0]*p_perp[0]+0.7071067811865475*bxbz[0]*p_parallel[0]; 
  Pxz[1] = (-0.7071067811865475*bxbz[0]*p_perp[1])+0.7071067811865475*bxbz[0]*p_parallel[1]-0.7071067811865475*p_perp[0]*bxbz[1]+0.7071067811865475*p_parallel[0]*bxbz[1]; 

  Pyy[0] = (-0.7071067811865475*byby[1]*p_perp[1])+0.7071067811865475*byby[1]*p_parallel[1]-0.7071067811865475*byby[0]*p_perp[0]+p_perp[0]+0.7071067811865475*byby[0]*p_parallel[0]; 
  Pyy[1] = (-0.7071067811865475*byby[0]*p_perp[1])+p_perp[1]+0.7071067811865475*byby[0]*p_parallel[1]-0.7071067811865475*p_perp[0]*byby[1]+0.7071067811865475*p_parallel[0]*byby[1]; 

  Pyz[0] = (-0.7071067811865475*bybz[1]*p_perp[1])+0.7071067811865475*bybz[1]*p_parallel[1]-0.7071067811865475*bybz[0]*p_perp[0]+0.7071067811865475*bybz[0]*p_parallel[0]; 
  Pyz[1] = (-0.7071067811865475*bybz[0]*p_perp[1])+0.7071067811865475*bybz[0]*p_parallel[1]-0.7071067811865475*p_perp[0]*bybz[1]+0.7071067811865475*p_parallel[0]*bybz[1]; 

  Pzz[0] = (-0.7071067811865475*bzbz[1]*p_perp[1])+0.7071067811865475*bzbz[1]*p_parallel[1]-0.7071067811865475*bzbz[0]*p_perp[0]+p_perp[0]+0.7071067811865475*bzbz[0]*p_parallel[0]; 
  Pzz[1] = (-0.7071067811865475*bzbz[0]*p_perp[1])+p_perp[1]+0.7071067811865475*bzbz[0]*p_parallel[1]-0.7071067811865475*p_perp[0]*bzbz[1]+0.7071067811865475*p_parallel[0]*bzbz[1]; 
  Txx[0] = 2.121320343559642*Pxx[1]*rho_inv[1]+2.121320343559642*Pxx[0]*rho_inv[0]; 
  Txx[1] = 2.121320343559642*Pxx[0]*rho_inv[1]+2.121320343559642*rho_inv[0]*Pxx[1]; 

  Tyy[0] = 2.121320343559642*Pyy[1]*rho_inv[1]+2.121320343559642*Pyy[0]*rho_inv[0]; 
  Tyy[1] = 2.121320343559642*Pyy[0]*rho_inv[1]+2.121320343559642*rho_inv[0]*Pyy[1]; 

  Tzz[0] = 2.121320343559642*Pzz[1]*rho_inv[1]+2.121320343559642*Pzz[0]*rho_inv[0]; 
  Tzz[1] = 2.121320343559642*Pzz[0]*rho_inv[1]+2.121320343559642*rho_inv[0]*Pzz[1]; 
} 
