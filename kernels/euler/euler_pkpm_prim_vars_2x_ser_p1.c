#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void euler_pkpm_prim_vars_2x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, 
  double* u_i, double* p_ij, double* T_ij) 
{ 
  // bvar:             Magnetic field unit vector and tensor.
  // vlasov_pkpm_moms: [rho, p_parallel, q_parallel], Moments computed from kinetic equation in pkpm model.
  // statevec:         [rho ux, rho uy, rho uz, p_perp], Fluid input state vector.
  // u_i:              Output flow velocity [ux, uy, uz].
  // p_ij:             Output pressure tensor, p_ij = (p_parallel - p_perp)bb + p_perp I.
  // T_ij:             Output Temperature tensor for penalization T_ij = 3.0*p_ij/rho.

  const double *rhoux = &statevec[0]; 
  const double *rhouy = &statevec[4]; 
  const double *rhouz = &statevec[8]; 
  const double *p_perp = &statevec[12]; 
  // Parallel pressure is first component of pkpm moment array and unit tensor are last six components of bvar array.
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_parallel = &vlasov_pkpm_moms[4]; 
  const double *bxbx = &bvar[12]; 
  const double *bxby = &bvar[16]; 
  const double *bxbz = &bvar[20]; 
  const double *byby = &bvar[24]; 
  const double *bybz = &bvar[28]; 
  const double *bzbz = &bvar[32]; 

  double *ux = &u_i[0]; 
  double *uy = &u_i[4]; 
  double *uz = &u_i[8]; 

  double *Pxx = &p_ij[0]; 
  double *Pxy = &p_ij[4]; 
  double *Pxz = &p_ij[8]; 
  double *Pyy = &p_ij[12]; 
  double *Pyz = &p_ij[16]; 
  double *Pzz = &p_ij[20]; 

  double *Txx = &T_ij[0]; 
  double *Tyy = &T_ij[12]; 
  double *Tzz = &T_ij[20]; 

  double rho_inv[4] = {0.0}; 
  ser_2x_p1_inv(rho, rho_inv); 
  ux[0] = 0.5*rho_inv[3]*rhoux[3]+0.5*rho_inv[2]*rhoux[2]+0.5*rho_inv[1]*rhoux[1]+0.5*rho_inv[0]*rhoux[0]; 
  ux[1] = 0.5*rho_inv[2]*rhoux[3]+0.5*rhoux[2]*rho_inv[3]+0.5*rho_inv[0]*rhoux[1]+0.5*rhoux[0]*rho_inv[1]; 
  ux[2] = 0.5*rho_inv[1]*rhoux[3]+0.5*rhoux[1]*rho_inv[3]+0.5*rho_inv[0]*rhoux[2]+0.5*rhoux[0]*rho_inv[2]; 
  ux[3] = 0.5*rho_inv[0]*rhoux[3]+0.5*rhoux[0]*rho_inv[3]+0.5*rho_inv[1]*rhoux[2]+0.5*rhoux[1]*rho_inv[2]; 

  uy[0] = 0.5*rho_inv[3]*rhouy[3]+0.5*rho_inv[2]*rhouy[2]+0.5*rho_inv[1]*rhouy[1]+0.5*rho_inv[0]*rhouy[0]; 
  uy[1] = 0.5*rho_inv[2]*rhouy[3]+0.5*rhouy[2]*rho_inv[3]+0.5*rho_inv[0]*rhouy[1]+0.5*rhouy[0]*rho_inv[1]; 
  uy[2] = 0.5*rho_inv[1]*rhouy[3]+0.5*rhouy[1]*rho_inv[3]+0.5*rho_inv[0]*rhouy[2]+0.5*rhouy[0]*rho_inv[2]; 
  uy[3] = 0.5*rho_inv[0]*rhouy[3]+0.5*rhouy[0]*rho_inv[3]+0.5*rho_inv[1]*rhouy[2]+0.5*rhouy[1]*rho_inv[2]; 

  uz[0] = 0.5*rho_inv[3]*rhouz[3]+0.5*rho_inv[2]*rhouz[2]+0.5*rho_inv[1]*rhouz[1]+0.5*rho_inv[0]*rhouz[0]; 
  uz[1] = 0.5*rho_inv[2]*rhouz[3]+0.5*rhouz[2]*rho_inv[3]+0.5*rho_inv[0]*rhouz[1]+0.5*rhouz[0]*rho_inv[1]; 
  uz[2] = 0.5*rho_inv[1]*rhouz[3]+0.5*rhouz[1]*rho_inv[3]+0.5*rho_inv[0]*rhouz[2]+0.5*rhouz[0]*rho_inv[2]; 
  uz[3] = 0.5*rho_inv[0]*rhouz[3]+0.5*rhouz[0]*rho_inv[3]+0.5*rho_inv[1]*rhouz[2]+0.5*rhouz[1]*rho_inv[2]; 

  Pxx[0] = (-0.5*bxbx[3]*p_perp[3])+0.5*bxbx[3]*p_parallel[3]-0.5*bxbx[2]*p_perp[2]+0.5*bxbx[2]*p_parallel[2]-0.5*bxbx[1]*p_perp[1]+0.5*bxbx[1]*p_parallel[1]-0.5*bxbx[0]*p_perp[0]+p_perp[0]+0.5*bxbx[0]*p_parallel[0]; 
  Pxx[1] = (-0.5*bxbx[2]*p_perp[3])+0.5*bxbx[2]*p_parallel[3]-0.5*p_perp[2]*bxbx[3]+0.5*p_parallel[2]*bxbx[3]-0.5*bxbx[0]*p_perp[1]+p_perp[1]+0.5*bxbx[0]*p_parallel[1]-0.5*p_perp[0]*bxbx[1]+0.5*p_parallel[0]*bxbx[1]; 
  Pxx[2] = (-0.5*bxbx[1]*p_perp[3])+0.5*bxbx[1]*p_parallel[3]-0.5*p_perp[1]*bxbx[3]+0.5*p_parallel[1]*bxbx[3]-0.5*bxbx[0]*p_perp[2]+p_perp[2]+0.5*bxbx[0]*p_parallel[2]-0.5*p_perp[0]*bxbx[2]+0.5*p_parallel[0]*bxbx[2]; 
  Pxx[3] = (-0.5*bxbx[0]*p_perp[3])+p_perp[3]+0.5*bxbx[0]*p_parallel[3]-0.5*p_perp[0]*bxbx[3]+0.5*p_parallel[0]*bxbx[3]-0.5*bxbx[1]*p_perp[2]+0.5*bxbx[1]*p_parallel[2]-0.5*p_perp[1]*bxbx[2]+0.5*p_parallel[1]*bxbx[2]; 

  Pxy[0] = (-0.5*bxby[3]*p_perp[3])+0.5*bxby[3]*p_parallel[3]-0.5*bxby[2]*p_perp[2]+0.5*bxby[2]*p_parallel[2]-0.5*bxby[1]*p_perp[1]+0.5*bxby[1]*p_parallel[1]-0.5*bxby[0]*p_perp[0]+0.5*bxby[0]*p_parallel[0]; 
  Pxy[1] = (-0.5*bxby[2]*p_perp[3])+0.5*bxby[2]*p_parallel[3]-0.5*p_perp[2]*bxby[3]+0.5*p_parallel[2]*bxby[3]-0.5*bxby[0]*p_perp[1]+0.5*bxby[0]*p_parallel[1]-0.5*p_perp[0]*bxby[1]+0.5*p_parallel[0]*bxby[1]; 
  Pxy[2] = (-0.5*bxby[1]*p_perp[3])+0.5*bxby[1]*p_parallel[3]-0.5*p_perp[1]*bxby[3]+0.5*p_parallel[1]*bxby[3]-0.5*bxby[0]*p_perp[2]+0.5*bxby[0]*p_parallel[2]-0.5*p_perp[0]*bxby[2]+0.5*p_parallel[0]*bxby[2]; 
  Pxy[3] = (-0.5*bxby[0]*p_perp[3])+0.5*bxby[0]*p_parallel[3]-0.5*p_perp[0]*bxby[3]+0.5*p_parallel[0]*bxby[3]-0.5*bxby[1]*p_perp[2]+0.5*bxby[1]*p_parallel[2]-0.5*p_perp[1]*bxby[2]+0.5*p_parallel[1]*bxby[2]; 

  Pxz[0] = (-0.5*bxbz[3]*p_perp[3])+0.5*bxbz[3]*p_parallel[3]-0.5*bxbz[2]*p_perp[2]+0.5*bxbz[2]*p_parallel[2]-0.5*bxbz[1]*p_perp[1]+0.5*bxbz[1]*p_parallel[1]-0.5*bxbz[0]*p_perp[0]+0.5*bxbz[0]*p_parallel[0]; 
  Pxz[1] = (-0.5*bxbz[2]*p_perp[3])+0.5*bxbz[2]*p_parallel[3]-0.5*p_perp[2]*bxbz[3]+0.5*p_parallel[2]*bxbz[3]-0.5*bxbz[0]*p_perp[1]+0.5*bxbz[0]*p_parallel[1]-0.5*p_perp[0]*bxbz[1]+0.5*p_parallel[0]*bxbz[1]; 
  Pxz[2] = (-0.5*bxbz[1]*p_perp[3])+0.5*bxbz[1]*p_parallel[3]-0.5*p_perp[1]*bxbz[3]+0.5*p_parallel[1]*bxbz[3]-0.5*bxbz[0]*p_perp[2]+0.5*bxbz[0]*p_parallel[2]-0.5*p_perp[0]*bxbz[2]+0.5*p_parallel[0]*bxbz[2]; 
  Pxz[3] = (-0.5*bxbz[0]*p_perp[3])+0.5*bxbz[0]*p_parallel[3]-0.5*p_perp[0]*bxbz[3]+0.5*p_parallel[0]*bxbz[3]-0.5*bxbz[1]*p_perp[2]+0.5*bxbz[1]*p_parallel[2]-0.5*p_perp[1]*bxbz[2]+0.5*p_parallel[1]*bxbz[2]; 

  Pyy[0] = (-0.5*byby[3]*p_perp[3])+0.5*byby[3]*p_parallel[3]-0.5*byby[2]*p_perp[2]+0.5*byby[2]*p_parallel[2]-0.5*byby[1]*p_perp[1]+0.5*byby[1]*p_parallel[1]-0.5*byby[0]*p_perp[0]+p_perp[0]+0.5*byby[0]*p_parallel[0]; 
  Pyy[1] = (-0.5*byby[2]*p_perp[3])+0.5*byby[2]*p_parallel[3]-0.5*p_perp[2]*byby[3]+0.5*p_parallel[2]*byby[3]-0.5*byby[0]*p_perp[1]+p_perp[1]+0.5*byby[0]*p_parallel[1]-0.5*p_perp[0]*byby[1]+0.5*p_parallel[0]*byby[1]; 
  Pyy[2] = (-0.5*byby[1]*p_perp[3])+0.5*byby[1]*p_parallel[3]-0.5*p_perp[1]*byby[3]+0.5*p_parallel[1]*byby[3]-0.5*byby[0]*p_perp[2]+p_perp[2]+0.5*byby[0]*p_parallel[2]-0.5*p_perp[0]*byby[2]+0.5*p_parallel[0]*byby[2]; 
  Pyy[3] = (-0.5*byby[0]*p_perp[3])+p_perp[3]+0.5*byby[0]*p_parallel[3]-0.5*p_perp[0]*byby[3]+0.5*p_parallel[0]*byby[3]-0.5*byby[1]*p_perp[2]+0.5*byby[1]*p_parallel[2]-0.5*p_perp[1]*byby[2]+0.5*p_parallel[1]*byby[2]; 

  Pyz[0] = (-0.5*bybz[3]*p_perp[3])+0.5*bybz[3]*p_parallel[3]-0.5*bybz[2]*p_perp[2]+0.5*bybz[2]*p_parallel[2]-0.5*bybz[1]*p_perp[1]+0.5*bybz[1]*p_parallel[1]-0.5*bybz[0]*p_perp[0]+0.5*bybz[0]*p_parallel[0]; 
  Pyz[1] = (-0.5*bybz[2]*p_perp[3])+0.5*bybz[2]*p_parallel[3]-0.5*p_perp[2]*bybz[3]+0.5*p_parallel[2]*bybz[3]-0.5*bybz[0]*p_perp[1]+0.5*bybz[0]*p_parallel[1]-0.5*p_perp[0]*bybz[1]+0.5*p_parallel[0]*bybz[1]; 
  Pyz[2] = (-0.5*bybz[1]*p_perp[3])+0.5*bybz[1]*p_parallel[3]-0.5*p_perp[1]*bybz[3]+0.5*p_parallel[1]*bybz[3]-0.5*bybz[0]*p_perp[2]+0.5*bybz[0]*p_parallel[2]-0.5*p_perp[0]*bybz[2]+0.5*p_parallel[0]*bybz[2]; 
  Pyz[3] = (-0.5*bybz[0]*p_perp[3])+0.5*bybz[0]*p_parallel[3]-0.5*p_perp[0]*bybz[3]+0.5*p_parallel[0]*bybz[3]-0.5*bybz[1]*p_perp[2]+0.5*bybz[1]*p_parallel[2]-0.5*p_perp[1]*bybz[2]+0.5*p_parallel[1]*bybz[2]; 

  Pzz[0] = (-0.5*bzbz[3]*p_perp[3])+0.5*bzbz[3]*p_parallel[3]-0.5*bzbz[2]*p_perp[2]+0.5*bzbz[2]*p_parallel[2]-0.5*bzbz[1]*p_perp[1]+0.5*bzbz[1]*p_parallel[1]-0.5*bzbz[0]*p_perp[0]+p_perp[0]+0.5*bzbz[0]*p_parallel[0]; 
  Pzz[1] = (-0.5*bzbz[2]*p_perp[3])+0.5*bzbz[2]*p_parallel[3]-0.5*p_perp[2]*bzbz[3]+0.5*p_parallel[2]*bzbz[3]-0.5*bzbz[0]*p_perp[1]+p_perp[1]+0.5*bzbz[0]*p_parallel[1]-0.5*p_perp[0]*bzbz[1]+0.5*p_parallel[0]*bzbz[1]; 
  Pzz[2] = (-0.5*bzbz[1]*p_perp[3])+0.5*bzbz[1]*p_parallel[3]-0.5*p_perp[1]*bzbz[3]+0.5*p_parallel[1]*bzbz[3]-0.5*bzbz[0]*p_perp[2]+p_perp[2]+0.5*bzbz[0]*p_parallel[2]-0.5*p_perp[0]*bzbz[2]+0.5*p_parallel[0]*bzbz[2]; 
  Pzz[3] = (-0.5*bzbz[0]*p_perp[3])+p_perp[3]+0.5*bzbz[0]*p_parallel[3]-0.5*p_perp[0]*bzbz[3]+0.5*p_parallel[0]*bzbz[3]-0.5*bzbz[1]*p_perp[2]+0.5*bzbz[1]*p_parallel[2]-0.5*p_perp[1]*bzbz[2]+0.5*p_parallel[1]*bzbz[2]; 
  Txx[0] = 1.5*Pxx[3]*rho_inv[3]+1.5*Pxx[2]*rho_inv[2]+1.5*Pxx[1]*rho_inv[1]+1.5*Pxx[0]*rho_inv[0]; 
  Txx[1] = 1.5*Pxx[2]*rho_inv[3]+1.5*rho_inv[2]*Pxx[3]+1.5*Pxx[0]*rho_inv[1]+1.5*rho_inv[0]*Pxx[1]; 
  Txx[2] = 1.5*Pxx[1]*rho_inv[3]+1.5*rho_inv[1]*Pxx[3]+1.5*Pxx[0]*rho_inv[2]+1.5*rho_inv[0]*Pxx[2]; 
  Txx[3] = 1.5*Pxx[0]*rho_inv[3]+1.5*rho_inv[0]*Pxx[3]+1.5*Pxx[1]*rho_inv[2]+1.5*rho_inv[1]*Pxx[2]; 

  Tyy[0] = 1.5*Pyy[3]*rho_inv[3]+1.5*Pyy[2]*rho_inv[2]+1.5*Pyy[1]*rho_inv[1]+1.5*Pyy[0]*rho_inv[0]; 
  Tyy[1] = 1.5*Pyy[2]*rho_inv[3]+1.5*rho_inv[2]*Pyy[3]+1.5*Pyy[0]*rho_inv[1]+1.5*rho_inv[0]*Pyy[1]; 
  Tyy[2] = 1.5*Pyy[1]*rho_inv[3]+1.5*rho_inv[1]*Pyy[3]+1.5*Pyy[0]*rho_inv[2]+1.5*rho_inv[0]*Pyy[2]; 
  Tyy[3] = 1.5*Pyy[0]*rho_inv[3]+1.5*rho_inv[0]*Pyy[3]+1.5*Pyy[1]*rho_inv[2]+1.5*rho_inv[1]*Pyy[2]; 

  Tzz[0] = 1.5*Pzz[3]*rho_inv[3]+1.5*Pzz[2]*rho_inv[2]+1.5*Pzz[1]*rho_inv[1]+1.5*Pzz[0]*rho_inv[0]; 
  Tzz[1] = 1.5*Pzz[2]*rho_inv[3]+1.5*rho_inv[2]*Pzz[3]+1.5*Pzz[0]*rho_inv[1]+1.5*rho_inv[0]*Pzz[1]; 
  Tzz[2] = 1.5*Pzz[1]*rho_inv[3]+1.5*rho_inv[1]*Pzz[3]+1.5*Pzz[0]*rho_inv[2]+1.5*rho_inv[0]*Pzz[2]; 
  Tzz[3] = 1.5*Pzz[0]*rho_inv[3]+1.5*rho_inv[0]*Pzz[3]+1.5*Pzz[1]*rho_inv[2]+1.5*rho_inv[1]*Pzz[2]; 
} 
