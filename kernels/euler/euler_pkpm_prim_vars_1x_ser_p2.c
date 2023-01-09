#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
GKYL_CU_DH void euler_pkpm_prim_vars_1x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, 
  double* u_i, double* p_ij, double* T_ij, 
  double* rho_inv, double* T_perp_over_m, double* T_perp_over_m_inv) 
{ 
  // bvar:              Magnetic field unit vector and tensor.
  // vlasov_pkpm_moms:  [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // statevec:          [rho ux, rho uy, rho uz], Fluid input state vector.
  // u_i:               Output flow velocity [ux, uy, uz].
  // p_ij:              Output pressure tensor, p_ij = (p_parallel - p_perp)bb + p_perp I.
  // T_ij:              Output Temperature tensor for penalization T_ij = 3.0*p_ij/rho.
  // rho_inv:           Output 1/rho.
  // T_perp_over_m:     Output p_perp/rho = T_perp/m.
  // T_perp_over_m_inv: Output (T_perp/m)^-1.

  const double *rhoux = &statevec[0]; 
  const double *rhouy = &statevec[3]; 
  const double *rhouz = &statevec[6]; 
  // Parallel pressure is first component of pkpm moment array and unit tensor are last six components of bvar array.
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_parallel = &vlasov_pkpm_moms[3]; 
  const double *p_perp = &vlasov_pkpm_moms[6]; 
  const double *bxbx = &bvar[9]; 
  const double *bxby = &bvar[12]; 
  const double *bxbz = &bvar[15]; 
  const double *byby = &bvar[18]; 
  const double *bybz = &bvar[21]; 
  const double *bzbz = &bvar[24]; 

  double *ux = &u_i[0]; 
  double *uy = &u_i[3]; 
  double *uz = &u_i[6]; 

  double *Pxx = &p_ij[0]; 
  double *Pxy = &p_ij[3]; 
  double *Pxz = &p_ij[6]; 
  double *Pyy = &p_ij[9]; 
  double *Pyz = &p_ij[12]; 
  double *Pzz = &p_ij[15]; 

  double *Txx = &T_ij[0]; 
  double *Tyy = &T_ij[9]; 
  double *Tzz = &T_ij[15]; 

  ser_1x_p2_inv(rho, rho_inv); 
  ux[0] = 0.7071067811865475*rho_inv[2]*rhoux[2]+0.7071067811865475*rho_inv[1]*rhoux[1]+0.7071067811865475*rho_inv[0]*rhoux[0]; 
  ux[1] = 0.6324555320336759*rho_inv[1]*rhoux[2]+0.6324555320336759*rhoux[1]*rho_inv[2]+0.7071067811865475*rho_inv[0]*rhoux[1]+0.7071067811865475*rhoux[0]*rho_inv[1]; 
  ux[2] = 0.4517539514526256*rho_inv[2]*rhoux[2]+0.7071067811865475*rho_inv[0]*rhoux[2]+0.7071067811865475*rhoux[0]*rho_inv[2]+0.6324555320336759*rho_inv[1]*rhoux[1]; 

  uy[0] = 0.7071067811865475*rho_inv[2]*rhouy[2]+0.7071067811865475*rho_inv[1]*rhouy[1]+0.7071067811865475*rho_inv[0]*rhouy[0]; 
  uy[1] = 0.6324555320336759*rho_inv[1]*rhouy[2]+0.6324555320336759*rhouy[1]*rho_inv[2]+0.7071067811865475*rho_inv[0]*rhouy[1]+0.7071067811865475*rhouy[0]*rho_inv[1]; 
  uy[2] = 0.4517539514526256*rho_inv[2]*rhouy[2]+0.7071067811865475*rho_inv[0]*rhouy[2]+0.7071067811865475*rhouy[0]*rho_inv[2]+0.6324555320336759*rho_inv[1]*rhouy[1]; 

  uz[0] = 0.7071067811865475*rho_inv[2]*rhouz[2]+0.7071067811865475*rho_inv[1]*rhouz[1]+0.7071067811865475*rho_inv[0]*rhouz[0]; 
  uz[1] = 0.6324555320336759*rho_inv[1]*rhouz[2]+0.6324555320336759*rhouz[1]*rho_inv[2]+0.7071067811865475*rho_inv[0]*rhouz[1]+0.7071067811865475*rhouz[0]*rho_inv[1]; 
  uz[2] = 0.4517539514526256*rho_inv[2]*rhouz[2]+0.7071067811865475*rho_inv[0]*rhouz[2]+0.7071067811865475*rhouz[0]*rho_inv[2]+0.6324555320336759*rho_inv[1]*rhouz[1]; 

  Pxx[0] = (-0.7071067811865475*bxbx[2]*p_perp[2])+0.7071067811865475*bxbx[2]*p_parallel[2]-0.7071067811865475*bxbx[1]*p_perp[1]+0.7071067811865475*bxbx[1]*p_parallel[1]-0.7071067811865475*bxbx[0]*p_perp[0]+p_perp[0]+0.7071067811865475*bxbx[0]*p_parallel[0]; 
  Pxx[1] = (-0.6324555320336759*bxbx[1]*p_perp[2])+0.6324555320336759*bxbx[1]*p_parallel[2]-0.6324555320336759*p_perp[1]*bxbx[2]+0.6324555320336759*p_parallel[1]*bxbx[2]-0.7071067811865475*bxbx[0]*p_perp[1]+p_perp[1]+0.7071067811865475*bxbx[0]*p_parallel[1]-0.7071067811865475*p_perp[0]*bxbx[1]+0.7071067811865475*p_parallel[0]*bxbx[1]; 
  Pxx[2] = (-0.4517539514526256*bxbx[2]*p_perp[2])-0.7071067811865475*bxbx[0]*p_perp[2]+p_perp[2]+0.4517539514526256*bxbx[2]*p_parallel[2]+0.7071067811865475*bxbx[0]*p_parallel[2]-0.7071067811865475*p_perp[0]*bxbx[2]+0.7071067811865475*p_parallel[0]*bxbx[2]-0.6324555320336759*bxbx[1]*p_perp[1]+0.6324555320336759*bxbx[1]*p_parallel[1]; 

  Pxy[0] = (-0.7071067811865475*bxby[2]*p_perp[2])+0.7071067811865475*bxby[2]*p_parallel[2]-0.7071067811865475*bxby[1]*p_perp[1]+0.7071067811865475*bxby[1]*p_parallel[1]-0.7071067811865475*bxby[0]*p_perp[0]+0.7071067811865475*bxby[0]*p_parallel[0]; 
  Pxy[1] = (-0.6324555320336759*bxby[1]*p_perp[2])+0.6324555320336759*bxby[1]*p_parallel[2]-0.6324555320336759*p_perp[1]*bxby[2]+0.6324555320336759*p_parallel[1]*bxby[2]-0.7071067811865475*bxby[0]*p_perp[1]+0.7071067811865475*bxby[0]*p_parallel[1]-0.7071067811865475*p_perp[0]*bxby[1]+0.7071067811865475*p_parallel[0]*bxby[1]; 
  Pxy[2] = (-0.4517539514526256*bxby[2]*p_perp[2])-0.7071067811865475*bxby[0]*p_perp[2]+0.4517539514526256*bxby[2]*p_parallel[2]+0.7071067811865475*bxby[0]*p_parallel[2]-0.7071067811865475*p_perp[0]*bxby[2]+0.7071067811865475*p_parallel[0]*bxby[2]-0.6324555320336759*bxby[1]*p_perp[1]+0.6324555320336759*bxby[1]*p_parallel[1]; 

  Pxz[0] = (-0.7071067811865475*bxbz[2]*p_perp[2])+0.7071067811865475*bxbz[2]*p_parallel[2]-0.7071067811865475*bxbz[1]*p_perp[1]+0.7071067811865475*bxbz[1]*p_parallel[1]-0.7071067811865475*bxbz[0]*p_perp[0]+0.7071067811865475*bxbz[0]*p_parallel[0]; 
  Pxz[1] = (-0.6324555320336759*bxbz[1]*p_perp[2])+0.6324555320336759*bxbz[1]*p_parallel[2]-0.6324555320336759*p_perp[1]*bxbz[2]+0.6324555320336759*p_parallel[1]*bxbz[2]-0.7071067811865475*bxbz[0]*p_perp[1]+0.7071067811865475*bxbz[0]*p_parallel[1]-0.7071067811865475*p_perp[0]*bxbz[1]+0.7071067811865475*p_parallel[0]*bxbz[1]; 
  Pxz[2] = (-0.4517539514526256*bxbz[2]*p_perp[2])-0.7071067811865475*bxbz[0]*p_perp[2]+0.4517539514526256*bxbz[2]*p_parallel[2]+0.7071067811865475*bxbz[0]*p_parallel[2]-0.7071067811865475*p_perp[0]*bxbz[2]+0.7071067811865475*p_parallel[0]*bxbz[2]-0.6324555320336759*bxbz[1]*p_perp[1]+0.6324555320336759*bxbz[1]*p_parallel[1]; 

  Pyy[0] = (-0.7071067811865475*byby[2]*p_perp[2])+0.7071067811865475*byby[2]*p_parallel[2]-0.7071067811865475*byby[1]*p_perp[1]+0.7071067811865475*byby[1]*p_parallel[1]-0.7071067811865475*byby[0]*p_perp[0]+p_perp[0]+0.7071067811865475*byby[0]*p_parallel[0]; 
  Pyy[1] = (-0.6324555320336759*byby[1]*p_perp[2])+0.6324555320336759*byby[1]*p_parallel[2]-0.6324555320336759*p_perp[1]*byby[2]+0.6324555320336759*p_parallel[1]*byby[2]-0.7071067811865475*byby[0]*p_perp[1]+p_perp[1]+0.7071067811865475*byby[0]*p_parallel[1]-0.7071067811865475*p_perp[0]*byby[1]+0.7071067811865475*p_parallel[0]*byby[1]; 
  Pyy[2] = (-0.4517539514526256*byby[2]*p_perp[2])-0.7071067811865475*byby[0]*p_perp[2]+p_perp[2]+0.4517539514526256*byby[2]*p_parallel[2]+0.7071067811865475*byby[0]*p_parallel[2]-0.7071067811865475*p_perp[0]*byby[2]+0.7071067811865475*p_parallel[0]*byby[2]-0.6324555320336759*byby[1]*p_perp[1]+0.6324555320336759*byby[1]*p_parallel[1]; 

  Pyz[0] = (-0.7071067811865475*bybz[2]*p_perp[2])+0.7071067811865475*bybz[2]*p_parallel[2]-0.7071067811865475*bybz[1]*p_perp[1]+0.7071067811865475*bybz[1]*p_parallel[1]-0.7071067811865475*bybz[0]*p_perp[0]+0.7071067811865475*bybz[0]*p_parallel[0]; 
  Pyz[1] = (-0.6324555320336759*bybz[1]*p_perp[2])+0.6324555320336759*bybz[1]*p_parallel[2]-0.6324555320336759*p_perp[1]*bybz[2]+0.6324555320336759*p_parallel[1]*bybz[2]-0.7071067811865475*bybz[0]*p_perp[1]+0.7071067811865475*bybz[0]*p_parallel[1]-0.7071067811865475*p_perp[0]*bybz[1]+0.7071067811865475*p_parallel[0]*bybz[1]; 
  Pyz[2] = (-0.4517539514526256*bybz[2]*p_perp[2])-0.7071067811865475*bybz[0]*p_perp[2]+0.4517539514526256*bybz[2]*p_parallel[2]+0.7071067811865475*bybz[0]*p_parallel[2]-0.7071067811865475*p_perp[0]*bybz[2]+0.7071067811865475*p_parallel[0]*bybz[2]-0.6324555320336759*bybz[1]*p_perp[1]+0.6324555320336759*bybz[1]*p_parallel[1]; 

  Pzz[0] = (-0.7071067811865475*bzbz[2]*p_perp[2])+0.7071067811865475*bzbz[2]*p_parallel[2]-0.7071067811865475*bzbz[1]*p_perp[1]+0.7071067811865475*bzbz[1]*p_parallel[1]-0.7071067811865475*bzbz[0]*p_perp[0]+p_perp[0]+0.7071067811865475*bzbz[0]*p_parallel[0]; 
  Pzz[1] = (-0.6324555320336759*bzbz[1]*p_perp[2])+0.6324555320336759*bzbz[1]*p_parallel[2]-0.6324555320336759*p_perp[1]*bzbz[2]+0.6324555320336759*p_parallel[1]*bzbz[2]-0.7071067811865475*bzbz[0]*p_perp[1]+p_perp[1]+0.7071067811865475*bzbz[0]*p_parallel[1]-0.7071067811865475*p_perp[0]*bzbz[1]+0.7071067811865475*p_parallel[0]*bzbz[1]; 
  Pzz[2] = (-0.4517539514526256*bzbz[2]*p_perp[2])-0.7071067811865475*bzbz[0]*p_perp[2]+p_perp[2]+0.4517539514526256*bzbz[2]*p_parallel[2]+0.7071067811865475*bzbz[0]*p_parallel[2]-0.7071067811865475*p_perp[0]*bzbz[2]+0.7071067811865475*p_parallel[0]*bzbz[2]-0.6324555320336759*bzbz[1]*p_perp[1]+0.6324555320336759*bzbz[1]*p_parallel[1]; 
  Txx[0] = 2.121320343559642*Pxx[2]*rho_inv[2]+2.121320343559642*Pxx[1]*rho_inv[1]+2.121320343559642*Pxx[0]*rho_inv[0]; 
  Txx[1] = 1.897366596101028*Pxx[1]*rho_inv[2]+1.897366596101028*rho_inv[1]*Pxx[2]+2.121320343559642*Pxx[0]*rho_inv[1]+2.121320343559642*rho_inv[0]*Pxx[1]; 
  Txx[2] = 1.355261854357877*Pxx[2]*rho_inv[2]+2.121320343559642*Pxx[0]*rho_inv[2]+2.121320343559642*rho_inv[0]*Pxx[2]+1.897366596101028*Pxx[1]*rho_inv[1]; 

  Tyy[0] = 2.121320343559642*Pyy[2]*rho_inv[2]+2.121320343559642*Pyy[1]*rho_inv[1]+2.121320343559642*Pyy[0]*rho_inv[0]; 
  Tyy[1] = 1.897366596101028*Pyy[1]*rho_inv[2]+1.897366596101028*rho_inv[1]*Pyy[2]+2.121320343559642*Pyy[0]*rho_inv[1]+2.121320343559642*rho_inv[0]*Pyy[1]; 
  Tyy[2] = 1.355261854357877*Pyy[2]*rho_inv[2]+2.121320343559642*Pyy[0]*rho_inv[2]+2.121320343559642*rho_inv[0]*Pyy[2]+1.897366596101028*Pyy[1]*rho_inv[1]; 

  Tzz[0] = 2.121320343559642*Pzz[2]*rho_inv[2]+2.121320343559642*Pzz[1]*rho_inv[1]+2.121320343559642*Pzz[0]*rho_inv[0]; 
  Tzz[1] = 1.897366596101028*Pzz[1]*rho_inv[2]+1.897366596101028*rho_inv[1]*Pzz[2]+2.121320343559642*Pzz[0]*rho_inv[1]+2.121320343559642*rho_inv[0]*Pzz[1]; 
  Tzz[2] = 1.355261854357877*Pzz[2]*rho_inv[2]+2.121320343559642*Pzz[0]*rho_inv[2]+2.121320343559642*rho_inv[0]*Pzz[2]+1.897366596101028*Pzz[1]*rho_inv[1]; 

  T_perp_over_m[0] = 0.7071067811865475*p_perp[2]*rho_inv[2]+0.7071067811865475*p_perp[1]*rho_inv[1]+0.7071067811865475*p_perp[0]*rho_inv[0]; 
  T_perp_over_m[1] = 0.6324555320336759*p_perp[1]*rho_inv[2]+0.6324555320336759*rho_inv[1]*p_perp[2]+0.7071067811865475*p_perp[0]*rho_inv[1]+0.7071067811865475*rho_inv[0]*p_perp[1]; 
  T_perp_over_m[2] = 0.4517539514526256*p_perp[2]*rho_inv[2]+0.7071067811865475*p_perp[0]*rho_inv[2]+0.7071067811865475*rho_inv[0]*p_perp[2]+0.6324555320336759*p_perp[1]*rho_inv[1]; 
  ser_1x_p2_inv(T_perp_over_m, T_perp_over_m_inv); 
} 
