#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p2(const double *qmem, const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *rhou_perp_i, const double *p_perp, double* out) 
{ 
  // qmem:             q/m*EM fields.
  // nu:               Collisionality.
  // nu_vth_sq:        nu*vth^2 = nu*T/m.
  // vlasov_pkpm_moms: [rho, p_parallel, q_parallel], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       Input fluid variables.
  // rhou_perp_i:      Input perpendicular momentum density, rhou - (rhou . b)b = [rhou_perp_x, rhou_perp_y, rhou_perp_z].
  // p_perp:           Input perpendicular pressure (E_perp = 1/2 rho u_perp^2 + p_perp).
  // out:              Output increment
  const double *rho = &vlasov_pkpm_moms[0]; 

  const double *Ex = &qmem[0]; 
  const double *Ey = &qmem[3]; 
  const double *Ez = &qmem[6]; 
  const double *Bx = &qmem[9]; 
  const double *By = &qmem[12]; 
  const double *Bz = &qmem[15]; 

  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[3]; 
  const double *rhouz = &euler_pkpm[6]; 

  const double *rhou_perp_x = &rhou_perp_i[0]; 
  const double *rhou_perp_y = &rhou_perp_i[3]; 
  const double *rhou_perp_z = &rhou_perp_i[6]; 

  double *outrhoux = &out[0]; 
  double *outrhouy = &out[3]; 
  double *outrhouz = &out[6]; 
  double *outE_perp = &out[9]; 

  outrhoux[0] += (-0.7071067811865475*By[2]*rhouz[2])+0.7071067811865475*Bz[2]*rhouy[2]+0.7071067811865475*Ex[2]*rho[2]-0.7071067811865475*By[1]*rhouz[1]+0.7071067811865475*Bz[1]*rhouy[1]+0.7071067811865475*Ex[1]*rho[1]-0.7071067811865475*By[0]*rhouz[0]+0.7071067811865475*Bz[0]*rhouy[0]+0.7071067811865475*Ex[0]*rho[0]; 
  outrhoux[1] += (-0.6324555320336759*By[1]*rhouz[2])+0.6324555320336759*Bz[1]*rhouy[2]+0.6324555320336759*Ex[1]*rho[2]+0.6324555320336759*rho[1]*Ex[2]+0.6324555320336759*rhouy[1]*Bz[2]-0.6324555320336759*rhouz[1]*By[2]-0.7071067811865475*By[0]*rhouz[1]+0.7071067811865475*Bz[0]*rhouy[1]+0.7071067811865475*Ex[0]*rho[1]+0.7071067811865475*rho[0]*Ex[1]+0.7071067811865475*rhouy[0]*Bz[1]-0.7071067811865475*rhouz[0]*By[1]; 
  outrhoux[2] += (-0.4517539514526256*By[2]*rhouz[2])-0.7071067811865475*By[0]*rhouz[2]+0.4517539514526256*Bz[2]*rhouy[2]+0.7071067811865475*Bz[0]*rhouy[2]+0.4517539514526256*Ex[2]*rho[2]+0.7071067811865475*Ex[0]*rho[2]+0.7071067811865475*rho[0]*Ex[2]+0.7071067811865475*rhouy[0]*Bz[2]-0.7071067811865475*rhouz[0]*By[2]-0.6324555320336759*By[1]*rhouz[1]+0.6324555320336759*Bz[1]*rhouy[1]+0.6324555320336759*Ex[1]*rho[1]; 

  outrhouy[0] += 0.7071067811865475*Bx[2]*rhouz[2]-0.7071067811865475*Bz[2]*rhoux[2]+0.7071067811865475*Ey[2]*rho[2]+0.7071067811865475*Bx[1]*rhouz[1]-0.7071067811865475*Bz[1]*rhoux[1]+0.7071067811865475*Ey[1]*rho[1]+0.7071067811865475*Bx[0]*rhouz[0]-0.7071067811865475*Bz[0]*rhoux[0]+0.7071067811865475*Ey[0]*rho[0]; 
  outrhouy[1] += 0.6324555320336759*Bx[1]*rhouz[2]-0.6324555320336759*Bz[1]*rhoux[2]+0.6324555320336759*Ey[1]*rho[2]+0.6324555320336759*rho[1]*Ey[2]-0.6324555320336759*rhoux[1]*Bz[2]+0.6324555320336759*rhouz[1]*Bx[2]+0.7071067811865475*Bx[0]*rhouz[1]-0.7071067811865475*Bz[0]*rhoux[1]+0.7071067811865475*Ey[0]*rho[1]+0.7071067811865475*rho[0]*Ey[1]-0.7071067811865475*rhoux[0]*Bz[1]+0.7071067811865475*rhouz[0]*Bx[1]; 
  outrhouy[2] += 0.4517539514526256*Bx[2]*rhouz[2]+0.7071067811865475*Bx[0]*rhouz[2]-0.4517539514526256*Bz[2]*rhoux[2]-0.7071067811865475*Bz[0]*rhoux[2]+0.4517539514526256*Ey[2]*rho[2]+0.7071067811865475*Ey[0]*rho[2]+0.7071067811865475*rho[0]*Ey[2]-0.7071067811865475*rhoux[0]*Bz[2]+0.7071067811865475*rhouz[0]*Bx[2]+0.6324555320336759*Bx[1]*rhouz[1]-0.6324555320336759*Bz[1]*rhoux[1]+0.6324555320336759*Ey[1]*rho[1]; 

  outrhouz[0] += (-0.7071067811865475*Bx[2]*rhouy[2])+0.7071067811865475*By[2]*rhoux[2]+0.7071067811865475*Ez[2]*rho[2]-0.7071067811865475*Bx[1]*rhouy[1]+0.7071067811865475*By[1]*rhoux[1]+0.7071067811865475*Ez[1]*rho[1]-0.7071067811865475*Bx[0]*rhouy[0]+0.7071067811865475*By[0]*rhoux[0]+0.7071067811865475*Ez[0]*rho[0]; 
  outrhouz[1] += (-0.6324555320336759*Bx[1]*rhouy[2])+0.6324555320336759*By[1]*rhoux[2]+0.6324555320336759*Ez[1]*rho[2]+0.6324555320336759*rho[1]*Ez[2]+0.6324555320336759*rhoux[1]*By[2]-0.6324555320336759*rhouy[1]*Bx[2]-0.7071067811865475*Bx[0]*rhouy[1]+0.7071067811865475*By[0]*rhoux[1]+0.7071067811865475*Ez[0]*rho[1]+0.7071067811865475*rho[0]*Ez[1]+0.7071067811865475*rhoux[0]*By[1]-0.7071067811865475*rhouy[0]*Bx[1]; 
  outrhouz[2] += (-0.4517539514526256*Bx[2]*rhouy[2])-0.7071067811865475*Bx[0]*rhouy[2]+0.4517539514526256*By[2]*rhoux[2]+0.7071067811865475*By[0]*rhoux[2]+0.4517539514526256*Ez[2]*rho[2]+0.7071067811865475*Ez[0]*rho[2]+0.7071067811865475*rho[0]*Ez[2]+0.7071067811865475*rhoux[0]*By[2]-0.7071067811865475*rhouy[0]*Bx[2]-0.6324555320336759*Bx[1]*rhouy[1]+0.6324555320336759*By[1]*rhoux[1]+0.6324555320336759*Ez[1]*rho[1]; 

  outE_perp[0] += 0.7071067811865475*Ez[2]*rhou_perp_z[2]+0.7071067811865475*Ey[2]*rhou_perp_y[2]+0.7071067811865475*Ex[2]*rhou_perp_x[2]+0.7071067811865475*nu_vth_sq[2]*rho[2]-0.7071067811865475*nu[2]*p_perp[2]+0.7071067811865475*Ez[1]*rhou_perp_z[1]+0.7071067811865475*Ey[1]*rhou_perp_y[1]+0.7071067811865475*Ex[1]*rhou_perp_x[1]+0.7071067811865475*nu_vth_sq[1]*rho[1]-0.7071067811865475*nu[1]*p_perp[1]+0.7071067811865475*Ez[0]*rhou_perp_z[0]+0.7071067811865475*Ey[0]*rhou_perp_y[0]+0.7071067811865475*Ex[0]*rhou_perp_x[0]+0.7071067811865475*nu_vth_sq[0]*rho[0]-0.7071067811865475*nu[0]*p_perp[0]; 
  outE_perp[1] += 0.6324555320336759*Ez[1]*rhou_perp_z[2]+0.6324555320336759*Ey[1]*rhou_perp_y[2]+0.6324555320336759*Ex[1]*rhou_perp_x[2]+0.6324555320336759*nu_vth_sq[1]*rho[2]-0.6324555320336759*nu[1]*p_perp[2]+0.6324555320336759*rho[1]*nu_vth_sq[2]-0.6324555320336759*p_perp[1]*nu[2]+0.6324555320336759*rhou_perp_z[1]*Ez[2]+0.6324555320336759*rhou_perp_y[1]*Ey[2]+0.6324555320336759*rhou_perp_x[1]*Ex[2]+0.7071067811865475*Ez[0]*rhou_perp_z[1]+0.7071067811865475*Ey[0]*rhou_perp_y[1]+0.7071067811865475*Ex[0]*rhou_perp_x[1]+0.7071067811865475*nu_vth_sq[0]*rho[1]-0.7071067811865475*nu[0]*p_perp[1]+0.7071067811865475*rho[0]*nu_vth_sq[1]-0.7071067811865475*p_perp[0]*nu[1]+0.7071067811865475*rhou_perp_z[0]*Ez[1]+0.7071067811865475*rhou_perp_y[0]*Ey[1]+0.7071067811865475*rhou_perp_x[0]*Ex[1]; 
  outE_perp[2] += 0.4517539514526256*Ez[2]*rhou_perp_z[2]+0.7071067811865475*Ez[0]*rhou_perp_z[2]+0.4517539514526256*Ey[2]*rhou_perp_y[2]+0.7071067811865475*Ey[0]*rhou_perp_y[2]+0.4517539514526256*Ex[2]*rhou_perp_x[2]+0.7071067811865475*Ex[0]*rhou_perp_x[2]+0.4517539514526256*nu_vth_sq[2]*rho[2]+0.7071067811865475*nu_vth_sq[0]*rho[2]-0.4517539514526256*nu[2]*p_perp[2]-0.7071067811865475*nu[0]*p_perp[2]+0.7071067811865475*rho[0]*nu_vth_sq[2]-0.7071067811865475*p_perp[0]*nu[2]+0.7071067811865475*rhou_perp_z[0]*Ez[2]+0.7071067811865475*rhou_perp_y[0]*Ey[2]+0.7071067811865475*rhou_perp_x[0]*Ex[2]+0.6324555320336759*Ez[1]*rhou_perp_z[1]+0.6324555320336759*Ey[1]*rhou_perp_y[1]+0.6324555320336759*Ex[1]*rhou_perp_x[1]+0.6324555320336759*nu_vth_sq[1]*rho[1]-0.6324555320336759*nu[1]*p_perp[1]; 

} 
