#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_source_2x_ser_p1(const double* qmem, const double* fluid, const double* p_ij, double* GKYL_RESTRICT out) 
{ 
  // qmem:  q/m*EM fields.
  // fluid: [rho, rho ux, rho uy, rho uz...], Fluid input state vector.
  // p_ij:  Input pressure tensor (only used by 10 moment).
  // out:   Output increment

  const double *Ex = &qmem[0]; 
  const double *Ey = &qmem[4]; 
  const double *Ez = &qmem[8]; 
  const double *Bx = &qmem[12]; 
  const double *By = &qmem[16]; 
  const double *Bz = &qmem[20]; 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[4]; 
  const double *rhouy = &fluid[8]; 
  const double *rhouz = &fluid[12]; 

  double *out_rhoux = &out[4]; 
  double *out_rhouy = &out[8]; 
  double *out_rhouz = &out[12]; 
  double *out_energy = &out[16]; 

  out_rhoux[0] += (-0.5*By[3]*rhouz[3])+0.5*Bz[3]*rhouy[3]+0.5*Ex[3]*rho[3]-0.5*By[2]*rhouz[2]+0.5*Bz[2]*rhouy[2]+0.5*Ex[2]*rho[2]-0.5*By[1]*rhouz[1]+0.5*Bz[1]*rhouy[1]+0.5*Ex[1]*rho[1]-0.5*By[0]*rhouz[0]+0.5*Bz[0]*rhouy[0]+0.5*Ex[0]*rho[0]; 
  out_rhoux[1] += (-0.5*By[2]*rhouz[3])+0.5*Bz[2]*rhouy[3]+0.5*Ex[2]*rho[3]+0.5*rho[2]*Ex[3]+0.5*rhouy[2]*Bz[3]-0.5*rhouz[2]*By[3]-0.5*By[0]*rhouz[1]+0.5*Bz[0]*rhouy[1]+0.5*Ex[0]*rho[1]+0.5*rho[0]*Ex[1]+0.5*rhouy[0]*Bz[1]-0.5*rhouz[0]*By[1]; 
  out_rhoux[2] += (-0.5*By[1]*rhouz[3])+0.5*Bz[1]*rhouy[3]+0.5*Ex[1]*rho[3]+0.5*rho[1]*Ex[3]+0.5*rhouy[1]*Bz[3]-0.5*rhouz[1]*By[3]-0.5*By[0]*rhouz[2]+0.5*Bz[0]*rhouy[2]+0.5*Ex[0]*rho[2]+0.5*rho[0]*Ex[2]+0.5*rhouy[0]*Bz[2]-0.5*rhouz[0]*By[2]; 
  out_rhoux[3] += (-0.5*By[0]*rhouz[3])+0.5*Bz[0]*rhouy[3]+0.5*Ex[0]*rho[3]+0.5*rho[0]*Ex[3]+0.5*rhouy[0]*Bz[3]-0.5*rhouz[0]*By[3]-0.5*By[1]*rhouz[2]+0.5*Bz[1]*rhouy[2]+0.5*Ex[1]*rho[2]+0.5*rho[1]*Ex[2]+0.5*rhouy[1]*Bz[2]-0.5*rhouz[1]*By[2]; 

  out_rhouy[0] += 0.5*Bx[3]*rhouz[3]-0.5*Bz[3]*rhoux[3]+0.5*Ey[3]*rho[3]+0.5*Bx[2]*rhouz[2]-0.5*Bz[2]*rhoux[2]+0.5*Ey[2]*rho[2]+0.5*Bx[1]*rhouz[1]-0.5*Bz[1]*rhoux[1]+0.5*Ey[1]*rho[1]+0.5*Bx[0]*rhouz[0]-0.5*Bz[0]*rhoux[0]+0.5*Ey[0]*rho[0]; 
  out_rhouy[1] += 0.5*Bx[2]*rhouz[3]-0.5*Bz[2]*rhoux[3]+0.5*Ey[2]*rho[3]+0.5*rho[2]*Ey[3]-0.5*rhoux[2]*Bz[3]+0.5*rhouz[2]*Bx[3]+0.5*Bx[0]*rhouz[1]-0.5*Bz[0]*rhoux[1]+0.5*Ey[0]*rho[1]+0.5*rho[0]*Ey[1]-0.5*rhoux[0]*Bz[1]+0.5*rhouz[0]*Bx[1]; 
  out_rhouy[2] += 0.5*Bx[1]*rhouz[3]-0.5*Bz[1]*rhoux[3]+0.5*Ey[1]*rho[3]+0.5*rho[1]*Ey[3]-0.5*rhoux[1]*Bz[3]+0.5*rhouz[1]*Bx[3]+0.5*Bx[0]*rhouz[2]-0.5*Bz[0]*rhoux[2]+0.5*Ey[0]*rho[2]+0.5*rho[0]*Ey[2]-0.5*rhoux[0]*Bz[2]+0.5*rhouz[0]*Bx[2]; 
  out_rhouy[3] += 0.5*Bx[0]*rhouz[3]-0.5*Bz[0]*rhoux[3]+0.5*Ey[0]*rho[3]+0.5*rho[0]*Ey[3]-0.5*rhoux[0]*Bz[3]+0.5*rhouz[0]*Bx[3]+0.5*Bx[1]*rhouz[2]-0.5*Bz[1]*rhoux[2]+0.5*Ey[1]*rho[2]+0.5*rho[1]*Ey[2]-0.5*rhoux[1]*Bz[2]+0.5*rhouz[1]*Bx[2]; 

  out_rhouz[0] += (-0.5*Bx[3]*rhouy[3])+0.5*By[3]*rhoux[3]+0.5*Ez[3]*rho[3]-0.5*Bx[2]*rhouy[2]+0.5*By[2]*rhoux[2]+0.5*Ez[2]*rho[2]-0.5*Bx[1]*rhouy[1]+0.5*By[1]*rhoux[1]+0.5*Ez[1]*rho[1]-0.5*Bx[0]*rhouy[0]+0.5*By[0]*rhoux[0]+0.5*Ez[0]*rho[0]; 
  out_rhouz[1] += (-0.5*Bx[2]*rhouy[3])+0.5*By[2]*rhoux[3]+0.5*Ez[2]*rho[3]+0.5*rho[2]*Ez[3]+0.5*rhoux[2]*By[3]-0.5*rhouy[2]*Bx[3]-0.5*Bx[0]*rhouy[1]+0.5*By[0]*rhoux[1]+0.5*Ez[0]*rho[1]+0.5*rho[0]*Ez[1]+0.5*rhoux[0]*By[1]-0.5*rhouy[0]*Bx[1]; 
  out_rhouz[2] += (-0.5*Bx[1]*rhouy[3])+0.5*By[1]*rhoux[3]+0.5*Ez[1]*rho[3]+0.5*rho[1]*Ez[3]+0.5*rhoux[1]*By[3]-0.5*rhouy[1]*Bx[3]-0.5*Bx[0]*rhouy[2]+0.5*By[0]*rhoux[2]+0.5*Ez[0]*rho[2]+0.5*rho[0]*Ez[2]+0.5*rhoux[0]*By[2]-0.5*rhouy[0]*Bx[2]; 
  out_rhouz[3] += (-0.5*Bx[0]*rhouy[3])+0.5*By[0]*rhoux[3]+0.5*Ez[0]*rho[3]+0.5*rho[0]*Ez[3]+0.5*rhoux[0]*By[3]-0.5*rhouy[0]*Bx[3]-0.5*Bx[1]*rhouy[2]+0.5*By[1]*rhoux[2]+0.5*Ez[1]*rho[2]+0.5*rho[1]*Ez[2]+0.5*rhoux[1]*By[2]-0.5*rhouy[1]*Bx[2]; 

  out_energy[0] += 0.5*Ez[3]*rhouz[3]+0.5*Ey[3]*rhouy[3]+0.5*Ex[3]*rhoux[3]+0.5*Ez[2]*rhouz[2]+0.5*Ey[2]*rhouy[2]+0.5*Ex[2]*rhoux[2]+0.5*Ez[1]*rhouz[1]+0.5*Ey[1]*rhouy[1]+0.5*Ex[1]*rhoux[1]+0.5*Ez[0]*rhouz[0]+0.5*Ey[0]*rhouy[0]+0.5*Ex[0]*rhoux[0]; 
  out_energy[1] += 0.5*Ez[2]*rhouz[3]+0.5*Ey[2]*rhouy[3]+0.5*Ex[2]*rhoux[3]+0.5*rhouz[2]*Ez[3]+0.5*rhouy[2]*Ey[3]+0.5*rhoux[2]*Ex[3]+0.5*Ez[0]*rhouz[1]+0.5*Ey[0]*rhouy[1]+0.5*Ex[0]*rhoux[1]+0.5*rhouz[0]*Ez[1]+0.5*rhouy[0]*Ey[1]+0.5*rhoux[0]*Ex[1]; 
  out_energy[2] += 0.5*Ez[1]*rhouz[3]+0.5*Ey[1]*rhouy[3]+0.5*Ex[1]*rhoux[3]+0.5*rhouz[1]*Ez[3]+0.5*rhouy[1]*Ey[3]+0.5*rhoux[1]*Ex[3]+0.5*Ez[0]*rhouz[2]+0.5*Ey[0]*rhouy[2]+0.5*Ex[0]*rhoux[2]+0.5*rhouz[0]*Ez[2]+0.5*rhouy[0]*Ey[2]+0.5*rhoux[0]*Ex[2]; 
  out_energy[3] += 0.5*Ez[0]*rhouz[3]+0.5*Ey[0]*rhouy[3]+0.5*Ex[0]*rhoux[3]+0.5*rhouz[0]*Ez[3]+0.5*rhouy[0]*Ey[3]+0.5*rhoux[0]*Ex[3]+0.5*Ez[1]*rhouz[2]+0.5*Ey[1]*rhouy[2]+0.5*Ex[1]*rhoux[2]+0.5*rhouz[1]*Ez[2]+0.5*rhouy[1]*Ey[2]+0.5*rhoux[1]*Ex[2]; 

} 
