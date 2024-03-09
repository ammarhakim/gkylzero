#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_source_1x_ser_p1(const double* qmem, const double* fluid, const double* p_ij, double* GKYL_RESTRICT out) 
{ 
  // qmem:  q/m*EM fields.
  // fluid: [rho, rho ux, rho uy, rho uz...], Fluid input state vector.
  // p_ij:  Input pressure tensor (only used by 10 moment).
  // out:   Output increment

  const double *Ex = &qmem[0]; 
  const double *Ey = &qmem[2]; 
  const double *Ez = &qmem[4]; 
  const double *Bx = &qmem[6]; 
  const double *By = &qmem[8]; 
  const double *Bz = &qmem[10]; 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[2]; 
  const double *rhouy = &fluid[4]; 
  const double *rhouz = &fluid[6]; 

  double *out_rhoux = &out[2]; 
  double *out_rhouy = &out[4]; 
  double *out_rhouz = &out[6]; 
  double *out_energy = &out[8]; 

  out_rhoux[0] += (-0.7071067811865475*By[1]*rhouz[1])+0.7071067811865475*Bz[1]*rhouy[1]+0.7071067811865475*Ex[1]*rho[1]-0.7071067811865475*By[0]*rhouz[0]+0.7071067811865475*Bz[0]*rhouy[0]+0.7071067811865475*Ex[0]*rho[0]; 
  out_rhoux[1] += (-0.7071067811865475*By[0]*rhouz[1])+0.7071067811865475*Bz[0]*rhouy[1]+0.7071067811865475*Ex[0]*rho[1]+0.7071067811865475*rho[0]*Ex[1]+0.7071067811865475*rhouy[0]*Bz[1]-0.7071067811865475*rhouz[0]*By[1]; 

  out_rhouy[0] += 0.7071067811865475*Bx[1]*rhouz[1]-0.7071067811865475*Bz[1]*rhoux[1]+0.7071067811865475*Ey[1]*rho[1]+0.7071067811865475*Bx[0]*rhouz[0]-0.7071067811865475*Bz[0]*rhoux[0]+0.7071067811865475*Ey[0]*rho[0]; 
  out_rhouy[1] += 0.7071067811865475*Bx[0]*rhouz[1]-0.7071067811865475*Bz[0]*rhoux[1]+0.7071067811865475*Ey[0]*rho[1]+0.7071067811865475*rho[0]*Ey[1]-0.7071067811865475*rhoux[0]*Bz[1]+0.7071067811865475*rhouz[0]*Bx[1]; 

  out_rhouz[0] += (-0.7071067811865475*Bx[1]*rhouy[1])+0.7071067811865475*By[1]*rhoux[1]+0.7071067811865475*Ez[1]*rho[1]-0.7071067811865475*Bx[0]*rhouy[0]+0.7071067811865475*By[0]*rhoux[0]+0.7071067811865475*Ez[0]*rho[0]; 
  out_rhouz[1] += (-0.7071067811865475*Bx[0]*rhouy[1])+0.7071067811865475*By[0]*rhoux[1]+0.7071067811865475*Ez[0]*rho[1]+0.7071067811865475*rho[0]*Ez[1]+0.7071067811865475*rhoux[0]*By[1]-0.7071067811865475*rhouy[0]*Bx[1]; 

  out_energy[0] += 0.7071067811865475*Ez[1]*rhouz[1]+0.7071067811865475*Ey[1]*rhouy[1]+0.7071067811865475*Ex[1]*rhoux[1]+0.7071067811865475*Ez[0]*rhouz[0]+0.7071067811865475*Ey[0]*rhouy[0]+0.7071067811865475*Ex[0]*rhoux[0]; 
  out_energy[1] += 0.7071067811865475*Ez[0]*rhouz[1]+0.7071067811865475*Ey[0]*rhouy[1]+0.7071067811865475*Ex[0]*rhoux[1]+0.7071067811865475*rhouz[0]*Ez[1]+0.7071067811865475*rhouy[0]*Ey[1]+0.7071067811865475*rhoux[0]*Ex[1]; 

} 
