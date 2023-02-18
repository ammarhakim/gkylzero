#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pressure_1x_ser_p2(const double gas_gamma, const double *uvar, const double *statevec, double* pressure) 
{ 
  // gas_gamma: Adiabatic index.
  // uvar: [ux, uy, uz], Fluid flow.
  // statevec: [rho, rho ux, rho uy, rho uz, energy], Fluid input state vector.
  // pressure: Output pressure variable, p = (gas_gamma - 1)*(Energy - 0.5*rho*u^2).

  const double *rho = &statevec[0]; 
  const double *rhou0 = &statevec[3]; 
  const double *rhou1 = &statevec[6]; 
  const double *rhou2 = &statevec[9]; 
  const double *energy = &statevec[12]; 
  const double *uvar0 = &uvar[0]; 
  const double *uvar1 = &uvar[3]; 
  const double *uvar2 = &uvar[6]; 

  double pressure_fac = (gas_gamma-1.0); 

  pressure[0] = ((-0.3535533905932737*rhou2[2]*uvar2[2])-0.3535533905932737*rhou1[2]*uvar1[2]-0.3535533905932737*rhou0[2]*uvar0[2]-0.3535533905932737*rhou2[1]*uvar2[1]-0.3535533905932737*rhou1[1]*uvar1[1]-0.3535533905932737*rhou0[1]*uvar0[1]-0.3535533905932737*rhou2[0]*uvar2[0]-0.3535533905932737*rhou1[0]*uvar1[0]-0.3535533905932737*rhou0[0]*uvar0[0]+energy[0])*pressure_fac; 
  pressure[1] = ((-0.3162277660168379*rhou2[1]*uvar2[2])-0.3162277660168379*rhou1[1]*uvar1[2]-0.3162277660168379*rhou0[1]*uvar0[2]-0.3162277660168379*uvar2[1]*rhou2[2]-0.3162277660168379*uvar1[1]*rhou1[2]-0.3162277660168379*uvar0[1]*rhou0[2]-0.3535533905932737*rhou2[0]*uvar2[1]-0.3535533905932737*rhou1[0]*uvar1[1]-0.3535533905932737*rhou0[0]*uvar0[1]-0.3535533905932737*uvar2[0]*rhou2[1]-0.3535533905932737*uvar1[0]*rhou1[1]-0.3535533905932737*uvar0[0]*rhou0[1]+energy[1])*pressure_fac; 
  pressure[2] = ((-0.2258769757263128*rhou2[2]*uvar2[2])-0.3535533905932737*rhou2[0]*uvar2[2]-0.2258769757263128*rhou1[2]*uvar1[2]-0.3535533905932737*rhou1[0]*uvar1[2]-0.2258769757263128*rhou0[2]*uvar0[2]-0.3535533905932737*rhou0[0]*uvar0[2]-0.3535533905932737*uvar2[0]*rhou2[2]-0.3535533905932737*uvar1[0]*rhou1[2]-0.3535533905932737*uvar0[0]*rhou0[2]+energy[2]-0.3162277660168379*rhou2[1]*uvar2[1]-0.3162277660168379*rhou1[1]*uvar1[1]-0.3162277660168379*rhou0[1]*uvar0[1])*pressure_fac; 
} 
