#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pressure_1x_ser_p1(const double gas_gamma, const double *uvar, const double *statevec, double* GKYL_RESTRICT pressure) 
{ 
  // gas_gamma: Adiabatic index.
  // uvar: [ux, uy, uz], Fluid flow.
  // statevec: [rho, rho ux, rho uy, rho uz, energy], Fluid input state vector.
  // pressure: Output pressure variable, p = (gas_gamma - 1)*(Energy - 0.5*rho*u^2).

  const double *rho = &statevec[0]; 
  const double *rhou0 = &statevec[2]; 
  const double *rhou1 = &statevec[4]; 
  const double *rhou2 = &statevec[6]; 
  const double *energy = &statevec[8]; 
  const double *uvar0 = &uvar[0]; 
  const double *uvar1 = &uvar[2]; 
  const double *uvar2 = &uvar[4]; 

  double pressure_fac = (gas_gamma-1.0); 

  pressure[0] = ((-0.3535533905932737*rhou2[1]*uvar2[1])-0.3535533905932737*rhou1[1]*uvar1[1]-0.3535533905932737*rhou0[1]*uvar0[1]-0.3535533905932737*rhou2[0]*uvar2[0]-0.3535533905932737*rhou1[0]*uvar1[0]-0.3535533905932737*rhou0[0]*uvar0[0]+energy[0])*pressure_fac; 
  pressure[1] = ((-0.3535533905932737*rhou2[0]*uvar2[1])-0.3535533905932737*rhou1[0]*uvar1[1]-0.3535533905932737*rhou0[0]*uvar0[1]-0.3535533905932737*uvar2[0]*rhou2[1]-0.3535533905932737*uvar1[0]*rhou1[1]-0.3535533905932737*uvar0[0]*rhou0[1]+energy[1])*pressure_fac; 
} 
