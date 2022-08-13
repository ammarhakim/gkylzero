#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_vol_1x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gas_gamma: Adiabatic index.
  // uvar: [ux, uy, uz], Fluid flow.
  // pvar: Fluid pressure.
  // statevec: [rho, rho ux, rho uy, rho uz, energy], Fluid input state vector.
  // out: Incremented output.

  const double *rho = &statevec[0]; 
  const double *rhou0 = &statevec[3]; 
  const double *rhou1 = &statevec[6]; 
  const double *rhou2 = &statevec[9]; 
  const double *energy = &statevec[12]; 
  const double *uvar0 = &uvar[0]; 
  const double *uvar1 = &uvar[3]; 
  const double *uvar2 = &uvar[6]; 
  double *outrho = &out[0]; 
  double *outrhou0 = &out[3]; 
  double *outrhou1 = &out[6]; 
  double *outrhou2 = &out[9]; 
  double *outenergy = &out[12]; 
  double dx10 = 2./dxv[0]; 

  double alpha_mid = 0.0; 
  alpha_mid += 0.5*dx10*(fabs(0.7071067811865475*uvar0[0]-0.7905694150420947*uvar0[2])+sqrt(((0.7071067811865475*pvar[0]-0.7905694150420947*pvar[2])*gas_gamma)/(0.7071067811865475*rho[0]-0.7905694150420947*rho[2]))); 

  outrho[1] += 1.732050807568877*rhou0[0]*dx10; 
  outrho[2] += 3.872983346207417*rhou0[1]*dx10; 

  outrhou0[1] += 1.224744871391589*rhou0[2]*uvar0[2]*dx10+1.224744871391589*rhou0[1]*uvar0[1]*dx10+1.224744871391589*rhou0[0]*uvar0[0]*dx10+1.732050807568877*pvar[0]*dx10; 
  outrhou0[2] += 2.449489742783178*rhou0[1]*uvar0[2]*dx10+2.449489742783178*uvar0[1]*rhou0[2]*dx10+2.738612787525831*rhou0[0]*uvar0[1]*dx10+2.738612787525831*uvar0[0]*rhou0[1]*dx10+3.872983346207417*pvar[1]*dx10; 

  outrhou1[1] += 1.224744871391589*rhou1[2]*uvar0[2]*dx10+1.224744871391589*rhou1[1]*uvar0[1]*dx10+1.224744871391589*rhou1[0]*uvar0[0]*dx10; 
  outrhou1[2] += 2.449489742783178*rhou1[1]*uvar0[2]*dx10+2.449489742783178*uvar0[1]*rhou1[2]*dx10+2.738612787525831*rhou1[0]*uvar0[1]*dx10+2.738612787525831*uvar0[0]*rhou1[1]*dx10; 

  outrhou2[1] += 1.224744871391589*rhou2[2]*uvar0[2]*dx10+1.224744871391589*rhou2[1]*uvar0[1]*dx10+1.224744871391589*rhou2[0]*uvar0[0]*dx10; 
  outrhou2[2] += 2.449489742783178*rhou2[1]*uvar0[2]*dx10+2.449489742783178*uvar0[1]*rhou2[2]*dx10+2.738612787525831*rhou2[0]*uvar0[1]*dx10+2.738612787525831*uvar0[0]*rhou2[1]*dx10; 

  outenergy[1] += 1.224744871391589*pvar[2]*uvar0[2]*dx10+1.224744871391589*energy[2]*uvar0[2]*dx10+1.224744871391589*pvar[1]*uvar0[1]*dx10+1.224744871391589*energy[1]*uvar0[1]*dx10+1.224744871391589*pvar[0]*uvar0[0]*dx10+1.224744871391589*energy[0]*uvar0[0]*dx10; 
  outenergy[2] += 2.449489742783178*pvar[1]*uvar0[2]*dx10+2.449489742783178*energy[1]*uvar0[2]*dx10+2.449489742783178*uvar0[1]*pvar[2]*dx10+2.449489742783178*uvar0[1]*energy[2]*dx10+2.738612787525831*pvar[0]*uvar0[1]*dx10+2.738612787525831*energy[0]*uvar0[1]*dx10+2.738612787525831*uvar0[0]*pvar[1]*dx10+2.738612787525831*uvar0[0]*energy[1]*dx10; 

  return alpha_mid; 
} 
