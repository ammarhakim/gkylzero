#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_vol_2x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gas_gamma: Adiabatic index.
  // uvar: [ux, uy, uz], Fluid flow.
  // pvar: Fluid pressure.
  // statevec: [rho, rho ux, rho uy, rho uz, energy], Fluid input state vector.
  // out: Incremented output.

  const double *rho = &statevec[0]; 
  const double *rhou0 = &statevec[4]; 
  const double *rhou1 = &statevec[8]; 
  const double *rhou2 = &statevec[12]; 
  const double *energy = &statevec[16]; 
  const double *uvar0 = &uvar[0]; 
  const double *uvar1 = &uvar[4]; 
  const double *uvar2 = &uvar[8]; 
  double *outrho = &out[0]; 
  double *outrhou0 = &out[4]; 
  double *outrhou1 = &out[8]; 
  double *outrhou2 = &out[12]; 
  double *outenergy = &out[16]; 
  double dx10 = 2./dxv[0]; 
  double dx11 = 2./dxv[1]; 

  double alpha_mid = 0.0; 
  alpha_mid += 0.5*dx10*(fabs(0.5*uvar0[0])+sqrt((1.0*pvar[0]*gas_gamma)/rho[0])); 
  alpha_mid += 0.5*dx11*(fabs(0.5*uvar1[0])+sqrt((1.0*pvar[0]*gas_gamma)/rho[0])); 

  outrho[1] += 1.732050807568877*rhou0[0]*dx10; 
  outrho[2] += 1.732050807568877*rhou1[0]*dx11; 
  outrho[3] += 1.732050807568877*rhou1[1]*dx11+1.732050807568877*rhou0[2]*dx10; 

  outrhou0[1] += 0.8660254037844386*rhou0[3]*uvar0[3]*dx10+0.8660254037844386*rhou0[2]*uvar0[2]*dx10+0.8660254037844386*rhou0[1]*uvar0[1]*dx10+0.8660254037844386*rhou0[0]*uvar0[0]*dx10+1.732050807568877*pvar[0]*dx10; 
  outrhou0[2] += 0.8660254037844386*rhou0[3]*uvar1[3]*dx11+0.8660254037844386*rhou0[2]*uvar1[2]*dx11+0.8660254037844386*rhou0[1]*uvar1[1]*dx11+0.8660254037844386*rhou0[0]*uvar1[0]*dx11; 
  outrhou0[3] += 0.8660254037844386*rhou0[2]*uvar1[3]*dx11+0.8660254037844386*uvar1[2]*rhou0[3]*dx11+0.8660254037844386*rhou0[0]*uvar1[1]*dx11+0.8660254037844386*uvar1[0]*rhou0[1]*dx11+0.8660254037844386*rhou0[1]*uvar0[3]*dx10+0.8660254037844386*uvar0[1]*rhou0[3]*dx10+0.8660254037844386*rhou0[0]*uvar0[2]*dx10+0.8660254037844386*uvar0[0]*rhou0[2]*dx10+1.732050807568877*pvar[2]*dx10; 

  outrhou1[1] += 0.8660254037844386*rhou1[3]*uvar0[3]*dx10+0.8660254037844386*rhou1[2]*uvar0[2]*dx10+0.8660254037844386*rhou1[1]*uvar0[1]*dx10+0.8660254037844386*rhou1[0]*uvar0[0]*dx10; 
  outrhou1[2] += 0.8660254037844386*rhou1[3]*uvar1[3]*dx11+0.8660254037844386*rhou1[2]*uvar1[2]*dx11+0.8660254037844386*rhou1[1]*uvar1[1]*dx11+0.8660254037844386*rhou1[0]*uvar1[0]*dx11+1.732050807568877*pvar[0]*dx11; 
  outrhou1[3] += 0.8660254037844386*rhou1[2]*uvar1[3]*dx11+0.8660254037844386*uvar1[2]*rhou1[3]*dx11+0.8660254037844386*rhou1[0]*uvar1[1]*dx11+0.8660254037844386*uvar1[0]*rhou1[1]*dx11+1.732050807568877*pvar[1]*dx11+0.8660254037844386*rhou1[1]*uvar0[3]*dx10+0.8660254037844386*uvar0[1]*rhou1[3]*dx10+0.8660254037844386*rhou1[0]*uvar0[2]*dx10+0.8660254037844386*uvar0[0]*rhou1[2]*dx10; 

  outrhou2[1] += 0.8660254037844386*rhou2[3]*uvar0[3]*dx10+0.8660254037844386*rhou2[2]*uvar0[2]*dx10+0.8660254037844386*rhou2[1]*uvar0[1]*dx10+0.8660254037844386*rhou2[0]*uvar0[0]*dx10; 
  outrhou2[2] += 0.8660254037844386*rhou2[3]*uvar1[3]*dx11+0.8660254037844386*rhou2[2]*uvar1[2]*dx11+0.8660254037844386*rhou2[1]*uvar1[1]*dx11+0.8660254037844386*rhou2[0]*uvar1[0]*dx11; 
  outrhou2[3] += 0.8660254037844386*rhou2[2]*uvar1[3]*dx11+0.8660254037844386*uvar1[2]*rhou2[3]*dx11+0.8660254037844386*rhou2[0]*uvar1[1]*dx11+0.8660254037844386*uvar1[0]*rhou2[1]*dx11+0.8660254037844386*rhou2[1]*uvar0[3]*dx10+0.8660254037844386*uvar0[1]*rhou2[3]*dx10+0.8660254037844386*rhou2[0]*uvar0[2]*dx10+0.8660254037844386*uvar0[0]*rhou2[2]*dx10; 

  outenergy[1] += 0.8660254037844386*pvar[3]*uvar0[3]*dx10+0.8660254037844386*energy[3]*uvar0[3]*dx10+0.8660254037844386*pvar[2]*uvar0[2]*dx10+0.8660254037844386*energy[2]*uvar0[2]*dx10+0.8660254037844386*pvar[1]*uvar0[1]*dx10+0.8660254037844386*energy[1]*uvar0[1]*dx10+0.8660254037844386*pvar[0]*uvar0[0]*dx10+0.8660254037844386*energy[0]*uvar0[0]*dx10; 
  outenergy[2] += 0.8660254037844386*pvar[3]*uvar1[3]*dx11+0.8660254037844386*energy[3]*uvar1[3]*dx11+0.8660254037844386*pvar[2]*uvar1[2]*dx11+0.8660254037844386*energy[2]*uvar1[2]*dx11+0.8660254037844386*pvar[1]*uvar1[1]*dx11+0.8660254037844386*energy[1]*uvar1[1]*dx11+0.8660254037844386*pvar[0]*uvar1[0]*dx11+0.8660254037844386*energy[0]*uvar1[0]*dx11; 
  outenergy[3] += 0.8660254037844386*pvar[2]*uvar1[3]*dx11+0.8660254037844386*energy[2]*uvar1[3]*dx11+0.8660254037844386*uvar1[2]*pvar[3]*dx11+0.8660254037844386*uvar1[2]*energy[3]*dx11+0.8660254037844386*pvar[0]*uvar1[1]*dx11+0.8660254037844386*energy[0]*uvar1[1]*dx11+0.8660254037844386*uvar1[0]*pvar[1]*dx11+0.8660254037844386*uvar1[0]*energy[1]*dx11+0.8660254037844386*pvar[1]*uvar0[3]*dx10+0.8660254037844386*energy[1]*uvar0[3]*dx10+0.8660254037844386*uvar0[1]*pvar[3]*dx10+0.8660254037844386*uvar0[1]*energy[3]*dx10+0.8660254037844386*pvar[0]*uvar0[2]*dx10+0.8660254037844386*energy[0]*uvar0[2]*dx10+0.8660254037844386*uvar0[0]*pvar[2]*dx10+0.8660254037844386*uvar0[0]*energy[2]*dx10; 

  return alpha_mid; 
} 
