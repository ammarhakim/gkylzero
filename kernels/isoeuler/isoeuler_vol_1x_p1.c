#include <gkyl_isoeuler_kernels.h> 
GKYL_CU_DH double isoeuler_vol_1x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // statevec: [rho, rho ux, rho uy, rho uz].
  // uvar: [ux, uy, uz].
  // out: Incremented output.

  const double *rho = &statevec[0]; 
  const double *rhou0 = &statevec[2]; 
  const double *rhou1 = &statevec[4]; 
  const double *rhou2 = &statevec[6]; 
  const double *uvar0 = &uvar[0]; 
  const double *uvar1 = &uvar[2]; 
  const double *uvar2 = &uvar[4]; 
  double *outrho = &out[0]; 
  double *outrhou0 = &out[2]; 
  double *outrhou1 = &out[4]; 
  double *outrhou2 = &out[6]; 
  double dx10 = 2./dxv[0]; 

  double vthsq = vth*vth; 
  double alpha_mid = 0.0; 
  alpha_mid += 0.5*dx10*(fabs(0.7071067811865475*uvar0[0])+sqrt(vthsq)); 

  outrho[1] += 1.224744871391589*rho[1]*rhou0[1]*dx10+1.224744871391589*rho[0]*rhou0[0]*dx10; 

  outrhou0[1] += 3.0*rho[1]*dx10*vthsq+1.224744871391589*rhou0[1]*uvar0[1]*dx10+1.224744871391589*rhou0[0]*uvar0[0]*dx10; 

  return alpha_mid; 
} 
