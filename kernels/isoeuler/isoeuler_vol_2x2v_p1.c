#include <gkyl_isoeuler_kernels.h> 
GKYL_CU_DH double isoeuler_vol_2x2v_ser_p1(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // statevec: [rho, rho ux, rho uy, rho uz].
  // uvar: [ux, uy, uz].
  // out: Incremented output.

  const double *rho = &statevec[0]; 
  const double *rhou0 = &statevec[4]; 
  const double *rhou1 = &statevec[8]; 
  const double *uvar0 = &uvar[0]; 
  const double *uvar1 = &uvar[4]; 
  double *outrho = &out[0]; 
  double *outrhou0 = &out[4]; 
  double *outrhou1 = &out[8]; 
  double dx10 = 2./dxv[0]; 
  double dx11 = 2./dxv[1]; 

  double vthsq = vth*vth; 
  outrho[0] += rho[3]*rhou1[3]*dx11+rho[2]*rhou1[2]*dx11+rho[1]*rhou1[1]*dx11+rho[0]*rhou1[0]*dx11+rho[3]*rhou0[3]*dx10+rho[2]*rhou0[2]*dx10+rho[1]*rhou0[1]*dx10+rho[0]*rhou0[0]*dx10; 

  outrhou0[0] += 1.732050807568877*rho[1]*vthsq+rhou0[3]*uvar1[3]*dx11+rhou0[2]*uvar1[2]*dx11+rhou0[1]*uvar1[1]*dx11+rhou0[0]*uvar1[0]*dx11+rhou0[3]*uvar0[3]*dx10+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 
  outrhou0[1] += rhou0[3]*uvar1[3]*dx11+rhou0[2]*uvar1[2]*dx11+rhou0[1]*uvar1[1]*dx11+rhou0[0]*uvar1[0]*dx11+rhou0[3]*uvar0[3]*dx10+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 
  outrhou0[2] += 1.732050807568877*rho[3]*vthsq+rhou0[3]*uvar1[3]*dx11+rhou0[2]*uvar1[2]*dx11+rhou0[1]*uvar1[1]*dx11+rhou0[0]*uvar1[0]*dx11+rhou0[3]*uvar0[3]*dx10+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 
  outrhou0[3] += rhou0[3]*uvar1[3]*dx11+rhou0[2]*uvar1[2]*dx11+rhou0[1]*uvar1[1]*dx11+rhou0[0]*uvar1[0]*dx11+rhou0[3]*uvar0[3]*dx10+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 

  outrhou1[0] += 1.732050807568877*rho[2]*vthsq+rhou1[3]*uvar1[3]*dx11+rhou1[2]*uvar1[2]*dx11+rhou1[1]*uvar1[1]*dx11+rhou1[0]*uvar1[0]*dx11+rhou1[3]*uvar0[3]*dx10+rhou1[2]*uvar0[2]*dx10+rhou1[1]*uvar0[1]*dx10+rhou1[0]*uvar0[0]*dx10; 
  outrhou1[1] += 1.732050807568877*rho[3]*vthsq+rhou1[3]*uvar1[3]*dx11+rhou1[2]*uvar1[2]*dx11+rhou1[1]*uvar1[1]*dx11+rhou1[0]*uvar1[0]*dx11+rhou1[3]*uvar0[3]*dx10+rhou1[2]*uvar0[2]*dx10+rhou1[1]*uvar0[1]*dx10+rhou1[0]*uvar0[0]*dx10; 
  outrhou1[2] += rhou1[3]*uvar1[3]*dx11+rhou1[2]*uvar1[2]*dx11+rhou1[1]*uvar1[1]*dx11+rhou1[0]*uvar1[0]*dx11+rhou1[3]*uvar0[3]*dx10+rhou1[2]*uvar0[2]*dx10+rhou1[1]*uvar0[1]*dx10+rhou1[0]*uvar0[0]*dx10; 
  outrhou1[3] += rhou1[3]*uvar1[3]*dx11+rhou1[2]*uvar1[2]*dx11+rhou1[1]*uvar1[1]*dx11+rhou1[0]*uvar1[0]*dx11+rhou1[3]*uvar0[3]*dx10+rhou1[2]*uvar0[2]*dx10+rhou1[1]*uvar0[1]*dx10+rhou1[0]*uvar0[0]*dx10; 

  return 0.; 
} 
