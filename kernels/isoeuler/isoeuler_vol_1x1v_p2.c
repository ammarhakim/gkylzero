#include <gkyl_isoeuler_kernels.h> 
GKYL_CU_DH double isoeuler_vol_1x1v_ser_p2(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // statevec: [rho, rho ux, rho uy, rho uz].
  // uvar: [ux, uy, uz].
  // out: Incremented output.

  const double *rho = &statevec[0]; 
  const double *rhou0 = &statevec[3]; 
  const double *uvar0 = &uvar[0]; 
  double *outrho = &out[0]; 
  double *outrhou0 = &out[3]; 
  double dx10 = 2./dxv[0]; 

  double vthsq = vth*vth; 
  outrho[0] += rho[2]*rhou0[2]*dx10+rho[1]*rhou0[1]*dx10+rho[0]*rhou0[0]*dx10; 

  outrhou0[0] += 1.732050807568877*rho[1]*vthsq+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 
  outrhou0[1] += 3.872983346207417*rho[2]*vthsq+rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 
  outrhou0[2] += rhou0[2]*uvar0[2]*dx10+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10; 

  return 0.; 
} 
