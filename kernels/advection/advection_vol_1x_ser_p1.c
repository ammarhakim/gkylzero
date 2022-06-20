#include <gkyl_advection_kernels.h> 
GKYL_CU_DH double advection_vol_1x_ser_p1(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u[NDIM]:   Advection velocity.
  // f:         Input function.
  // out:       Incremented output.
  const double dx2 = 2.0/dxv[0]; 
  double alpha_mid = 0.0; 
  alpha_mid += fabs(0.3535533905932737*u[0])/dxv[0]; 

  out[1] += 1.224744871391589*(f[1]*u[1]+f[0]*u[0])*dx2; 

  return alpha_mid; 
} 
