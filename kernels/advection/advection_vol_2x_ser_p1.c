#include <gkyl_advection_kernels.h> 
GKYL_CU_DH double advection_vol_2x_ser_p1(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u[NDIM]:   Advection velocity.
  // f:         Input function.
  // out:       Incremented output.
  const double dx2 = 2.0/dxv[0]; 
  const double dy2 = 2.0/dxv[1]; 
  double alpha_mid = 0.0; 
  alpha_mid += fabs(0.25*u[0])/dxv[0]; 
  alpha_mid += fabs(0.25*u[4])/dxv[1]; 

  out[1] += 0.8660254037844386*(f[3]*u[3]+f[2]*u[2]+f[1]*u[1]+f[0]*u[0])*dx2; 
  out[2] += 0.8660254037844386*(f[3]*u[7]+f[2]*u[6]+f[1]*u[5]+f[0]*u[4])*dy2; 
  out[3] += 0.8660254037844386*((f[2]*u[7]+f[3]*u[6]+f[0]*u[5]+f[1]*u[4])*dy2+(f[1]*u[3]+u[1]*f[3]+f[0]*u[2]+u[0]*f[2])*dx2); 

  return alpha_mid; 
} 
