#include <gkyl_dg_diffusion_kernels.h> 
GKYL_CU_DH double dg_diffusion_vol_1x_ser_p2(const double* w, const double* dx, double D, const double *q, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion coefficient
  // q: Input field
  // out: Incremented output

  const double dx0 = 2.0/dx[0]; 
  double J = 0.0; 
  J += pow(dx0, 2.0); 

  return 9.0*D*J; 
} 
