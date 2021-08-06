#include <gkyl_const_diffusion_kernels.h>

GKYL_CU_DH double
const_diffusion_vol_1x_ser_p1(const double* w, const double* dx,
  const double* D, const double* f, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Constant diffusion coefficient
  // f: Input distribution function
  // out: Incremented output

  // Volume term for p=1 does have no contribution.

  const double J = 4/dx[0]/dx[0];
  
  return D[0]*J;
}
