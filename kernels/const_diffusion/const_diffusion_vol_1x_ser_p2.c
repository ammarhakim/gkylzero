#include <gkyl_const_diffusion_kernels.h>

GKYL_CU_DH double
const_diffusion_vol_1x_ser_p2(const double* w, const double* dx,
  const double* D, const double* f, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Constant diffusion coefficient
  // f: Input distribution function
  // out: Incremented output

  const double J = 4/dx[0]/dx[0];

  out[2] += 6.708203932499365*D[0]*f[0]*J;
  return D[0]*J;
}
