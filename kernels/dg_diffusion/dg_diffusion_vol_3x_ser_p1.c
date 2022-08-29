#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double
dg_diffusion_vol_3x_ser_p1(const double* w, const double* dx,
  const double* D, const double* q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion tensor
  // q: Input field
  // out: Incremented output

  const double J[3] = {16.0/dx[0]/dx[0], 16.0/dx[1]/dx[1], 16.0/dx[2]/dx[2]};


  return D[0]*J[0] + D[8]*J[1] + D[16]*J[2];
}
