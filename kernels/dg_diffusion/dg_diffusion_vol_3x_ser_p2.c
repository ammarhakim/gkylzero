#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double
dg_diffusion_vol_3x_ser_p2(const double* w, const double* dx,
  const double* D, const double* q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion tensor
  // q: Input field
  // out: Incremented output

  const double J[3] = {4/dx[0]/dx[0], 4/dx[1]/dx[1], 4/dx[2]/dx[2]};


  return D[0]*J[0] + D[20]*J[1] + D[40]*J[2];
}
