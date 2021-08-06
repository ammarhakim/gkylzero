#include <gkyl_const_diffusion_kernels.h>

GKYL_CU_DH void
const_diffusion_surfx_1x_ser_p2(const double* w, const double* dx,
  const double* D, const double* fl, const double* fc, const double* fr,
  double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Constant diffusion coefficient
  // fl: Input distribution function in the left cell
  // fc: Input distribution function in the center cell
  // fr: Input distribution function in the right cell
  // out: Incremented output

  const double J = 4/dx[0]/dx[0];

  out[0] += (0.6708203932499367*fr[2]+0.6708203932499367*fl[2]-1.341640786499873*fc[2]-1.190784930203603*fr[1]+1.190784930203603*fl[1]+0.9374999999999994*fr[0]+0.9374999999999994*fl[0]-1.874999999999999*fc[0])*D[0]*J;
  out[1] += (0.738287450370789*fr[2]-0.738287450370789*fl[2]-1.453125*fr[1]-1.453125*fl[1]-5.343749999999997*fc[1]+1.190784930203602*fr[0]-1.190784930203602*fl[0])*D[0]*J;
  out[2] += (-0.1406249999999993*fr[2]-0.1406249999999993*fl[2]-6.281249999999997*fc[2]-0.3025768239224553*fr[1]+0.3025768239224553*fl[1]+0.4192627457812099*fr[0]+0.4192627457812099*fl[0]-7.546729424061787*fc[0])*D[0]*J;
}
