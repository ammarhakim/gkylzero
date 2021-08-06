#include <gkyl_const_diffusion_kernels.h>

GKYL_CU_DH void
const_diffusion_surfx_1x_ser_p1(const double* w, const double* dx,
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

  out[0] += (-0.5412658773652739*fr[1]+0.5412658773652739*fl[1]+0.5624999999999997*fr[0]+0.5624999999999997*fl[0]-1.124999999999999*fc[0])*D[0]*J;
  out[1] += (-0.4374999999999995*fr[1]-0.4374999999999995*fl[1]-2.874999999999999*fc[1]+0.5412658773652738*fr[0]-0.5412658773652738*fl[0])*D[0]*J;
}
