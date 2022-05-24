#include <gkyl_diffusion_kernels.h>

GKYL_CU_DH void
diffusion_surfx_1x_ser_p1(const double* w, const double* dx,
  const double* Dl, const double* Dc, const double* Dr,
  const double* ql, const double* qc, const double* qr,
  double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Constant diffusion coefficient
  // ql: Input field in the left cell
  // qc: Input field function in the center cell
  // qr: Input field function in the right cell
  // out: Incremented output

  const double J = 4/dx[0]/dx[0];

  out[0] += (-0.5412658773652739*qr[1]+0.5412658773652739*ql[1]+0.5624999999999997*qr[0]+0.5624999999999997*ql[0]-1.124999999999999*qc[0])*Dc[0]*J;
  out[1] += (-0.4374999999999995*qr[1]-0.4374999999999995*ql[1]-2.874999999999999*qc[1]+0.5412658773652738*qr[0]-0.5412658773652738*ql[0])*Dc[0]*J;
}
