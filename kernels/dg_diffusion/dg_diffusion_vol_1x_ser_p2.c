#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double
dg_diffusion_vol_1x_ser_p2(const double* w, const double* dx,
  const double* D, const double* q, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion tensor
  // q: Input field
  // out: Incremented output

  const double J = 4/dx[0]/dx[0];

  out[0] += (0.0)*J;
  out[1] += (4.743416490252569*q[1]*D[2]+2.121320343559642*q[0]*D[1])*J;
  out[2] += (14.23024947075771*D[2]*q[2]+10.60660171779821*q[0]*D[2]+9.48683298050514*D[1]*q[1]+4.743416490252569*D[0]*q[0])*J;

  return D[0]*J;
}
