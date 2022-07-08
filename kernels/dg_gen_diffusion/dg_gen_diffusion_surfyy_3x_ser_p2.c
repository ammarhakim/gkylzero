#include <gkyl_dg_gen_diffusion_kernels.h>

GKYL_CU_DH void
dg_gen_diffusion_surfyy_3x_ser_p2(const double* w, const double* dx,
  const double* Dij, const double* q[], double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // Dij: Diffusion coefficient in the center cell
  // q: Input field in the left cell
  // out: Incremented output

  const double Jyy = 4/dx[1]/dx[1];
  const double* Dyy = &Dij[60];

  const double* qclc = q[10];
  const double* qccc = q[13];
  const double* qcuc = q[16];

  out[0] += Jyy*(0.0);
  out[1] += Jyy*(0.0);
  out[2] += Jyy*(0.0);
  out[3] += Jyy*(0.0);
  out[4] += Jyy*(0.0);
  out[5] += Jyy*(0.0);
  out[6] += Jyy*(0.0);
  out[7] += Jyy*(0.0);
  out[8] += Jyy*(0.0);
  out[9] += Jyy*(0.0);
  out[10] += Jyy*(0.0);
  out[11] += Jyy*(0.0);
  out[12] += Jyy*(0.0);
  out[13] += Jyy*(0.0);
  out[14] += Jyy*(0.0);
  out[15] += Jyy*(0.0);
  out[16] += Jyy*(0.0);
  out[17] += Jyy*(0.0);
  out[18] += Jyy*(0.0);
  out[19] += Jyy*(0.0);
}
