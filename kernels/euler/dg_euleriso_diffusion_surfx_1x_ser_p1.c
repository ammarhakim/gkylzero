#include <dg_euleriso_diffusion_kernels.h>

GKYL_CU_DH void
dg_euleriso_diffusion_surfx_1x_ser_p1(const double* w, const double* dx,
  const double* D_in,
  const double* uvarl, const double* uvarc, const double* uvarr,
  const double* statevecl, const double* statevecc, const double* statevecr,
  double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // uvarl: Input velocity in the left cell
  // uvarc: Input velocity in the center cell
  // uvarr: Input velocity in the right cell
  // statevecl: Input field in the left cell
  // statevecc: Input field in the center cell
  // statevecr: Input field in the right cell
  // out: Incremented output

  const double J = 4/dx[0]/dx[0];
  const double *D = &D_in[2]; 
  double mu = D[0]; 
  const double *uvarxl = &uvarl[0]; 
  const double *uvaryl = &uvarl[2]; 
  const double *uvarzl = &uvarl[4]; 
  const double *uvarxc = &uvarc[0]; 
  const double *uvaryc = &uvarc[2]; 
  const double *uvarzc = &uvarc[4]; 
  const double *uvarxr = &uvarr[0]; 
  const double *uvaryr = &uvarr[2]; 
  const double *uvarzr = &uvarr[4]; 
  out[2] += J*((-0.5412658773652741*uvarxr[1]*mu)+0.5412658773652741*uvarxl[1]*mu+0.5625*uvarxr[0]*mu+0.5625*uvarxl[0]*mu-1.125*uvarxc[0]*mu);
  out[3] += J*((-0.4375*uvarxr[1]*mu)-0.4375*uvarxl[1]*mu-2.875*uvarxc[1]*mu+0.5412658773652741*uvarxr[0]*mu-0.5412658773652741*uvarxl[0]*mu);
}
