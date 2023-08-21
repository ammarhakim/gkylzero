#include <dg_euleriso_diffusion_kernels.h>

GKYL_CU_DH double
dg_euleriso_diffusion_surfy_2x_ser_p1(const double* w, const double* dx,
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

  const double J = 4/dx[1]/dx[1];
  const double *D = &D_in[4]; 
  double mu = D[0]; 
  const double *uvarxl = &uvarl[0]; 
  const double *uvaryl = &uvarl[4]; 
  const double *uvarzl = &uvarl[8]; 
  const double *uvarxc = &uvarc[0]; 
  const double *uvaryc = &uvarc[4]; 
  const double *uvarzc = &uvarc[8]; 
  const double *uvarxr = &uvarr[0]; 
  const double *uvaryr = &uvarr[4]; 
  const double *uvarzr = &uvarr[8]; 
  out[4] += J*((-0.5412658773652741*uvaryr[2]*mu)+0.5412658773652741*uvaryl[2]*mu-0.5412658773652741*uvarxr[1]*mu+0.5412658773652741*uvarxl[1]*mu+0.5625*uvaryr[0]*mu+0.5625*uvaryl[0]*mu-1.125*uvaryc[0]*mu+0.5625*uvarxr[0]*mu+0.5625*uvarxl[0]*mu-1.125*uvarxc[0]*mu);
  out[5] += J*((-0.5412658773652741*uvaryr[3]*mu)+0.5412658773652741*uvaryl[3]*mu+0.5625*uvaryr[1]*mu+0.5625*uvaryl[1]*mu-1.125*uvaryc[1]*mu-0.4375*uvarxr[1]*mu-0.4375*uvarxl[1]*mu-2.875*uvarxc[1]*mu+0.5412658773652741*uvarxr[0]*mu-0.5412658773652741*uvarxl[0]*mu);
  out[6] += J*((-0.5412658773652741*uvarxr[3]*mu)+0.5412658773652741*uvarxl[3]*mu-0.4375*uvaryr[2]*mu-0.4375*uvaryl[2]*mu-2.875*uvaryc[2]*mu+0.5625*uvarxr[2]*mu+0.5625*uvarxl[2]*mu-1.125*uvarxc[2]*mu+0.5412658773652741*uvaryr[0]*mu-0.5412658773652741*uvaryl[0]*mu);
  out[7] += J*((-0.4375*uvaryr[3]*mu)-0.4375*uvaryl[3]*mu-2.875*uvaryc[3]*mu-0.4375*uvarxr[3]*mu-0.4375*uvarxl[3]*mu-2.875*uvarxc[3]*mu+0.5412658773652741*uvarxr[2]*mu-0.5412658773652741*uvarxl[2]*mu+0.5412658773652741*uvaryr[1]*mu-0.5412658773652741*uvaryl[1]*mu);
  out[8] += J*((-0.5412658773652741*uvaryr[2]*mu)+0.5412658773652741*uvaryl[2]*mu-0.5412658773652741*uvarxr[1]*mu+0.5412658773652741*uvarxl[1]*mu+0.5625*uvaryr[0]*mu+0.5625*uvaryl[0]*mu-1.125*uvaryc[0]*mu+0.5625*uvarxr[0]*mu+0.5625*uvarxl[0]*mu-1.125*uvarxc[0]*mu);
  out[9] += J*((-0.5412658773652741*uvaryr[3]*mu)+0.5412658773652741*uvaryl[3]*mu+0.5625*uvaryr[1]*mu+0.5625*uvaryl[1]*mu-1.125*uvaryc[1]*mu-0.4375*uvarxr[1]*mu-0.4375*uvarxl[1]*mu-2.875*uvarxc[1]*mu+0.5412658773652741*uvarxr[0]*mu-0.5412658773652741*uvarxl[0]*mu);
  out[10] += J*((-0.5412658773652741*uvarxr[3]*mu)+0.5412658773652741*uvarxl[3]*mu-0.4375*uvaryr[2]*mu-0.4375*uvaryl[2]*mu-2.875*uvaryc[2]*mu+0.5625*uvarxr[2]*mu+0.5625*uvarxl[2]*mu-1.125*uvarxc[2]*mu+0.5412658773652741*uvaryr[0]*mu-0.5412658773652741*uvaryl[0]*mu);
  out[11] += J*((-0.4375*uvaryr[3]*mu)-0.4375*uvaryl[3]*mu-2.875*uvaryc[3]*mu-0.4375*uvarxr[3]*mu-0.4375*uvarxl[3]*mu-2.875*uvarxc[3]*mu+0.5412658773652741*uvarxr[2]*mu-0.5412658773652741*uvarxl[2]*mu+0.5412658773652741*uvaryr[1]*mu-0.5412658773652741*uvaryl[1]*mu);
  return 0.;

}
