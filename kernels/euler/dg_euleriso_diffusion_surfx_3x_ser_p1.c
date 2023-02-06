#include <dg_euleriso_diffusion_kernels.h>

GKYL_CU_DH void
dg_euleriso_diffusion_surfx_3x_ser_p1(const double* w, const double* dx,
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
  const double *D = &D_in[8]; 
  double mu = D[0]; 
  const double *uvarxl = &uvarl[0]; 
  const double *uvaryl = &uvarl[8]; 
  const double *uvarzl = &uvarl[16]; 
  const double *uvarxc = &uvarc[0]; 
  const double *uvaryc = &uvarc[8]; 
  const double *uvarzc = &uvarc[16]; 
  const double *uvarxr = &uvarr[0]; 
  const double *uvaryr = &uvarr[8]; 
  const double *uvarzr = &uvarr[16]; 
  out[8] += J*((-0.5412658773652741*uvarzr[3]*mu)+0.5412658773652741*uvarzl[3]*mu-0.5412658773652741*uvaryr[2]*mu+0.5412658773652741*uvaryl[2]*mu-0.5412658773652741*uvarxr[1]*mu+0.5412658773652741*uvarxl[1]*mu+0.5625*uvarzr[0]*mu+0.5625*uvarzl[0]*mu-1.125*uvarzc[0]*mu+0.5625*uvaryr[0]*mu+0.5625*uvaryl[0]*mu-1.125*uvaryc[0]*mu+0.5625*uvarxr[0]*mu+0.5625*uvarxl[0]*mu-1.125*uvarxc[0]*mu);
  out[9] += J*((-0.5412658773652741*uvarzr[5]*mu)+0.5412658773652741*uvarzl[5]*mu-0.5412658773652741*uvaryr[4]*mu+0.5412658773652741*uvaryl[4]*mu+0.5625*uvarzr[1]*mu+0.5625*uvarzl[1]*mu-1.125*uvarzc[1]*mu+0.5625*uvaryr[1]*mu+0.5625*uvaryl[1]*mu-1.125*uvaryc[1]*mu-0.4375*uvarxr[1]*mu-0.4375*uvarxl[1]*mu-2.875*uvarxc[1]*mu+0.5412658773652741*uvarxr[0]*mu-0.5412658773652741*uvarxl[0]*mu);
  out[10] += J*((-0.5412658773652741*uvarzr[6]*mu)+0.5412658773652741*uvarzl[6]*mu-0.5412658773652741*uvarxr[4]*mu+0.5412658773652741*uvarxl[4]*mu+0.5625*uvarzr[2]*mu+0.5625*uvarzl[2]*mu-1.125*uvarzc[2]*mu-0.4375*uvaryr[2]*mu-0.4375*uvaryl[2]*mu-2.875*uvaryc[2]*mu+0.5625*uvarxr[2]*mu+0.5625*uvarxl[2]*mu-1.125*uvarxc[2]*mu+0.5412658773652741*uvaryr[0]*mu-0.5412658773652741*uvaryl[0]*mu);
  out[11] += J*((-0.5412658773652741*uvaryr[6]*mu)+0.5412658773652741*uvaryl[6]*mu-0.5412658773652741*uvarxr[5]*mu+0.5412658773652741*uvarxl[5]*mu-0.4375*uvarzr[3]*mu-0.4375*uvarzl[3]*mu-2.875*uvarzc[3]*mu+0.5625*uvaryr[3]*mu+0.5625*uvaryl[3]*mu-1.125*uvaryc[3]*mu+0.5625*uvarxr[3]*mu+0.5625*uvarxl[3]*mu-1.125*uvarxc[3]*mu+0.5412658773652741*uvarzr[0]*mu-0.5412658773652741*uvarzl[0]*mu);
  out[12] += J*((-0.5412658773652741*uvarzr[7]*mu)+0.5412658773652741*uvarzl[7]*mu+0.5625*uvarzr[4]*mu+0.5625*uvarzl[4]*mu-1.125*uvarzc[4]*mu-0.4375*uvaryr[4]*mu-0.4375*uvaryl[4]*mu-2.875*uvaryc[4]*mu-0.4375*uvarxr[4]*mu-0.4375*uvarxl[4]*mu-2.875*uvarxc[4]*mu+0.5412658773652741*uvarxr[2]*mu-0.5412658773652741*uvarxl[2]*mu+0.5412658773652741*uvaryr[1]*mu-0.5412658773652741*uvaryl[1]*mu);
  out[13] += J*((-0.5412658773652741*uvaryr[7]*mu)+0.5412658773652741*uvaryl[7]*mu-0.4375*uvarzr[5]*mu-0.4375*uvarzl[5]*mu-2.875*uvarzc[5]*mu+0.5625*uvaryr[5]*mu+0.5625*uvaryl[5]*mu-1.125*uvaryc[5]*mu-0.4375*uvarxr[5]*mu-0.4375*uvarxl[5]*mu-2.875*uvarxc[5]*mu+0.5412658773652741*uvarxr[3]*mu-0.5412658773652741*uvarxl[3]*mu+0.5412658773652741*uvarzr[1]*mu-0.5412658773652741*uvarzl[1]*mu);
  out[14] += J*((-0.5412658773652741*uvarxr[7]*mu)+0.5412658773652741*uvarxl[7]*mu-0.4375*uvarzr[6]*mu-0.4375*uvarzl[6]*mu-2.875*uvarzc[6]*mu-0.4375*uvaryr[6]*mu-0.4375*uvaryl[6]*mu-2.875*uvaryc[6]*mu+0.5625*uvarxr[6]*mu+0.5625*uvarxl[6]*mu-1.125*uvarxc[6]*mu+0.5412658773652741*uvaryr[3]*mu-0.5412658773652741*uvaryl[3]*mu+0.5412658773652741*uvarzr[2]*mu-0.5412658773652741*uvarzl[2]*mu);
  out[15] += J*((-0.4375*uvarzr[7]*mu)-0.4375*uvarzl[7]*mu-2.875*uvarzc[7]*mu-0.4375*uvaryr[7]*mu-0.4375*uvaryl[7]*mu-2.875*uvaryc[7]*mu-0.4375*uvarxr[7]*mu-0.4375*uvarxl[7]*mu-2.875*uvarxc[7]*mu+0.5412658773652741*uvarxr[6]*mu-0.5412658773652741*uvarxl[6]*mu+0.5412658773652741*uvaryr[5]*mu-0.5412658773652741*uvaryl[5]*mu+0.5412658773652741*uvarzr[4]*mu-0.5412658773652741*uvarzl[4]*mu);
  out[16] += J*((-0.5412658773652741*uvarzr[3]*mu)+0.5412658773652741*uvarzl[3]*mu-0.5412658773652741*uvaryr[2]*mu+0.5412658773652741*uvaryl[2]*mu-0.5412658773652741*uvarxr[1]*mu+0.5412658773652741*uvarxl[1]*mu+0.5625*uvarzr[0]*mu+0.5625*uvarzl[0]*mu-1.125*uvarzc[0]*mu+0.5625*uvaryr[0]*mu+0.5625*uvaryl[0]*mu-1.125*uvaryc[0]*mu+0.5625*uvarxr[0]*mu+0.5625*uvarxl[0]*mu-1.125*uvarxc[0]*mu);
  out[17] += J*((-0.5412658773652741*uvarzr[5]*mu)+0.5412658773652741*uvarzl[5]*mu-0.5412658773652741*uvaryr[4]*mu+0.5412658773652741*uvaryl[4]*mu+0.5625*uvarzr[1]*mu+0.5625*uvarzl[1]*mu-1.125*uvarzc[1]*mu+0.5625*uvaryr[1]*mu+0.5625*uvaryl[1]*mu-1.125*uvaryc[1]*mu-0.4375*uvarxr[1]*mu-0.4375*uvarxl[1]*mu-2.875*uvarxc[1]*mu+0.5412658773652741*uvarxr[0]*mu-0.5412658773652741*uvarxl[0]*mu);
  out[18] += J*((-0.5412658773652741*uvarzr[6]*mu)+0.5412658773652741*uvarzl[6]*mu-0.5412658773652741*uvarxr[4]*mu+0.5412658773652741*uvarxl[4]*mu+0.5625*uvarzr[2]*mu+0.5625*uvarzl[2]*mu-1.125*uvarzc[2]*mu-0.4375*uvaryr[2]*mu-0.4375*uvaryl[2]*mu-2.875*uvaryc[2]*mu+0.5625*uvarxr[2]*mu+0.5625*uvarxl[2]*mu-1.125*uvarxc[2]*mu+0.5412658773652741*uvaryr[0]*mu-0.5412658773652741*uvaryl[0]*mu);
  out[19] += J*((-0.5412658773652741*uvaryr[6]*mu)+0.5412658773652741*uvaryl[6]*mu-0.5412658773652741*uvarxr[5]*mu+0.5412658773652741*uvarxl[5]*mu-0.4375*uvarzr[3]*mu-0.4375*uvarzl[3]*mu-2.875*uvarzc[3]*mu+0.5625*uvaryr[3]*mu+0.5625*uvaryl[3]*mu-1.125*uvaryc[3]*mu+0.5625*uvarxr[3]*mu+0.5625*uvarxl[3]*mu-1.125*uvarxc[3]*mu+0.5412658773652741*uvarzr[0]*mu-0.5412658773652741*uvarzl[0]*mu);
  out[20] += J*((-0.5412658773652741*uvarzr[7]*mu)+0.5412658773652741*uvarzl[7]*mu+0.5625*uvarzr[4]*mu+0.5625*uvarzl[4]*mu-1.125*uvarzc[4]*mu-0.4375*uvaryr[4]*mu-0.4375*uvaryl[4]*mu-2.875*uvaryc[4]*mu-0.4375*uvarxr[4]*mu-0.4375*uvarxl[4]*mu-2.875*uvarxc[4]*mu+0.5412658773652741*uvarxr[2]*mu-0.5412658773652741*uvarxl[2]*mu+0.5412658773652741*uvaryr[1]*mu-0.5412658773652741*uvaryl[1]*mu);
  out[21] += J*((-0.5412658773652741*uvaryr[7]*mu)+0.5412658773652741*uvaryl[7]*mu-0.4375*uvarzr[5]*mu-0.4375*uvarzl[5]*mu-2.875*uvarzc[5]*mu+0.5625*uvaryr[5]*mu+0.5625*uvaryl[5]*mu-1.125*uvaryc[5]*mu-0.4375*uvarxr[5]*mu-0.4375*uvarxl[5]*mu-2.875*uvarxc[5]*mu+0.5412658773652741*uvarxr[3]*mu-0.5412658773652741*uvarxl[3]*mu+0.5412658773652741*uvarzr[1]*mu-0.5412658773652741*uvarzl[1]*mu);
  out[22] += J*((-0.5412658773652741*uvarxr[7]*mu)+0.5412658773652741*uvarxl[7]*mu-0.4375*uvarzr[6]*mu-0.4375*uvarzl[6]*mu-2.875*uvarzc[6]*mu-0.4375*uvaryr[6]*mu-0.4375*uvaryl[6]*mu-2.875*uvaryc[6]*mu+0.5625*uvarxr[6]*mu+0.5625*uvarxl[6]*mu-1.125*uvarxc[6]*mu+0.5412658773652741*uvaryr[3]*mu-0.5412658773652741*uvaryl[3]*mu+0.5412658773652741*uvarzr[2]*mu-0.5412658773652741*uvarzl[2]*mu);
  out[23] += J*((-0.4375*uvarzr[7]*mu)-0.4375*uvarzl[7]*mu-2.875*uvarzc[7]*mu-0.4375*uvaryr[7]*mu-0.4375*uvaryl[7]*mu-2.875*uvaryc[7]*mu-0.4375*uvarxr[7]*mu-0.4375*uvarxl[7]*mu-2.875*uvarxc[7]*mu+0.5412658773652741*uvarxr[6]*mu-0.5412658773652741*uvarxl[6]*mu+0.5412658773652741*uvaryr[5]*mu-0.5412658773652741*uvaryl[5]*mu+0.5412658773652741*uvarzr[4]*mu-0.5412658773652741*uvarzl[4]*mu);
  out[24] += J*((-0.5412658773652741*uvarzr[3]*mu)+0.5412658773652741*uvarzl[3]*mu-0.5412658773652741*uvaryr[2]*mu+0.5412658773652741*uvaryl[2]*mu-0.5412658773652741*uvarxr[1]*mu+0.5412658773652741*uvarxl[1]*mu+0.5625*uvarzr[0]*mu+0.5625*uvarzl[0]*mu-1.125*uvarzc[0]*mu+0.5625*uvaryr[0]*mu+0.5625*uvaryl[0]*mu-1.125*uvaryc[0]*mu+0.5625*uvarxr[0]*mu+0.5625*uvarxl[0]*mu-1.125*uvarxc[0]*mu);
  out[25] += J*((-0.5412658773652741*uvarzr[5]*mu)+0.5412658773652741*uvarzl[5]*mu-0.5412658773652741*uvaryr[4]*mu+0.5412658773652741*uvaryl[4]*mu+0.5625*uvarzr[1]*mu+0.5625*uvarzl[1]*mu-1.125*uvarzc[1]*mu+0.5625*uvaryr[1]*mu+0.5625*uvaryl[1]*mu-1.125*uvaryc[1]*mu-0.4375*uvarxr[1]*mu-0.4375*uvarxl[1]*mu-2.875*uvarxc[1]*mu+0.5412658773652741*uvarxr[0]*mu-0.5412658773652741*uvarxl[0]*mu);
  out[26] += J*((-0.5412658773652741*uvarzr[6]*mu)+0.5412658773652741*uvarzl[6]*mu-0.5412658773652741*uvarxr[4]*mu+0.5412658773652741*uvarxl[4]*mu+0.5625*uvarzr[2]*mu+0.5625*uvarzl[2]*mu-1.125*uvarzc[2]*mu-0.4375*uvaryr[2]*mu-0.4375*uvaryl[2]*mu-2.875*uvaryc[2]*mu+0.5625*uvarxr[2]*mu+0.5625*uvarxl[2]*mu-1.125*uvarxc[2]*mu+0.5412658773652741*uvaryr[0]*mu-0.5412658773652741*uvaryl[0]*mu);
  out[27] += J*((-0.5412658773652741*uvaryr[6]*mu)+0.5412658773652741*uvaryl[6]*mu-0.5412658773652741*uvarxr[5]*mu+0.5412658773652741*uvarxl[5]*mu-0.4375*uvarzr[3]*mu-0.4375*uvarzl[3]*mu-2.875*uvarzc[3]*mu+0.5625*uvaryr[3]*mu+0.5625*uvaryl[3]*mu-1.125*uvaryc[3]*mu+0.5625*uvarxr[3]*mu+0.5625*uvarxl[3]*mu-1.125*uvarxc[3]*mu+0.5412658773652741*uvarzr[0]*mu-0.5412658773652741*uvarzl[0]*mu);
  out[28] += J*((-0.5412658773652741*uvarzr[7]*mu)+0.5412658773652741*uvarzl[7]*mu+0.5625*uvarzr[4]*mu+0.5625*uvarzl[4]*mu-1.125*uvarzc[4]*mu-0.4375*uvaryr[4]*mu-0.4375*uvaryl[4]*mu-2.875*uvaryc[4]*mu-0.4375*uvarxr[4]*mu-0.4375*uvarxl[4]*mu-2.875*uvarxc[4]*mu+0.5412658773652741*uvarxr[2]*mu-0.5412658773652741*uvarxl[2]*mu+0.5412658773652741*uvaryr[1]*mu-0.5412658773652741*uvaryl[1]*mu);
  out[29] += J*((-0.5412658773652741*uvaryr[7]*mu)+0.5412658773652741*uvaryl[7]*mu-0.4375*uvarzr[5]*mu-0.4375*uvarzl[5]*mu-2.875*uvarzc[5]*mu+0.5625*uvaryr[5]*mu+0.5625*uvaryl[5]*mu-1.125*uvaryc[5]*mu-0.4375*uvarxr[5]*mu-0.4375*uvarxl[5]*mu-2.875*uvarxc[5]*mu+0.5412658773652741*uvarxr[3]*mu-0.5412658773652741*uvarxl[3]*mu+0.5412658773652741*uvarzr[1]*mu-0.5412658773652741*uvarzl[1]*mu);
  out[30] += J*((-0.5412658773652741*uvarxr[7]*mu)+0.5412658773652741*uvarxl[7]*mu-0.4375*uvarzr[6]*mu-0.4375*uvarzl[6]*mu-2.875*uvarzc[6]*mu-0.4375*uvaryr[6]*mu-0.4375*uvaryl[6]*mu-2.875*uvaryc[6]*mu+0.5625*uvarxr[6]*mu+0.5625*uvarxl[6]*mu-1.125*uvarxc[6]*mu+0.5412658773652741*uvaryr[3]*mu-0.5412658773652741*uvaryl[3]*mu+0.5412658773652741*uvarzr[2]*mu-0.5412658773652741*uvarzl[2]*mu);
  out[31] += J*((-0.4375*uvarzr[7]*mu)-0.4375*uvarzl[7]*mu-2.875*uvarzc[7]*mu-0.4375*uvaryr[7]*mu-0.4375*uvaryl[7]*mu-2.875*uvaryc[7]*mu-0.4375*uvarxr[7]*mu-0.4375*uvarxl[7]*mu-2.875*uvarxc[7]*mu+0.5412658773652741*uvarxr[6]*mu-0.5412658773652741*uvarxl[6]*mu+0.5412658773652741*uvaryr[5]*mu-0.5412658773652741*uvaryl[5]*mu+0.5412658773652741*uvarzr[4]*mu-0.5412658773652741*uvarzl[4]*mu);
}
