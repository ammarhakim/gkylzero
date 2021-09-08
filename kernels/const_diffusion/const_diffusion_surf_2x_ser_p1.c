#include <gkyl_const_diffusion_kernels.h>

GKYL_CU_DH void
const_diffusion_surfx_2x_ser_p1(const double* w, const double* dx,
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

  out[0] += (-0.5412658773652741*fr[1]+0.5412658773652741*fl[1]+0.5625*fr[0]+0.5625*fl[0]-1.125*fc[0])*D[0]*J;
  out[1] += (-0.4375*fr[1]-0.4375*fl[1]-2.875*fc[1]+0.5412658773652741*fr[0]-0.5412658773652741*fl[0])*D[0]*J;
  out[2] += (-0.5412658773652741*fr[3]+0.5412658773652741*fl[3]+0.5625*fr[2]+0.5625*fl[2]-1.125*fc[2])*D[0]*J;
  out[3] += (-0.4375*fr[3]-0.4375*fl[3]-2.875*fc[3]+0.5412658773652741*fr[2]-0.5412658773652741*fl[2])*D[0]*J;
}

GKYL_CU_DH void
const_diffusion_surfy_2x_ser_p1(const double* w, const double* dx,
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

  const double J = 4/dx[1]/dx[1];

  out[0] += (-0.5412658773652741*fr[2]+0.5412658773652741*fl[2]+0.5625*fr[0]+0.5625*fl[0]-1.125*fc[0])*D[1]*J;
  out[1] += (-0.5412658773652741*fr[3]+0.5412658773652741*fl[3]+0.5625*fr[1]+0.5625*fl[1]-1.125*fc[1])*D[1]*J;
  out[2] += (-0.4375*fr[2]-0.4375*fl[2]-2.875*fc[2]+0.5412658773652741*fr[0]-0.5412658773652741*fl[0])*D[1]*J;
  out[3] += (-0.4375*fr[3]-0.4375*fl[3]-2.875*fc[3]+0.5412658773652741*fr[1]-0.5412658773652741*fl[1])*D[1]*J;
}

