#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[4] = {0.0}; 

  double edgeSurf_incr[4] = {0.0}; 
  double boundSurf_incr[4] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 1.6237976320958223*coeff[0]*fSkin[1]+1.6237976320958223*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]+1.6237976320958223*coeff[0]*fSkin[0]-1.6237976320958223*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 1.6237976320958223*coeff[0]*fSkin[3]+1.6237976320958223*coeff[0]*fEdge[3]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 3.5625*coeff[0]*fSkin[3]+2.0625*coeff[0]*fEdge[3]+1.6237976320958223*coeff[0]*fSkin[2]-1.6237976320958223*coeff[0]*fEdge[2]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 1.5*coeff[0]*fSkin[3]; 

  } else { 

  edgeSurf_incr[0] = -(1.6237976320958223*coeff[0]*fSkin[1])-1.6237976320958223*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]-1.6237976320958223*coeff[0]*fSkin[0]+1.6237976320958223*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(1.6237976320958223*coeff[0]*fSkin[3])-1.6237976320958223*coeff[0]*fEdge[3]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 3.5625*coeff[0]*fSkin[3]+2.0625*coeff[0]*fEdge[3]-1.6237976320958223*coeff[0]*fSkin[2]+1.6237976320958223*coeff[0]*fEdge[2]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 1.5*coeff[0]*fSkin[3]; 

  }

  out[0] += -(1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac); 
  out[1] += -(1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac); 
  out[2] += -(1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac); 
  out[3] += -(1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac); 

  return 0.;
}

