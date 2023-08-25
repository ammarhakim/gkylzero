#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order4_boundary_surfx_1x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[2] = {0.0}; 

  double edgeSurf_incr[2] = {0.0}; 
  double boundSurf_incr[2] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 1.623797632095822*coeff[0]*fSkin[1]+1.623797632095822*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]+1.623797632095822*coeff[0]*fSkin[0]-1.623797632095822*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 

  } else { 

  edgeSurf_incr[0] = (-1.623797632095822*coeff[0]*fSkin[1])-1.623797632095822*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]-1.623797632095822*coeff[0]*fSkin[0]+1.623797632095822*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_order4_boundary_surfx_1x_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[2] = {0.0}; 

  double edgeSurf_incr[2] = {0.0}; 
  double boundSurf_incr[2] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.5303300858899105*coeff[1]*fSkin[1])+1.148198316929614*coeff[0]*fSkin[1]+0.5303300858899105*coeff[1]*fEdge[1]+1.148198316929614*coeff[0]*fEdge[1]+0.6629126073623879*coeff[0]*fSkin[0]-0.6629126073623879*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-0.9185586535436913*coeff[1]*fSkin[1])+2.519067907977075*coeff[0]*fSkin[1]+0.9185586535436913*coeff[1]*fEdge[1]+1.458407736197253*coeff[0]*fEdge[1]+1.148198316929614*coeff[0]*fSkin[0]-1.148198316929614*coeff[0]*fEdge[0]; 

  boundSurf_incr[0] = -1.060660171779821*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = 1.837117307087383*coeff[1]*fSkin[1]+1.060660171779821*coeff[0]*fSkin[1]; 

  } else { 

  edgeSurf_incr[0] = (-0.5303300858899105*coeff[1]*fSkin[1])-1.148198316929614*coeff[0]*fSkin[1]+0.5303300858899105*coeff[1]*fEdge[1]-1.148198316929614*coeff[0]*fEdge[1]+0.6629126073623879*coeff[0]*fSkin[0]-0.6629126073623879*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.9185586535436913*coeff[1]*fSkin[1]+2.519067907977075*coeff[0]*fSkin[1]-0.9185586535436913*coeff[1]*fEdge[1]+1.458407736197253*coeff[0]*fEdge[1]-1.148198316929614*coeff[0]*fSkin[0]+1.148198316929614*coeff[0]*fEdge[0]; 

  boundSurf_incr[0] = -1.060660171779821*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = 1.060660171779821*coeff[0]*fSkin[1]-1.837117307087383*coeff[1]*fSkin[1]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 

  }

  return 0.;
}

