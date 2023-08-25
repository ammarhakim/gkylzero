#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order4_boundary_surfx_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
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

  edgeSurf_incr[0] = 1.623797632095822*coeff[0]*fSkin[1]+1.623797632095822*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]+1.623797632095822*coeff[0]*fSkin[0]-1.623797632095822*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 1.623797632095822*coeff[0]*fSkin[3]+1.623797632095822*coeff[0]*fEdge[3]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 3.5625*coeff[0]*fSkin[3]+2.0625*coeff[0]*fEdge[3]+1.623797632095822*coeff[0]*fSkin[2]-1.623797632095822*coeff[0]*fEdge[2]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 1.5*coeff[0]*fSkin[3]; 

  } else { 

  edgeSurf_incr[0] = (-1.623797632095822*coeff[0]*fSkin[1])-1.623797632095822*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]-1.623797632095822*coeff[0]*fSkin[0]+1.623797632095822*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-1.623797632095822*coeff[0]*fSkin[3])-1.623797632095822*coeff[0]*fEdge[3]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 3.5625*coeff[0]*fSkin[3]+2.0625*coeff[0]*fEdge[3]-1.623797632095822*coeff[0]*fSkin[2]+1.623797632095822*coeff[0]*fEdge[2]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 1.5*coeff[0]*fSkin[3]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += -1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_order4_boundary_surfx_2x_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
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

  edgeSurf_incr[0] = (-0.375*coeff[3]*fSkin[3])+0.8118988160479111*coeff[2]*fSkin[3]+0.375*coeff[3]*fEdge[3]+0.8118988160479111*coeff[2]*fEdge[3]+0.46875*coeff[2]*fSkin[2]-0.46875*coeff[2]*fEdge[2]-0.375*coeff[1]*fSkin[1]+0.8118988160479111*coeff[0]*fSkin[1]+0.375*coeff[1]*fEdge[1]+0.8118988160479111*coeff[0]*fEdge[1]+0.46875*coeff[0]*fSkin[0]-0.46875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-0.6495190528383289*coeff[3]*fSkin[3])+1.78125*coeff[2]*fSkin[3]+0.6495190528383289*coeff[3]*fEdge[3]+1.03125*coeff[2]*fEdge[3]+0.8118988160479111*coeff[2]*fSkin[2]-0.8118988160479111*coeff[2]*fEdge[2]-0.6495190528383289*coeff[1]*fSkin[1]+1.78125*coeff[0]*fSkin[1]+0.6495190528383289*coeff[1]*fEdge[1]+1.03125*coeff[0]*fEdge[1]+0.8118988160479111*coeff[0]*fSkin[0]-0.8118988160479111*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-0.375*coeff[1]*fSkin[3])+0.8118988160479111*coeff[0]*fSkin[3]+0.375*coeff[1]*fEdge[3]+0.8118988160479111*coeff[0]*fEdge[3]-0.375*fSkin[1]*coeff[3]+0.375*fEdge[1]*coeff[3]+0.46875*coeff[0]*fSkin[2]-0.46875*coeff[0]*fEdge[2]+0.8118988160479111*fSkin[1]*coeff[2]+0.8118988160479111*fEdge[1]*coeff[2]+0.46875*fSkin[0]*coeff[2]-0.46875*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = (-0.6495190528383289*coeff[1]*fSkin[3])+1.78125*coeff[0]*fSkin[3]+0.6495190528383289*coeff[1]*fEdge[3]+1.03125*coeff[0]*fEdge[3]-0.6495190528383289*fSkin[1]*coeff[3]+0.6495190528383289*fEdge[1]*coeff[3]+0.8118988160479111*coeff[0]*fSkin[2]-0.8118988160479111*coeff[0]*fEdge[2]+1.78125*fSkin[1]*coeff[2]+1.03125*fEdge[1]*coeff[2]+0.8118988160479111*fSkin[0]*coeff[2]-0.8118988160479111*fEdge[0]*coeff[2]; 

  boundSurf_incr[0] = (-0.75*coeff[3]*fSkin[3])-0.75*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = 1.299038105676658*coeff[3]*fSkin[3]+0.75*coeff[2]*fSkin[3]+1.299038105676658*coeff[1]*fSkin[1]+0.75*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = (-0.75*coeff[1]*fSkin[3])-0.75*fSkin[1]*coeff[3]; 
  boundSurf_incr[3] = 1.299038105676658*coeff[1]*fSkin[3]+0.75*coeff[0]*fSkin[3]+1.299038105676658*fSkin[1]*coeff[3]+0.75*fSkin[1]*coeff[2]; 

  } else { 

  edgeSurf_incr[0] = (-0.375*coeff[3]*fSkin[3])-0.8118988160479111*coeff[2]*fSkin[3]+0.375*coeff[3]*fEdge[3]-0.8118988160479111*coeff[2]*fEdge[3]+0.46875*coeff[2]*fSkin[2]-0.46875*coeff[2]*fEdge[2]-0.375*coeff[1]*fSkin[1]-0.8118988160479111*coeff[0]*fSkin[1]+0.375*coeff[1]*fEdge[1]-0.8118988160479111*coeff[0]*fEdge[1]+0.46875*coeff[0]*fSkin[0]-0.46875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.6495190528383289*coeff[3]*fSkin[3]+1.78125*coeff[2]*fSkin[3]-0.6495190528383289*coeff[3]*fEdge[3]+1.03125*coeff[2]*fEdge[3]-0.8118988160479111*coeff[2]*fSkin[2]+0.8118988160479111*coeff[2]*fEdge[2]+0.6495190528383289*coeff[1]*fSkin[1]+1.78125*coeff[0]*fSkin[1]-0.6495190528383289*coeff[1]*fEdge[1]+1.03125*coeff[0]*fEdge[1]-0.8118988160479111*coeff[0]*fSkin[0]+0.8118988160479111*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-0.375*coeff[1]*fSkin[3])-0.8118988160479111*coeff[0]*fSkin[3]+0.375*coeff[1]*fEdge[3]-0.8118988160479111*coeff[0]*fEdge[3]-0.375*fSkin[1]*coeff[3]+0.375*fEdge[1]*coeff[3]+0.46875*coeff[0]*fSkin[2]-0.46875*coeff[0]*fEdge[2]-0.8118988160479111*fSkin[1]*coeff[2]-0.8118988160479111*fEdge[1]*coeff[2]+0.46875*fSkin[0]*coeff[2]-0.46875*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = 0.6495190528383289*coeff[1]*fSkin[3]+1.78125*coeff[0]*fSkin[3]-0.6495190528383289*coeff[1]*fEdge[3]+1.03125*coeff[0]*fEdge[3]+0.6495190528383289*fSkin[1]*coeff[3]-0.6495190528383289*fEdge[1]*coeff[3]-0.8118988160479111*coeff[0]*fSkin[2]+0.8118988160479111*coeff[0]*fEdge[2]+1.78125*fSkin[1]*coeff[2]+1.03125*fEdge[1]*coeff[2]-0.8118988160479111*fSkin[0]*coeff[2]+0.8118988160479111*fEdge[0]*coeff[2]; 

  boundSurf_incr[0] = (-0.75*coeff[3]*fSkin[3])-0.75*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = (-1.299038105676658*coeff[3]*fSkin[3])+0.75*coeff[2]*fSkin[3]-1.299038105676658*coeff[1]*fSkin[1]+0.75*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = (-0.75*coeff[1]*fSkin[3])-0.75*fSkin[1]*coeff[3]; 
  boundSurf_incr[3] = (-1.299038105676658*coeff[1]*fSkin[3])+0.75*coeff[0]*fSkin[3]-1.299038105676658*fSkin[1]*coeff[3]+0.75*fSkin[1]*coeff[2]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += -1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 

  }

  return 0.;
}

