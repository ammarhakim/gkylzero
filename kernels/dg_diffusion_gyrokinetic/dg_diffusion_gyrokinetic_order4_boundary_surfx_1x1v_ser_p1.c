#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_1x1v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],4.);

  double fSkin[6];
  fSkin[0] = 0.7071067811865476*(jacobgeo_inv[1]*qSkin[1]+jacobgeo_inv[0]*qSkin[0]); 
  fSkin[1] = 0.7071067811865476*(jacobgeo_inv[0]*qSkin[1]+qSkin[0]*jacobgeo_inv[1]); 
  fSkin[2] = 0.7071067811865476*(jacobgeo_inv[1]*qSkin[3]+jacobgeo_inv[0]*qSkin[2]); 
  fSkin[3] = 0.7071067811865476*(jacobgeo_inv[0]*qSkin[3]+jacobgeo_inv[1]*qSkin[2]); 
  fSkin[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qSkin[5]+21.21320343559643*jacobgeo_inv[0]*qSkin[4]); 
  fSkin[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[0]*qSkin[5]+21.21320343559643*jacobgeo_inv[1]*qSkin[4]); 

  double fEdge[6];
  fEdge[0] = 0.7071067811865476*(jacobgeo_inv[1]*qEdge[1]+jacobgeo_inv[0]*qEdge[0]); 
  fEdge[1] = 0.7071067811865476*(jacobgeo_inv[0]*qEdge[1]+qEdge[0]*jacobgeo_inv[1]); 
  fEdge[2] = 0.7071067811865476*(jacobgeo_inv[1]*qEdge[3]+jacobgeo_inv[0]*qEdge[2]); 
  fEdge[3] = 0.7071067811865476*(jacobgeo_inv[0]*qEdge[3]+jacobgeo_inv[1]*qEdge[2]); 
  fEdge[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qEdge[5]+21.21320343559643*jacobgeo_inv[0]*qEdge[4]); 
  fEdge[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[0]*qEdge[5]+21.21320343559643*jacobgeo_inv[1]*qEdge[4]); 

  double vol_incr[6] = {0.0}; 

  double edgeSurf_incr[6] = {0.0}; 
  double boundSurf_incr[6] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 1.623797632095822*coeff[0]*fSkin[1]+1.623797632095822*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]+1.623797632095822*coeff[0]*fSkin[0]-1.623797632095822*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 1.623797632095822*coeff[0]*fSkin[3]+1.623797632095822*coeff[0]*fEdge[3]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 3.5625*coeff[0]*fSkin[3]+2.0625*coeff[0]*fEdge[3]+1.623797632095822*coeff[0]*fSkin[2]-1.623797632095822*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = 1.623797632095823*coeff[0]*fSkin[5]+1.623797632095823*coeff[0]*fEdge[5]+0.9375*coeff[0]*fSkin[4]-0.9375*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*fSkin[5]+2.0625*coeff[0]*fEdge[5]+1.623797632095823*coeff[0]*fSkin[4]-1.623797632095823*coeff[0]*fEdge[4]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 1.5*coeff[0]*fSkin[3]; 
  boundSurf_incr[5] = 1.5*coeff[0]*fSkin[5]; 

  } else { 

  edgeSurf_incr[0] = (-1.623797632095822*coeff[0]*fSkin[1])-1.623797632095822*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]-1.623797632095822*coeff[0]*fSkin[0]+1.623797632095822*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-1.623797632095822*coeff[0]*fSkin[3])-1.623797632095822*coeff[0]*fEdge[3]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 3.5625*coeff[0]*fSkin[3]+2.0625*coeff[0]*fEdge[3]-1.623797632095822*coeff[0]*fSkin[2]+1.623797632095822*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = (-1.623797632095823*coeff[0]*fSkin[5])-1.623797632095823*coeff[0]*fEdge[5]+0.9375*coeff[0]*fSkin[4]-0.9375*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*fSkin[5]+2.0625*coeff[0]*fEdge[5]-1.623797632095823*coeff[0]*fSkin[4]+1.623797632095823*coeff[0]*fEdge[4]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 1.5*coeff[0]*fSkin[3]; 
  boundSurf_incr[5] = 1.5*coeff[0]*fSkin[5]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += -1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += -1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += -1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 

  }

  return 0.;
}

