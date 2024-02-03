#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_1x1v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],4.);

  double fSkin[8];
  fSkin[0] = 0.7071067811865476*(jacobgeo_inv[2]*qSkin[4]+jacobgeo_inv[1]*qSkin[1]+jacobgeo_inv[0]*qSkin[0]); 
  fSkin[1] = 0.1*(6.324555320336761*(jacobgeo_inv[1]*qSkin[4]+qSkin[1]*jacobgeo_inv[2])+7.071067811865476*(jacobgeo_inv[0]*qSkin[1]+qSkin[0]*jacobgeo_inv[1])); 
  fSkin[2] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[2]*qSkin[6]+21.21320343559643*(jacobgeo_inv[1]*qSkin[3]+jacobgeo_inv[0]*qSkin[2])); 
  fSkin[3] = 0.03333333333333333*(18.97366596101028*jacobgeo_inv[1]*qSkin[6]+(18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*qSkin[3]+21.21320343559643*jacobgeo_inv[1]*qSkin[2]); 
  fSkin[4] = 0.01428571428571429*((31.62277660168381*jacobgeo_inv[2]+49.49747468305833*jacobgeo_inv[0])*qSkin[4]+49.49747468305833*qSkin[0]*jacobgeo_inv[2]+44.27188724235732*jacobgeo_inv[1]*qSkin[1]); 
  fSkin[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qSkin[7]+21.21320343559643*jacobgeo_inv[0]*qSkin[5]); 
  fSkin[6] = 0.004761904761904762*((94.86832980505142*jacobgeo_inv[2]+148.492424049175*jacobgeo_inv[0])*qSkin[6]+132.815661727072*jacobgeo_inv[1]*qSkin[3]+148.492424049175*jacobgeo_inv[2]*qSkin[2]); 
  fSkin[7] = 0.03333333333333333*((18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*qSkin[7]+21.21320343559643*jacobgeo_inv[1]*qSkin[5]); 

  double fEdge[8];
  fEdge[0] = 0.7071067811865476*(jacobgeo_inv[2]*qEdge[4]+jacobgeo_inv[1]*qEdge[1]+jacobgeo_inv[0]*qEdge[0]); 
  fEdge[1] = 0.1*(6.324555320336761*(jacobgeo_inv[1]*qEdge[4]+qEdge[1]*jacobgeo_inv[2])+7.071067811865476*(jacobgeo_inv[0]*qEdge[1]+qEdge[0]*jacobgeo_inv[1])); 
  fEdge[2] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[2]*qEdge[6]+21.21320343559643*(jacobgeo_inv[1]*qEdge[3]+jacobgeo_inv[0]*qEdge[2])); 
  fEdge[3] = 0.03333333333333333*(18.97366596101028*jacobgeo_inv[1]*qEdge[6]+(18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*qEdge[3]+21.21320343559643*jacobgeo_inv[1]*qEdge[2]); 
  fEdge[4] = 0.01428571428571429*((31.62277660168381*jacobgeo_inv[2]+49.49747468305833*jacobgeo_inv[0])*qEdge[4]+49.49747468305833*qEdge[0]*jacobgeo_inv[2]+44.27188724235732*jacobgeo_inv[1]*qEdge[1]); 
  fEdge[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qEdge[7]+21.21320343559643*jacobgeo_inv[0]*qEdge[5]); 
  fEdge[6] = 0.004761904761904762*((94.86832980505142*jacobgeo_inv[2]+148.492424049175*jacobgeo_inv[0])*qEdge[6]+132.815661727072*jacobgeo_inv[1]*qEdge[3]+148.492424049175*jacobgeo_inv[2]*qEdge[2]); 
  fEdge[7] = 0.03333333333333333*((18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*qEdge[7]+21.21320343559643*jacobgeo_inv[1]*qEdge[5]); 

  double vol_incr[8] = {0.0}; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[4]-6.708203932499369*coeff[0]*fEdge[4]+8.11898816047911*coeff[0]*fSkin[1]+8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 14.16059535957086*coeff[0]*fSkin[4]-9.077304717673634*coeff[0]*fEdge[4]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]+8.11898816047911*coeff[0]*fSkin[0]-8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.708203932499369*coeff[0]*fSkin[6]-6.708203932499369*coeff[0]*fEdge[6]+8.11898816047911*coeff[0]*fSkin[3]+8.11898816047911*coeff[0]*fEdge[3]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 14.16059535957087*coeff[0]*fSkin[6]-9.077304717673634*coeff[0]*fEdge[6]+15.46875*coeff[0]*fSkin[3]+12.65625*coeff[0]*fEdge[3]+8.11898816047911*coeff[0]*fSkin[2]-8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = 20.34375*coeff[0]*fSkin[4]-0.65625*coeff[0]*fEdge[4]+15.61296411439865*coeff[0]*fSkin[1]+4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = 8.118988160479114*coeff[0]*fSkin[7]+8.118988160479114*coeff[0]*fEdge[7]+4.6875*coeff[0]*fSkin[5]-4.6875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = 20.34375*coeff[0]*fSkin[6]-0.65625*coeff[0]*fEdge[6]+15.61296411439865*coeff[0]*fSkin[3]+4.72019845319029*coeff[0]*fEdge[3]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = 15.46875*coeff[0]*fSkin[7]+12.65625*coeff[0]*fEdge[7]+8.118988160479114*coeff[0]*fSkin[5]-8.118988160479114*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = 2.8125*coeff[0]*fSkin[1]-5.083290641897234*coeff[0]*fSkin[4]; 
  boundSurf_incr[3] = 2.8125*coeff[0]*fSkin[3]-5.083290641897235*coeff[0]*fSkin[6]; 
  boundSurf_incr[4] = 19.6875*coeff[0]*fSkin[4]-10.89276566120836*coeff[0]*fSkin[1]; 
  boundSurf_incr[6] = 19.6875*coeff[0]*fSkin[6]-10.89276566120836*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = 2.8125*coeff[0]*fSkin[7]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[4]-6.708203932499369*coeff[0]*fEdge[4]-8.11898816047911*coeff[0]*fSkin[1]-8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-14.16059535957086*coeff[0]*fSkin[4])+9.077304717673634*coeff[0]*fEdge[4]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]-8.11898816047911*coeff[0]*fSkin[0]+8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.708203932499369*coeff[0]*fSkin[6]-6.708203932499369*coeff[0]*fEdge[6]-8.11898816047911*coeff[0]*fSkin[3]-8.11898816047911*coeff[0]*fEdge[3]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-14.16059535957087*coeff[0]*fSkin[6])+9.077304717673634*coeff[0]*fEdge[6]+15.46875*coeff[0]*fSkin[3]+12.65625*coeff[0]*fEdge[3]-8.11898816047911*coeff[0]*fSkin[2]+8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = 20.34375*coeff[0]*fSkin[4]-0.65625*coeff[0]*fEdge[4]-15.61296411439865*coeff[0]*fSkin[1]-4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = (-8.118988160479114*coeff[0]*fSkin[7])-8.118988160479114*coeff[0]*fEdge[7]+4.6875*coeff[0]*fSkin[5]-4.6875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = 20.34375*coeff[0]*fSkin[6]-0.65625*coeff[0]*fEdge[6]-15.61296411439865*coeff[0]*fSkin[3]-4.72019845319029*coeff[0]*fEdge[3]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = 15.46875*coeff[0]*fSkin[7]+12.65625*coeff[0]*fEdge[7]-8.118988160479114*coeff[0]*fSkin[5]+8.118988160479114*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = 5.083290641897234*coeff[0]*fSkin[4]+2.8125*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 5.083290641897235*coeff[0]*fSkin[6]+2.8125*coeff[0]*fSkin[3]; 
  boundSurf_incr[4] = 19.6875*coeff[0]*fSkin[4]+10.89276566120836*coeff[0]*fSkin[1]; 
  boundSurf_incr[6] = 19.6875*coeff[0]*fSkin[6]+10.89276566120836*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = 2.8125*coeff[0]*fSkin[7]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += -1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += -1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += -1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 
  out[6] += -1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq; 
  out[7] += -1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq; 

  }

  return 0.;
}

