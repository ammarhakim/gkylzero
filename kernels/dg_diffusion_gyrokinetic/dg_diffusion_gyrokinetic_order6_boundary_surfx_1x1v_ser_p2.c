#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order6_boundary_surfx_1x1v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],6.);

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

  edgeSurf_incr[0] = (-35.21807064562169*coeff[0]*fSkin[4])+35.21807064562169*coeff[0]*fEdge[4]-34.09975027401226*coeff[0]*fSkin[1]-34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-70.53065765632411*coeff[0]*fSkin[4])+51.46831774920949*coeff[0]*fEdge[4]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]-34.09975027401226*coeff[0]*fSkin[0]+34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-35.21807064562168*coeff[0]*fSkin[6])+35.21807064562168*coeff[0]*fEdge[6]-34.09975027401226*coeff[0]*fSkin[3]-34.09975027401226*coeff[0]*fEdge[3]-19.6875*coeff[0]*fSkin[2]+19.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-70.53065765632414*coeff[0]*fSkin[6])+51.4683177492095*coeff[0]*fEdge[6]-61.5234375*coeff[0]*fSkin[3]-56.6015625*coeff[0]*fEdge[3]-34.09975027401226*coeff[0]*fSkin[2]+34.09975027401226*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = (-70.6640625*coeff[0]*fSkin[4])-3.1640625*coeff[0]*fEdge[4]-31.316701275974*coeff[0]*fSkin[1]-12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = (-34.09975027401227*coeff[0]*fSkin[7])-34.09975027401227*coeff[0]*fEdge[7]-19.6875*coeff[0]*fSkin[5]+19.6875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = (-70.6640625*coeff[0]*fSkin[6])-3.1640625*coeff[0]*fEdge[6]-31.31670127597403*coeff[0]*fSkin[3]-12.2543613688594*coeff[0]*fEdge[3]-12.57788237343632*coeff[0]*fSkin[2]+12.57788237343632*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = (-61.5234375*coeff[0]*fSkin[7])-56.6015625*coeff[0]*fEdge[7]-34.09975027401227*coeff[0]*fSkin[5]+34.09975027401227*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = 19.06233990711463*coeff[0]*fSkin[4]-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = 19.06233990711463*coeff[0]*fSkin[6]-4.921875*coeff[0]*fSkin[3]; 
  boundSurf_incr[4] = 19.06233990711463*coeff[0]*fSkin[1]-73.828125*coeff[0]*fSkin[4]; 
  boundSurf_incr[6] = 19.06233990711463*coeff[0]*fSkin[3]-73.828125*coeff[0]*fSkin[6]; 
  boundSurf_incr[7] = -4.921875*coeff[0]*fSkin[7]; 

  } else { 

  edgeSurf_incr[0] = (-35.21807064562169*coeff[0]*fSkin[4])+35.21807064562169*coeff[0]*fEdge[4]+34.09975027401226*coeff[0]*fSkin[1]+34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 70.53065765632411*coeff[0]*fSkin[4]-51.46831774920949*coeff[0]*fEdge[4]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]+34.09975027401226*coeff[0]*fSkin[0]-34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-35.21807064562168*coeff[0]*fSkin[6])+35.21807064562168*coeff[0]*fEdge[6]+34.09975027401226*coeff[0]*fSkin[3]+34.09975027401226*coeff[0]*fEdge[3]-19.6875*coeff[0]*fSkin[2]+19.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 70.53065765632414*coeff[0]*fSkin[6]-51.4683177492095*coeff[0]*fEdge[6]-61.5234375*coeff[0]*fSkin[3]-56.6015625*coeff[0]*fEdge[3]+34.09975027401226*coeff[0]*fSkin[2]-34.09975027401226*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = (-70.6640625*coeff[0]*fSkin[4])-3.1640625*coeff[0]*fEdge[4]+31.316701275974*coeff[0]*fSkin[1]+12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = 34.09975027401227*coeff[0]*fSkin[7]+34.09975027401227*coeff[0]*fEdge[7]-19.6875*coeff[0]*fSkin[5]+19.6875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = (-70.6640625*coeff[0]*fSkin[6])-3.1640625*coeff[0]*fEdge[6]+31.31670127597403*coeff[0]*fSkin[3]+12.2543613688594*coeff[0]*fEdge[3]-12.57788237343632*coeff[0]*fSkin[2]+12.57788237343632*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = (-61.5234375*coeff[0]*fSkin[7])-56.6015625*coeff[0]*fEdge[7]+34.09975027401227*coeff[0]*fSkin[5]-34.09975027401227*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = (-19.06233990711463*coeff[0]*fSkin[4])-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = (-19.06233990711463*coeff[0]*fSkin[6])-4.921875*coeff[0]*fSkin[3]; 
  boundSurf_incr[4] = (-73.828125*coeff[0]*fSkin[4])-19.06233990711463*coeff[0]*fSkin[1]; 
  boundSurf_incr[6] = (-73.828125*coeff[0]*fSkin[6])-19.06233990711463*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = -4.921875*coeff[0]*fSkin[7]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq; 

  }

  return 0.;
}

