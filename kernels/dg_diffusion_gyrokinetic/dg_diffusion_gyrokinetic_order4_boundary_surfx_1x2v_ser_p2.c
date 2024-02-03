#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_1x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],4.);

  double fSkin[20];
  fSkin[0] = 0.7071067811865475*(jacobgeo_inv[2]*qSkin[7]+jacobgeo_inv[1]*qSkin[1]+jacobgeo_inv[0]*qSkin[0]); 
  fSkin[1] = 0.1414213562373095*(4.47213595499958*(jacobgeo_inv[1]*qSkin[7]+qSkin[1]*jacobgeo_inv[2])+5.0*(jacobgeo_inv[0]*qSkin[1]+qSkin[0]*jacobgeo_inv[1])); 
  fSkin[2] = 0.04714045207910316*(15.0*jacobgeo_inv[2]*qSkin[11]+15.0*(jacobgeo_inv[1]*qSkin[4]+jacobgeo_inv[0]*qSkin[2])); 
  fSkin[3] = 0.04714045207910316*(15.0*jacobgeo_inv[2]*qSkin[13]+15.0*(jacobgeo_inv[1]*qSkin[5]+jacobgeo_inv[0]*qSkin[3])); 
  fSkin[4] = 0.04714045207910316*(13.41640786499874*jacobgeo_inv[1]*qSkin[11]+(13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qSkin[4]+15.0*jacobgeo_inv[1]*qSkin[2]); 
  fSkin[5] = 0.04714045207910316*(13.41640786499874*jacobgeo_inv[1]*qSkin[13]+(13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qSkin[5]+15.0*jacobgeo_inv[1]*qSkin[3]); 
  fSkin[6] = 0.7071067811865475*(jacobgeo_inv[2]*qSkin[17]+jacobgeo_inv[1]*qSkin[10]+jacobgeo_inv[0]*qSkin[6]); 
  fSkin[7] = 0.02020305089104421*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*qSkin[7]+35.0*qSkin[0]*jacobgeo_inv[2]+31.30495168499706*jacobgeo_inv[1]*qSkin[1]); 
  fSkin[8] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qSkin[12]+15.0*jacobgeo_inv[0]*qSkin[8]); 
  fSkin[9] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qSkin[15]+15.0*jacobgeo_inv[0]*qSkin[9]); 
  fSkin[10] = 0.1414213562373095*(4.47213595499958*jacobgeo_inv[1]*qSkin[17]+(4.47213595499958*jacobgeo_inv[2]+5.0*jacobgeo_inv[0])*qSkin[10]+5.0*jacobgeo_inv[1]*qSkin[6]); 
  fSkin[11] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*qSkin[11]+93.91485505499116*jacobgeo_inv[1]*qSkin[4]+105.0*jacobgeo_inv[2]*qSkin[2]); 
  fSkin[12] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qSkin[12]+15.0*jacobgeo_inv[1]*qSkin[8]); 
  fSkin[13] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*qSkin[13]+93.91485505499116*jacobgeo_inv[1]*qSkin[5]+105.0*jacobgeo_inv[2]*qSkin[3]); 
  fSkin[14] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qSkin[18]+15.0*jacobgeo_inv[0]*qSkin[14]); 
  fSkin[15] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qSkin[15]+15.0*jacobgeo_inv[1]*qSkin[9]); 
  fSkin[16] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qSkin[19]+15.0*jacobgeo_inv[0]*qSkin[16]); 
  fSkin[17] = 0.02020305089104421*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*qSkin[17]+31.30495168499706*jacobgeo_inv[1]*qSkin[10]+35.0*jacobgeo_inv[2]*qSkin[6]); 
  fSkin[18] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qSkin[18]+15.0*jacobgeo_inv[1]*qSkin[14]); 
  fSkin[19] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qSkin[19]+15.0*jacobgeo_inv[1]*qSkin[16]); 

  double fEdge[20];
  fEdge[0] = 0.7071067811865475*(jacobgeo_inv[2]*qEdge[7]+jacobgeo_inv[1]*qEdge[1]+jacobgeo_inv[0]*qEdge[0]); 
  fEdge[1] = 0.1414213562373095*(4.47213595499958*(jacobgeo_inv[1]*qEdge[7]+qEdge[1]*jacobgeo_inv[2])+5.0*(jacobgeo_inv[0]*qEdge[1]+qEdge[0]*jacobgeo_inv[1])); 
  fEdge[2] = 0.04714045207910316*(15.0*jacobgeo_inv[2]*qEdge[11]+15.0*(jacobgeo_inv[1]*qEdge[4]+jacobgeo_inv[0]*qEdge[2])); 
  fEdge[3] = 0.04714045207910316*(15.0*jacobgeo_inv[2]*qEdge[13]+15.0*(jacobgeo_inv[1]*qEdge[5]+jacobgeo_inv[0]*qEdge[3])); 
  fEdge[4] = 0.04714045207910316*(13.41640786499874*jacobgeo_inv[1]*qEdge[11]+(13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qEdge[4]+15.0*jacobgeo_inv[1]*qEdge[2]); 
  fEdge[5] = 0.04714045207910316*(13.41640786499874*jacobgeo_inv[1]*qEdge[13]+(13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qEdge[5]+15.0*jacobgeo_inv[1]*qEdge[3]); 
  fEdge[6] = 0.7071067811865475*(jacobgeo_inv[2]*qEdge[17]+jacobgeo_inv[1]*qEdge[10]+jacobgeo_inv[0]*qEdge[6]); 
  fEdge[7] = 0.02020305089104421*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*qEdge[7]+35.0*qEdge[0]*jacobgeo_inv[2]+31.30495168499706*jacobgeo_inv[1]*qEdge[1]); 
  fEdge[8] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qEdge[12]+15.0*jacobgeo_inv[0]*qEdge[8]); 
  fEdge[9] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qEdge[15]+15.0*jacobgeo_inv[0]*qEdge[9]); 
  fEdge[10] = 0.1414213562373095*(4.47213595499958*jacobgeo_inv[1]*qEdge[17]+(4.47213595499958*jacobgeo_inv[2]+5.0*jacobgeo_inv[0])*qEdge[10]+5.0*jacobgeo_inv[1]*qEdge[6]); 
  fEdge[11] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*qEdge[11]+93.91485505499116*jacobgeo_inv[1]*qEdge[4]+105.0*jacobgeo_inv[2]*qEdge[2]); 
  fEdge[12] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qEdge[12]+15.0*jacobgeo_inv[1]*qEdge[8]); 
  fEdge[13] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*qEdge[13]+93.91485505499116*jacobgeo_inv[1]*qEdge[5]+105.0*jacobgeo_inv[2]*qEdge[3]); 
  fEdge[14] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qEdge[18]+15.0*jacobgeo_inv[0]*qEdge[14]); 
  fEdge[15] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qEdge[15]+15.0*jacobgeo_inv[1]*qEdge[9]); 
  fEdge[16] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qEdge[19]+15.0*jacobgeo_inv[0]*qEdge[16]); 
  fEdge[17] = 0.02020305089104421*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*qEdge[17]+31.30495168499706*jacobgeo_inv[1]*qEdge[10]+35.0*jacobgeo_inv[2]*qEdge[6]); 
  fEdge[18] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qEdge[18]+15.0*jacobgeo_inv[1]*qEdge[14]); 
  fEdge[19] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qEdge[19]+15.0*jacobgeo_inv[1]*qEdge[16]); 

  double vol_incr[20] = {0.0}; 

  double edgeSurf_incr[20] = {0.0}; 
  double boundSurf_incr[20] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[7]-6.708203932499369*coeff[0]*fEdge[7]+8.11898816047911*coeff[0]*fSkin[1]+8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 14.16059535957086*coeff[0]*fSkin[7]-9.077304717673634*coeff[0]*fEdge[7]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]+8.11898816047911*coeff[0]*fSkin[0]-8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.708203932499369*coeff[0]*fSkin[11]-6.708203932499369*coeff[0]*fEdge[11]+8.11898816047911*coeff[0]*fSkin[4]+8.11898816047911*coeff[0]*fEdge[4]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 6.708203932499369*coeff[0]*fSkin[13]-6.708203932499369*coeff[0]*fEdge[13]+8.11898816047911*coeff[0]*fSkin[5]+8.11898816047911*coeff[0]*fEdge[5]+4.6875*coeff[0]*fSkin[3]-4.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 14.16059535957087*coeff[0]*fSkin[11]-9.077304717673634*coeff[0]*fEdge[11]+15.46875*coeff[0]*fSkin[4]+12.65625*coeff[0]*fEdge[4]+8.11898816047911*coeff[0]*fSkin[2]-8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 14.16059535957087*coeff[0]*fSkin[13]-9.077304717673634*coeff[0]*fEdge[13]+15.46875*coeff[0]*fSkin[5]+12.65625*coeff[0]*fEdge[5]+8.11898816047911*coeff[0]*fSkin[3]-8.11898816047911*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 6.708203932499369*coeff[0]*fSkin[17]-6.708203932499369*coeff[0]*fEdge[17]+8.11898816047911*coeff[0]*fSkin[10]+8.11898816047911*coeff[0]*fEdge[10]+4.6875*coeff[0]*fSkin[6]-4.6875*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 20.34375*coeff[0]*fSkin[7]-0.65625*coeff[0]*fEdge[7]+15.61296411439865*coeff[0]*fSkin[1]+4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = 8.118988160479114*coeff[0]*fSkin[12]+8.118988160479114*coeff[0]*fEdge[12]+4.6875*coeff[0]*fSkin[8]-4.6875*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 8.118988160479114*coeff[0]*fSkin[15]+8.118988160479114*coeff[0]*fEdge[15]+4.6875*coeff[0]*fSkin[9]-4.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 14.16059535957086*coeff[0]*fSkin[17]-9.077304717673634*coeff[0]*fEdge[17]+15.46875*coeff[0]*fSkin[10]+12.65625*coeff[0]*fEdge[10]+8.11898816047911*coeff[0]*fSkin[6]-8.11898816047911*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*fSkin[11]-0.65625*coeff[0]*fEdge[11]+15.61296411439865*coeff[0]*fSkin[4]+4.72019845319029*coeff[0]*fEdge[4]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = 15.46875*coeff[0]*fSkin[12]+12.65625*coeff[0]*fEdge[12]+8.118988160479114*coeff[0]*fSkin[8]-8.118988160479114*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = 20.34375*coeff[0]*fSkin[13]-0.65625*coeff[0]*fEdge[13]+15.61296411439865*coeff[0]*fSkin[5]+4.72019845319029*coeff[0]*fEdge[5]+4.192627457812107*coeff[0]*fSkin[3]-4.192627457812107*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = 8.118988160479114*coeff[0]*fSkin[18]+8.118988160479114*coeff[0]*fEdge[18]+4.6875*coeff[0]*fSkin[14]-4.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 15.46875*coeff[0]*fSkin[15]+12.65625*coeff[0]*fEdge[15]+8.118988160479114*coeff[0]*fSkin[9]-8.118988160479114*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = 8.118988160479114*coeff[0]*fSkin[19]+8.118988160479114*coeff[0]*fEdge[19]+4.6875*coeff[0]*fSkin[16]-4.6875*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = 20.34375*coeff[0]*fSkin[17]-0.65625*coeff[0]*fEdge[17]+15.61296411439865*coeff[0]*fSkin[10]+4.720198453190292*coeff[0]*fEdge[10]+4.192627457812107*coeff[0]*fSkin[6]-4.192627457812107*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = 15.46875*coeff[0]*fSkin[18]+12.65625*coeff[0]*fEdge[18]+8.118988160479114*coeff[0]*fSkin[14]-8.118988160479114*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = 15.46875*coeff[0]*fSkin[19]+12.65625*coeff[0]*fEdge[19]+8.118988160479114*coeff[0]*fSkin[16]-8.118988160479114*coeff[0]*fEdge[16]; 

  boundSurf_incr[1] = 2.8125*coeff[0]*fSkin[1]-5.083290641897234*coeff[0]*fSkin[7]; 
  boundSurf_incr[4] = 2.8125*coeff[0]*fSkin[4]-5.083290641897235*coeff[0]*fSkin[11]; 
  boundSurf_incr[5] = 2.8125*coeff[0]*fSkin[5]-5.083290641897235*coeff[0]*fSkin[13]; 
  boundSurf_incr[7] = 19.6875*coeff[0]*fSkin[7]-10.89276566120836*coeff[0]*fSkin[1]; 
  boundSurf_incr[10] = 2.8125*coeff[0]*fSkin[10]-5.083290641897234*coeff[0]*fSkin[17]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*fSkin[11]-10.89276566120836*coeff[0]*fSkin[4]; 
  boundSurf_incr[12] = 2.8125*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 19.6875*coeff[0]*fSkin[13]-10.89276566120836*coeff[0]*fSkin[5]; 
  boundSurf_incr[15] = 2.8125*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = 19.6875*coeff[0]*fSkin[17]-10.89276566120836*coeff[0]*fSkin[10]; 
  boundSurf_incr[18] = 2.8125*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = 2.8125*coeff[0]*fSkin[19]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[7]-6.708203932499369*coeff[0]*fEdge[7]-8.11898816047911*coeff[0]*fSkin[1]-8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-14.16059535957086*coeff[0]*fSkin[7])+9.077304717673634*coeff[0]*fEdge[7]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]-8.11898816047911*coeff[0]*fSkin[0]+8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.708203932499369*coeff[0]*fSkin[11]-6.708203932499369*coeff[0]*fEdge[11]-8.11898816047911*coeff[0]*fSkin[4]-8.11898816047911*coeff[0]*fEdge[4]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 6.708203932499369*coeff[0]*fSkin[13]-6.708203932499369*coeff[0]*fEdge[13]-8.11898816047911*coeff[0]*fSkin[5]-8.11898816047911*coeff[0]*fEdge[5]+4.6875*coeff[0]*fSkin[3]-4.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-14.16059535957087*coeff[0]*fSkin[11])+9.077304717673634*coeff[0]*fEdge[11]+15.46875*coeff[0]*fSkin[4]+12.65625*coeff[0]*fEdge[4]-8.11898816047911*coeff[0]*fSkin[2]+8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = (-14.16059535957087*coeff[0]*fSkin[13])+9.077304717673634*coeff[0]*fEdge[13]+15.46875*coeff[0]*fSkin[5]+12.65625*coeff[0]*fEdge[5]-8.11898816047911*coeff[0]*fSkin[3]+8.11898816047911*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 6.708203932499369*coeff[0]*fSkin[17]-6.708203932499369*coeff[0]*fEdge[17]-8.11898816047911*coeff[0]*fSkin[10]-8.11898816047911*coeff[0]*fEdge[10]+4.6875*coeff[0]*fSkin[6]-4.6875*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 20.34375*coeff[0]*fSkin[7]-0.65625*coeff[0]*fEdge[7]-15.61296411439865*coeff[0]*fSkin[1]-4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = (-8.118988160479114*coeff[0]*fSkin[12])-8.118988160479114*coeff[0]*fEdge[12]+4.6875*coeff[0]*fSkin[8]-4.6875*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = (-8.118988160479114*coeff[0]*fSkin[15])-8.118988160479114*coeff[0]*fEdge[15]+4.6875*coeff[0]*fSkin[9]-4.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-14.16059535957086*coeff[0]*fSkin[17])+9.077304717673634*coeff[0]*fEdge[17]+15.46875*coeff[0]*fSkin[10]+12.65625*coeff[0]*fEdge[10]-8.11898816047911*coeff[0]*fSkin[6]+8.11898816047911*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*fSkin[11]-0.65625*coeff[0]*fEdge[11]-15.61296411439865*coeff[0]*fSkin[4]-4.72019845319029*coeff[0]*fEdge[4]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = 15.46875*coeff[0]*fSkin[12]+12.65625*coeff[0]*fEdge[12]-8.118988160479114*coeff[0]*fSkin[8]+8.118988160479114*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = 20.34375*coeff[0]*fSkin[13]-0.65625*coeff[0]*fEdge[13]-15.61296411439865*coeff[0]*fSkin[5]-4.72019845319029*coeff[0]*fEdge[5]+4.192627457812107*coeff[0]*fSkin[3]-4.192627457812107*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = (-8.118988160479114*coeff[0]*fSkin[18])-8.118988160479114*coeff[0]*fEdge[18]+4.6875*coeff[0]*fSkin[14]-4.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 15.46875*coeff[0]*fSkin[15]+12.65625*coeff[0]*fEdge[15]-8.118988160479114*coeff[0]*fSkin[9]+8.118988160479114*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = (-8.118988160479114*coeff[0]*fSkin[19])-8.118988160479114*coeff[0]*fEdge[19]+4.6875*coeff[0]*fSkin[16]-4.6875*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = 20.34375*coeff[0]*fSkin[17]-0.65625*coeff[0]*fEdge[17]-15.61296411439865*coeff[0]*fSkin[10]-4.720198453190292*coeff[0]*fEdge[10]+4.192627457812107*coeff[0]*fSkin[6]-4.192627457812107*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = 15.46875*coeff[0]*fSkin[18]+12.65625*coeff[0]*fEdge[18]-8.118988160479114*coeff[0]*fSkin[14]+8.118988160479114*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = 15.46875*coeff[0]*fSkin[19]+12.65625*coeff[0]*fEdge[19]-8.118988160479114*coeff[0]*fSkin[16]+8.118988160479114*coeff[0]*fEdge[16]; 

  boundSurf_incr[1] = 5.083290641897234*coeff[0]*fSkin[7]+2.8125*coeff[0]*fSkin[1]; 
  boundSurf_incr[4] = 5.083290641897235*coeff[0]*fSkin[11]+2.8125*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 5.083290641897235*coeff[0]*fSkin[13]+2.8125*coeff[0]*fSkin[5]; 
  boundSurf_incr[7] = 19.6875*coeff[0]*fSkin[7]+10.89276566120836*coeff[0]*fSkin[1]; 
  boundSurf_incr[10] = 5.083290641897234*coeff[0]*fSkin[17]+2.8125*coeff[0]*fSkin[10]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*fSkin[11]+10.89276566120836*coeff[0]*fSkin[4]; 
  boundSurf_incr[12] = 2.8125*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 19.6875*coeff[0]*fSkin[13]+10.89276566120836*coeff[0]*fSkin[5]; 
  boundSurf_incr[15] = 2.8125*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = 19.6875*coeff[0]*fSkin[17]+10.89276566120836*coeff[0]*fSkin[10]; 
  boundSurf_incr[18] = 2.8125*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = 2.8125*coeff[0]*fSkin[19]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += -1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += -1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += -1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 
  out[6] += -1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq; 
  out[7] += -1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq; 
  out[8] += -1.0*(vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*rdx2Sq; 
  out[9] += -1.0*(vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*rdx2Sq; 
  out[10] += -1.0*(vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*rdx2Sq; 
  out[11] += -1.0*(vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*rdx2Sq; 
  out[12] += -1.0*(vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*rdx2Sq; 
  out[13] += -1.0*(vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*rdx2Sq; 
  out[14] += -1.0*(vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*rdx2Sq; 
  out[15] += -1.0*(vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*rdx2Sq; 
  out[16] += -1.0*(vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*rdx2Sq; 
  out[17] += -1.0*(vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*rdx2Sq; 
  out[18] += -1.0*(vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*rdx2Sq; 
  out[19] += -1.0*(vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*rdx2Sq; 

  }

  return 0.;
}

