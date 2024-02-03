#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order6_boundary_surfx_1x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],6.);

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

  edgeSurf_incr[0] = (-35.21807064562169*coeff[0]*fSkin[7])+35.21807064562169*coeff[0]*fEdge[7]-34.09975027401226*coeff[0]*fSkin[1]-34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-70.53065765632411*coeff[0]*fSkin[7])+51.46831774920949*coeff[0]*fEdge[7]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]-34.09975027401226*coeff[0]*fSkin[0]+34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-35.21807064562168*coeff[0]*fSkin[11])+35.21807064562168*coeff[0]*fEdge[11]-34.09975027401226*coeff[0]*fSkin[4]-34.09975027401226*coeff[0]*fEdge[4]-19.6875*coeff[0]*fSkin[2]+19.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-35.21807064562168*coeff[0]*fSkin[13])+35.21807064562168*coeff[0]*fEdge[13]-34.09975027401226*coeff[0]*fSkin[5]-34.09975027401226*coeff[0]*fEdge[5]-19.6875*coeff[0]*fSkin[3]+19.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-70.53065765632414*coeff[0]*fSkin[11])+51.4683177492095*coeff[0]*fEdge[11]-61.5234375*coeff[0]*fSkin[4]-56.6015625*coeff[0]*fEdge[4]-34.09975027401226*coeff[0]*fSkin[2]+34.09975027401226*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = (-70.53065765632414*coeff[0]*fSkin[13])+51.4683177492095*coeff[0]*fEdge[13]-61.5234375*coeff[0]*fSkin[5]-56.6015625*coeff[0]*fEdge[5]-34.09975027401226*coeff[0]*fSkin[3]+34.09975027401226*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-35.21807064562169*coeff[0]*fSkin[17])+35.21807064562169*coeff[0]*fEdge[17]-34.09975027401226*coeff[0]*fSkin[10]-34.09975027401226*coeff[0]*fEdge[10]-19.6875*coeff[0]*fSkin[6]+19.6875*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-70.6640625*coeff[0]*fSkin[7])-3.1640625*coeff[0]*fEdge[7]-31.316701275974*coeff[0]*fSkin[1]-12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = (-34.09975027401227*coeff[0]*fSkin[12])-34.09975027401227*coeff[0]*fEdge[12]-19.6875*coeff[0]*fSkin[8]+19.6875*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = (-34.09975027401227*coeff[0]*fSkin[15])-34.09975027401227*coeff[0]*fEdge[15]-19.6875*coeff[0]*fSkin[9]+19.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-70.53065765632411*coeff[0]*fSkin[17])+51.46831774920949*coeff[0]*fEdge[17]-61.5234375*coeff[0]*fSkin[10]-56.6015625*coeff[0]*fEdge[10]-34.09975027401226*coeff[0]*fSkin[6]+34.09975027401226*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = (-70.6640625*coeff[0]*fSkin[11])-3.1640625*coeff[0]*fEdge[11]-31.31670127597403*coeff[0]*fSkin[4]-12.2543613688594*coeff[0]*fEdge[4]-12.57788237343632*coeff[0]*fSkin[2]+12.57788237343632*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = (-61.5234375*coeff[0]*fSkin[12])-56.6015625*coeff[0]*fEdge[12]-34.09975027401227*coeff[0]*fSkin[8]+34.09975027401227*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = (-70.6640625*coeff[0]*fSkin[13])-3.1640625*coeff[0]*fEdge[13]-31.31670127597403*coeff[0]*fSkin[5]-12.2543613688594*coeff[0]*fEdge[5]-12.57788237343632*coeff[0]*fSkin[3]+12.57788237343632*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = (-34.09975027401227*coeff[0]*fSkin[18])-34.09975027401227*coeff[0]*fEdge[18]-19.6875*coeff[0]*fSkin[14]+19.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-61.5234375*coeff[0]*fSkin[15])-56.6015625*coeff[0]*fEdge[15]-34.09975027401227*coeff[0]*fSkin[9]+34.09975027401227*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = (-34.09975027401227*coeff[0]*fSkin[19])-34.09975027401227*coeff[0]*fEdge[19]-19.6875*coeff[0]*fSkin[16]+19.6875*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = (-70.6640625*coeff[0]*fSkin[17])-3.1640625*coeff[0]*fEdge[17]-31.316701275974*coeff[0]*fSkin[10]-12.2543613688594*coeff[0]*fEdge[10]-12.57788237343632*coeff[0]*fSkin[6]+12.57788237343632*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = (-61.5234375*coeff[0]*fSkin[18])-56.6015625*coeff[0]*fEdge[18]-34.09975027401227*coeff[0]*fSkin[14]+34.09975027401227*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = (-61.5234375*coeff[0]*fSkin[19])-56.6015625*coeff[0]*fEdge[19]-34.09975027401227*coeff[0]*fSkin[16]+34.09975027401227*coeff[0]*fEdge[16]; 

  boundSurf_incr[1] = 19.06233990711463*coeff[0]*fSkin[7]-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[4] = 19.06233990711463*coeff[0]*fSkin[11]-4.921875*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 19.06233990711463*coeff[0]*fSkin[13]-4.921875*coeff[0]*fSkin[5]; 
  boundSurf_incr[7] = 19.06233990711463*coeff[0]*fSkin[1]-73.828125*coeff[0]*fSkin[7]; 
  boundSurf_incr[10] = 19.06233990711463*coeff[0]*fSkin[17]-4.921875*coeff[0]*fSkin[10]; 
  boundSurf_incr[11] = 19.06233990711463*coeff[0]*fSkin[4]-73.828125*coeff[0]*fSkin[11]; 
  boundSurf_incr[12] = -4.921875*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 19.06233990711463*coeff[0]*fSkin[5]-73.828125*coeff[0]*fSkin[13]; 
  boundSurf_incr[15] = -4.921875*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = 19.06233990711463*coeff[0]*fSkin[10]-73.828125*coeff[0]*fSkin[17]; 
  boundSurf_incr[18] = -4.921875*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = -4.921875*coeff[0]*fSkin[19]; 

  } else { 

  edgeSurf_incr[0] = (-35.21807064562169*coeff[0]*fSkin[7])+35.21807064562169*coeff[0]*fEdge[7]+34.09975027401226*coeff[0]*fSkin[1]+34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 70.53065765632411*coeff[0]*fSkin[7]-51.46831774920949*coeff[0]*fEdge[7]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]+34.09975027401226*coeff[0]*fSkin[0]-34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-35.21807064562168*coeff[0]*fSkin[11])+35.21807064562168*coeff[0]*fEdge[11]+34.09975027401226*coeff[0]*fSkin[4]+34.09975027401226*coeff[0]*fEdge[4]-19.6875*coeff[0]*fSkin[2]+19.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-35.21807064562168*coeff[0]*fSkin[13])+35.21807064562168*coeff[0]*fEdge[13]+34.09975027401226*coeff[0]*fSkin[5]+34.09975027401226*coeff[0]*fEdge[5]-19.6875*coeff[0]*fSkin[3]+19.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 70.53065765632414*coeff[0]*fSkin[11]-51.4683177492095*coeff[0]*fEdge[11]-61.5234375*coeff[0]*fSkin[4]-56.6015625*coeff[0]*fEdge[4]+34.09975027401226*coeff[0]*fSkin[2]-34.09975027401226*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 70.53065765632414*coeff[0]*fSkin[13]-51.4683177492095*coeff[0]*fEdge[13]-61.5234375*coeff[0]*fSkin[5]-56.6015625*coeff[0]*fEdge[5]+34.09975027401226*coeff[0]*fSkin[3]-34.09975027401226*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-35.21807064562169*coeff[0]*fSkin[17])+35.21807064562169*coeff[0]*fEdge[17]+34.09975027401226*coeff[0]*fSkin[10]+34.09975027401226*coeff[0]*fEdge[10]-19.6875*coeff[0]*fSkin[6]+19.6875*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-70.6640625*coeff[0]*fSkin[7])-3.1640625*coeff[0]*fEdge[7]+31.316701275974*coeff[0]*fSkin[1]+12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = 34.09975027401227*coeff[0]*fSkin[12]+34.09975027401227*coeff[0]*fEdge[12]-19.6875*coeff[0]*fSkin[8]+19.6875*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 34.09975027401227*coeff[0]*fSkin[15]+34.09975027401227*coeff[0]*fEdge[15]-19.6875*coeff[0]*fSkin[9]+19.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 70.53065765632411*coeff[0]*fSkin[17]-51.46831774920949*coeff[0]*fEdge[17]-61.5234375*coeff[0]*fSkin[10]-56.6015625*coeff[0]*fEdge[10]+34.09975027401226*coeff[0]*fSkin[6]-34.09975027401226*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = (-70.6640625*coeff[0]*fSkin[11])-3.1640625*coeff[0]*fEdge[11]+31.31670127597403*coeff[0]*fSkin[4]+12.2543613688594*coeff[0]*fEdge[4]-12.57788237343632*coeff[0]*fSkin[2]+12.57788237343632*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = (-61.5234375*coeff[0]*fSkin[12])-56.6015625*coeff[0]*fEdge[12]+34.09975027401227*coeff[0]*fSkin[8]-34.09975027401227*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = (-70.6640625*coeff[0]*fSkin[13])-3.1640625*coeff[0]*fEdge[13]+31.31670127597403*coeff[0]*fSkin[5]+12.2543613688594*coeff[0]*fEdge[5]-12.57788237343632*coeff[0]*fSkin[3]+12.57788237343632*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = 34.09975027401227*coeff[0]*fSkin[18]+34.09975027401227*coeff[0]*fEdge[18]-19.6875*coeff[0]*fSkin[14]+19.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-61.5234375*coeff[0]*fSkin[15])-56.6015625*coeff[0]*fEdge[15]+34.09975027401227*coeff[0]*fSkin[9]-34.09975027401227*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = 34.09975027401227*coeff[0]*fSkin[19]+34.09975027401227*coeff[0]*fEdge[19]-19.6875*coeff[0]*fSkin[16]+19.6875*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = (-70.6640625*coeff[0]*fSkin[17])-3.1640625*coeff[0]*fEdge[17]+31.316701275974*coeff[0]*fSkin[10]+12.2543613688594*coeff[0]*fEdge[10]-12.57788237343632*coeff[0]*fSkin[6]+12.57788237343632*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = (-61.5234375*coeff[0]*fSkin[18])-56.6015625*coeff[0]*fEdge[18]+34.09975027401227*coeff[0]*fSkin[14]-34.09975027401227*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = (-61.5234375*coeff[0]*fSkin[19])-56.6015625*coeff[0]*fEdge[19]+34.09975027401227*coeff[0]*fSkin[16]-34.09975027401227*coeff[0]*fEdge[16]; 

  boundSurf_incr[1] = (-19.06233990711463*coeff[0]*fSkin[7])-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[4] = (-19.06233990711463*coeff[0]*fSkin[11])-4.921875*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = (-19.06233990711463*coeff[0]*fSkin[13])-4.921875*coeff[0]*fSkin[5]; 
  boundSurf_incr[7] = (-73.828125*coeff[0]*fSkin[7])-19.06233990711463*coeff[0]*fSkin[1]; 
  boundSurf_incr[10] = (-19.06233990711463*coeff[0]*fSkin[17])-4.921875*coeff[0]*fSkin[10]; 
  boundSurf_incr[11] = (-73.828125*coeff[0]*fSkin[11])-19.06233990711463*coeff[0]*fSkin[4]; 
  boundSurf_incr[12] = -4.921875*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = (-73.828125*coeff[0]*fSkin[13])-19.06233990711463*coeff[0]*fSkin[5]; 
  boundSurf_incr[15] = -4.921875*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = (-73.828125*coeff[0]*fSkin[17])-19.06233990711463*coeff[0]*fSkin[10]; 
  boundSurf_incr[18] = -4.921875*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = -4.921875*coeff[0]*fSkin[19]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*rdx2Sq; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*rdx2Sq; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*rdx2Sq; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*rdx2Sq; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*rdx2Sq; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*rdx2Sq; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*rdx2Sq; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*rdx2Sq; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*rdx2Sq; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*rdx2Sq; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*rdx2Sq; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*rdx2Sq; 

  }

  return 0.;
}

