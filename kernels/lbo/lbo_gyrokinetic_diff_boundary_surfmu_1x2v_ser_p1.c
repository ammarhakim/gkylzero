#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_boundary_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[4]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fskin/edge: Distribution function in cells 
  // out: Incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  double surfVar_l = w[2]-0.5*dxv[2];
  double surfVar_r = w[2]+0.5*dxv[2];

  double facDiff[2]; 
  // Expand diffusion coefficient in conf basis.
  facDiff[0] = 1.414213562373095*(bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0])*m_; 
  facDiff[1] = 1.414213562373095*(bmag_inv[0]*nuVtSqSum[1]+nuVtSqSum[0]*bmag_inv[1])*m_; 

  double bFacFskin[12] = {0.0}; 
  bFacFskin[3] = 0.8660254037844386*fskin[0]*dxv[2]; 
  bFacFskin[5] = 0.8660254037844386*fskin[1]*dxv[2]; 
  bFacFskin[6] = 0.8660254037844386*dxv[2]*fskin[2]; 
  bFacFskin[7] = 0.8660254037844386*dxv[2]*fskin[4]; 
  bFacFskin[10] = 0.8660254037844387*dxv[2]*fskin[8]; 
  bFacFskin[11] = 0.8660254037844387*dxv[2]*fskin[9]; 

  double vol_incr[12] = {0.0};
  vol_incr[3] = 0.7071067811865475*(facDiff[1]*bFacFskin[5]+facDiff[0]*bFacFskin[3]); 
  vol_incr[5] = 0.7071067811865475*(facDiff[0]*bFacFskin[5]+facDiff[1]*bFacFskin[3]); 
  vol_incr[6] = 0.7071067811865475*(facDiff[1]*bFacFskin[7]+facDiff[0]*bFacFskin[6]); 
  vol_incr[7] = 0.7071067811865475*(facDiff[0]*bFacFskin[7]+facDiff[1]*bFacFskin[6]); 
  vol_incr[10] = 0.7071067811865474*(facDiff[1]*bFacFskin[11]+facDiff[0]*bFacFskin[10]); 
  vol_incr[11] = 0.105409255338946*(6.708203932499369*facDiff[0]*bFacFskin[11]+6.708203932499369*facDiff[1]*bFacFskin[10]); 

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[12] = {0.0}; 
  edgeSurf[0] = (-0.5412658773652741*fskin[3]*surfVar_r)-0.5412658773652741*fedge[3]*surfVar_r-0.5625*fskin[0]*surfVar_r+0.5625*fedge[0]*surfVar_r; 
  edgeSurf[1] = (-0.5412658773652741*fskin[5]*surfVar_r)-0.5412658773652741*fedge[5]*surfVar_r-0.5625*fskin[1]*surfVar_r+0.5625*fedge[1]*surfVar_r; 
  edgeSurf[2] = (-0.5412658773652741*fskin[6]*surfVar_r)-0.5412658773652741*fedge[6]*surfVar_r-0.5625*fskin[2]*surfVar_r+0.5625*fedge[2]*surfVar_r; 
  edgeSurf[3] = (-1.4375*fskin[3]*surfVar_r)-0.4375*fedge[3]*surfVar_r-1.407291281149713*fskin[0]*surfVar_r+0.5412658773652741*fedge[0]*surfVar_r; 
  edgeSurf[4] = (-0.5412658773652741*fskin[7]*surfVar_r)-0.5412658773652741*fedge[7]*surfVar_r-0.5625*fskin[4]*surfVar_r+0.5625*fedge[4]*surfVar_r; 
  edgeSurf[5] = (-1.4375*fskin[5]*surfVar_r)-0.4375*fedge[5]*surfVar_r-1.407291281149713*fskin[1]*surfVar_r+0.5412658773652741*fedge[1]*surfVar_r; 
  edgeSurf[6] = (-1.4375*fskin[6]*surfVar_r)-0.4375*fedge[6]*surfVar_r-1.407291281149713*fskin[2]*surfVar_r+0.5412658773652741*fedge[2]*surfVar_r; 
  edgeSurf[7] = (-1.4375*fskin[7]*surfVar_r)-0.4375*fedge[7]*surfVar_r-1.407291281149713*fskin[4]*surfVar_r+0.5412658773652741*fedge[4]*surfVar_r; 
  edgeSurf[8] = (-0.5412658773652742*fskin[10]*surfVar_r)-0.5412658773652742*fedge[10]*surfVar_r-0.5625*fskin[8]*surfVar_r+0.5625*fedge[8]*surfVar_r; 
  edgeSurf[9] = (-0.5412658773652742*fskin[11]*surfVar_r)-0.5412658773652742*fedge[11]*surfVar_r-0.5625*fskin[9]*surfVar_r+0.5625*fedge[9]*surfVar_r; 
  edgeSurf[10] = (-1.4375*fskin[10]*surfVar_r)-0.4375*fedge[10]*surfVar_r-1.407291281149713*fskin[8]*surfVar_r+0.5412658773652742*fedge[8]*surfVar_r; 
  edgeSurf[11] = (-1.4375*fskin[11]*surfVar_r)-0.4375*fedge[11]*surfVar_r-1.407291281149713*fskin[9]*surfVar_r+0.5412658773652742*fedge[9]*surfVar_r; 

  double boundSurf[12] = {0.0}; 
  boundSurf[3] = 0.8660254037844386*fskin[0]*surfVar_l-1.5*fskin[3]*surfVar_l; 
  boundSurf[5] = 0.8660254037844386*fskin[1]*surfVar_l-1.5*fskin[5]*surfVar_l; 
  boundSurf[6] = 0.8660254037844386*fskin[2]*surfVar_l-1.5*fskin[6]*surfVar_l; 
  boundSurf[7] = 0.8660254037844386*fskin[4]*surfVar_l-1.5*fskin[7]*surfVar_l; 
  boundSurf[10] = 0.8660254037844387*fskin[8]*surfVar_l-1.5*fskin[10]*surfVar_l; 
  boundSurf[11] = 0.8660254037844387*fskin[9]*surfVar_l-1.5*fskin[11]*surfVar_l; 

  edgeSurf_incr[0] = 0.7071067811865475*(edgeSurf[1]*facDiff[1]+edgeSurf[0]*facDiff[0]); 
  edgeSurf_incr[1] = 0.7071067811865475*(edgeSurf[0]*facDiff[1]+facDiff[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.7071067811865475*(facDiff[1]*edgeSurf[4]+facDiff[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.7071067811865475*(facDiff[1]*edgeSurf[5]+facDiff[0]*edgeSurf[3]); 
  edgeSurf_incr[4] = 0.7071067811865475*(facDiff[0]*edgeSurf[4]+facDiff[1]*edgeSurf[2]); 
  edgeSurf_incr[5] = 0.7071067811865475*(facDiff[0]*edgeSurf[5]+facDiff[1]*edgeSurf[3]); 
  edgeSurf_incr[6] = 0.7071067811865475*(facDiff[1]*edgeSurf[7]+facDiff[0]*edgeSurf[6]); 
  edgeSurf_incr[7] = 0.7071067811865475*(facDiff[0]*edgeSurf[7]+facDiff[1]*edgeSurf[6]); 
  edgeSurf_incr[8] = 0.04714045207910316*(15.0*facDiff[1]*edgeSurf[9]+15.0*facDiff[0]*edgeSurf[8]); 
  edgeSurf_incr[9] = 0.04714045207910316*(15.0*facDiff[0]*edgeSurf[9]+15.0*facDiff[1]*edgeSurf[8]); 
  edgeSurf_incr[10] = 0.04714045207910316*(15.0*facDiff[1]*edgeSurf[11]+15.0*facDiff[0]*edgeSurf[10]); 
  edgeSurf_incr[11] = 0.04714045207910316*(15.0*facDiff[0]*edgeSurf[11]+15.0*facDiff[1]*edgeSurf[10]); 

  boundSurf_incr[0] = 0.7071067811865475*(boundSurf[1]*facDiff[1]+boundSurf[0]*facDiff[0]); 
  boundSurf_incr[1] = 0.7071067811865475*(boundSurf[0]*facDiff[1]+facDiff[0]*boundSurf[1]); 
  boundSurf_incr[2] = 0.7071067811865475*(facDiff[1]*boundSurf[4]+facDiff[0]*boundSurf[2]); 
  boundSurf_incr[3] = 0.7071067811865475*(facDiff[1]*boundSurf[5]+facDiff[0]*boundSurf[3]); 
  boundSurf_incr[4] = 0.7071067811865475*(facDiff[0]*boundSurf[4]+facDiff[1]*boundSurf[2]); 
  boundSurf_incr[5] = 0.7071067811865475*(facDiff[0]*boundSurf[5]+facDiff[1]*boundSurf[3]); 
  boundSurf_incr[6] = 0.7071067811865475*(facDiff[1]*boundSurf[7]+facDiff[0]*boundSurf[6]); 
  boundSurf_incr[7] = 0.7071067811865475*(facDiff[0]*boundSurf[7]+facDiff[1]*boundSurf[6]); 
  boundSurf_incr[8] = 0.04714045207910316*(15.0*facDiff[1]*boundSurf[9]+15.0*facDiff[0]*boundSurf[8]); 
  boundSurf_incr[9] = 0.04714045207910316*(15.0*facDiff[0]*boundSurf[9]+15.0*facDiff[1]*boundSurf[8]); 
  boundSurf_incr[10] = 0.04714045207910316*(15.0*facDiff[1]*boundSurf[11]+15.0*facDiff[0]*boundSurf[10]); 
  boundSurf_incr[11] = 0.04714045207910316*(15.0*facDiff[0]*boundSurf[11]+15.0*facDiff[1]*boundSurf[10]); 


  } else { 

  double edgeSurf[12] = {0.0}; 
  edgeSurf[0] = 0.0625*(8.660254037844386*(fskin[3]+fedge[3])-9.0*fskin[0]+9.0*fedge[0])*surfVar_l; 
  edgeSurf[1] = 0.0625*(8.660254037844386*(fskin[5]+fedge[5])-9.0*fskin[1]+9.0*fedge[1])*surfVar_l; 
  edgeSurf[2] = 0.0625*(8.660254037844386*(fskin[6]+fedge[6])-9.0*fskin[2]+9.0*fedge[2])*surfVar_l; 
  edgeSurf[3] = -0.0625*(23.0*fskin[3]+7.0*fedge[3]-22.5166604983954*fskin[0]+8.660254037844386*fedge[0])*surfVar_l; 
  edgeSurf[4] = 0.0625*(8.660254037844386*(fskin[7]+fedge[7])-9.0*fskin[4]+9.0*fedge[4])*surfVar_l; 
  edgeSurf[5] = -0.0625*(23.0*fskin[5]+7.0*fedge[5]-22.5166604983954*fskin[1]+8.660254037844386*fedge[1])*surfVar_l; 
  edgeSurf[6] = -0.0625*(23.0*fskin[6]+7.0*fedge[6]-22.5166604983954*fskin[2]+8.660254037844386*fedge[2])*surfVar_l; 
  edgeSurf[7] = -0.0625*(23.0*fskin[7]+7.0*fedge[7]-22.5166604983954*fskin[4]+8.660254037844386*fedge[4])*surfVar_l; 
  edgeSurf[8] = 0.0625*(8.660254037844387*(fskin[10]+fedge[10])-9.0*fskin[8]+9.0*fedge[8])*surfVar_l; 
  edgeSurf[9] = 0.0625*(8.660254037844387*(fskin[11]+fedge[11])-9.0*fskin[9]+9.0*fedge[9])*surfVar_l; 
  edgeSurf[10] = -0.0125*(115.0*fskin[10]+35.0*fedge[10]-112.583302491977*fskin[8]+43.30127018922195*fedge[8])*surfVar_l; 
  edgeSurf[11] = -0.0125*(115.0*fskin[11]+35.0*fedge[11]-112.583302491977*fskin[9]+43.30127018922195*fedge[9])*surfVar_l; 

  double boundSurf[12] = {0.0}; 
  boundSurf[3] = -0.5*(3.0*fskin[3]+1.732050807568877*fskin[0])*surfVar_r; 
  boundSurf[5] = -0.5*(3.0*fskin[5]+1.732050807568877*fskin[1])*surfVar_r; 
  boundSurf[6] = -0.5*(3.0*fskin[6]+1.732050807568877*fskin[2])*surfVar_r; 
  boundSurf[7] = -0.5*(3.0*fskin[7]+1.732050807568877*fskin[4])*surfVar_r; 
  boundSurf[10] = -0.1*(15.0*fskin[10]+8.660254037844387*fskin[8])*surfVar_r; 
  boundSurf[11] = -0.1*(15.0*fskin[11]+8.660254037844387*fskin[9])*surfVar_r; 

  edgeSurf_incr[0] = 0.7071067811865475*(edgeSurf[1]*facDiff[1]+edgeSurf[0]*facDiff[0]); 
  edgeSurf_incr[1] = 0.7071067811865475*(edgeSurf[0]*facDiff[1]+facDiff[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.7071067811865475*(facDiff[1]*edgeSurf[4]+facDiff[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.7071067811865475*(facDiff[1]*edgeSurf[5]+facDiff[0]*edgeSurf[3]); 
  edgeSurf_incr[4] = 0.7071067811865475*(facDiff[0]*edgeSurf[4]+facDiff[1]*edgeSurf[2]); 
  edgeSurf_incr[5] = 0.7071067811865475*(facDiff[0]*edgeSurf[5]+facDiff[1]*edgeSurf[3]); 
  edgeSurf_incr[6] = 0.7071067811865475*(facDiff[1]*edgeSurf[7]+facDiff[0]*edgeSurf[6]); 
  edgeSurf_incr[7] = 0.7071067811865475*(facDiff[0]*edgeSurf[7]+facDiff[1]*edgeSurf[6]); 
  edgeSurf_incr[8] = 0.04714045207910316*(15.0*facDiff[1]*edgeSurf[9]+15.0*facDiff[0]*edgeSurf[8]); 
  edgeSurf_incr[9] = 0.04714045207910316*(15.0*facDiff[0]*edgeSurf[9]+15.0*facDiff[1]*edgeSurf[8]); 
  edgeSurf_incr[10] = 0.04714045207910316*(15.0*facDiff[1]*edgeSurf[11]+15.0*facDiff[0]*edgeSurf[10]); 
  edgeSurf_incr[11] = 0.04714045207910316*(15.0*facDiff[0]*edgeSurf[11]+15.0*facDiff[1]*edgeSurf[10]); 

  boundSurf_incr[0] = 0.7071067811865475*(boundSurf[1]*facDiff[1]+boundSurf[0]*facDiff[0]); 
  boundSurf_incr[1] = 0.7071067811865475*(boundSurf[0]*facDiff[1]+facDiff[0]*boundSurf[1]); 
  boundSurf_incr[2] = 0.7071067811865475*(facDiff[1]*boundSurf[4]+facDiff[0]*boundSurf[2]); 
  boundSurf_incr[3] = 0.7071067811865475*(facDiff[1]*boundSurf[5]+facDiff[0]*boundSurf[3]); 
  boundSurf_incr[4] = 0.7071067811865475*(facDiff[0]*boundSurf[4]+facDiff[1]*boundSurf[2]); 
  boundSurf_incr[5] = 0.7071067811865475*(facDiff[0]*boundSurf[5]+facDiff[1]*boundSurf[3]); 
  boundSurf_incr[6] = 0.7071067811865475*(facDiff[1]*boundSurf[7]+facDiff[0]*boundSurf[6]); 
  boundSurf_incr[7] = 0.7071067811865475*(facDiff[0]*boundSurf[7]+facDiff[1]*boundSurf[6]); 
  boundSurf_incr[8] = 0.04714045207910316*(15.0*facDiff[1]*boundSurf[9]+15.0*facDiff[0]*boundSurf[8]); 
  boundSurf_incr[9] = 0.04714045207910316*(15.0*facDiff[0]*boundSurf[9]+15.0*facDiff[1]*boundSurf[8]); 
  boundSurf_incr[10] = 0.04714045207910316*(15.0*facDiff[1]*boundSurf[11]+15.0*facDiff[0]*boundSurf[10]); 
  boundSurf_incr[11] = 0.04714045207910316*(15.0*facDiff[0]*boundSurf[11]+15.0*facDiff[1]*boundSurf[10]); 

  } 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdvSq4; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdvSq4; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdvSq4; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdvSq4; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdvSq4; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdvSq4; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdvSq4; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdvSq4; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*rdvSq4; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*rdvSq4; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*rdvSq4; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*rdvSq4; 
} 
