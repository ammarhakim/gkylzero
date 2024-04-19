#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfmu_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[4]: Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fskin/edge: Distribution function in cells 
  // out: Incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 

  double surfVar_l = w[3]-0.5*dxv[3];
  double surfVar_r = w[3]+0.5*dxv[3];

  double facDiff[4]; 
  // Expand diffusion coefficient in conf basis.
  facDiff[0] = (bmag_inv[3]*nuVtSqSum[3]+bmag_inv[2]*nuVtSqSum[2]+bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0])*m_; 
  facDiff[1] = (bmag_inv[2]*nuVtSqSum[3]+nuVtSqSum[2]*bmag_inv[3]+bmag_inv[0]*nuVtSqSum[1]+nuVtSqSum[0]*bmag_inv[1])*m_; 
  facDiff[2] = (bmag_inv[1]*nuVtSqSum[3]+nuVtSqSum[1]*bmag_inv[3]+bmag_inv[0]*nuVtSqSum[2]+nuVtSqSum[0]*bmag_inv[2])*m_; 
  facDiff[3] = (bmag_inv[0]*nuVtSqSum[3]+nuVtSqSum[0]*bmag_inv[3]+bmag_inv[1]*nuVtSqSum[2]+nuVtSqSum[1]*bmag_inv[2])*m_; 

  double bFacFskin[24] = {0.0}; 
  bFacFskin[4] = 0.8660254037844386*fskin[0]*dxv[3]; 
  bFacFskin[8] = 0.8660254037844386*fskin[1]*dxv[3]; 
  bFacFskin[9] = 0.8660254037844386*fskin[2]*dxv[3]; 
  bFacFskin[10] = 0.8660254037844386*dxv[3]*fskin[3]; 
  bFacFskin[12] = 0.8660254037844386*dxv[3]*fskin[5]; 
  bFacFskin[13] = 0.8660254037844386*dxv[3]*fskin[6]; 
  bFacFskin[14] = 0.8660254037844386*dxv[3]*fskin[7]; 
  bFacFskin[15] = 0.8660254037844386*dxv[3]*fskin[11]; 
  bFacFskin[19] = 0.8660254037844387*dxv[3]*fskin[16]; 
  bFacFskin[21] = 0.8660254037844387*dxv[3]*fskin[17]; 
  bFacFskin[22] = 0.8660254037844387*dxv[3]*fskin[18]; 
  bFacFskin[23] = 0.8660254037844387*dxv[3]*fskin[20]; 

  double vol_incr[24] = {0.0};
  vol_incr[4] = 0.5*(facDiff[3]*bFacFskin[12]+facDiff[2]*bFacFskin[9]+facDiff[1]*bFacFskin[8]+facDiff[0]*bFacFskin[4]); 
  vol_incr[8] = 0.5*(facDiff[2]*bFacFskin[12]+facDiff[3]*bFacFskin[9]+facDiff[0]*bFacFskin[8]+facDiff[1]*bFacFskin[4]); 
  vol_incr[9] = 0.5*(facDiff[1]*bFacFskin[12]+facDiff[0]*bFacFskin[9]+facDiff[3]*bFacFskin[8]+facDiff[2]*bFacFskin[4]); 
  vol_incr[10] = 0.5*(facDiff[3]*bFacFskin[15]+facDiff[2]*bFacFskin[14]+facDiff[1]*bFacFskin[13]+facDiff[0]*bFacFskin[10]); 
  vol_incr[12] = 0.5*(facDiff[0]*bFacFskin[12]+facDiff[1]*bFacFskin[9]+facDiff[2]*bFacFskin[8]+facDiff[3]*bFacFskin[4]); 
  vol_incr[13] = 0.5*(facDiff[2]*bFacFskin[15]+facDiff[3]*bFacFskin[14]+facDiff[0]*bFacFskin[13]+facDiff[1]*bFacFskin[10]); 
  vol_incr[14] = 0.5*(facDiff[1]*bFacFskin[15]+facDiff[0]*bFacFskin[14]+facDiff[3]*bFacFskin[13]+facDiff[2]*bFacFskin[10]); 
  vol_incr[15] = 0.5*(facDiff[0]*bFacFskin[15]+facDiff[1]*bFacFskin[14]+facDiff[2]*bFacFskin[13]+facDiff[3]*bFacFskin[10]); 
  vol_incr[19] = 0.4999999999999999*(facDiff[3]*bFacFskin[23]+facDiff[2]*bFacFskin[22]+facDiff[1]*bFacFskin[21]+facDiff[0]*bFacFskin[19]); 
  vol_incr[21] = 0.07453559924999298*(6.708203932499369*facDiff[2]*bFacFskin[23]+6.708203932499369*(facDiff[3]*bFacFskin[22]+facDiff[0]*bFacFskin[21])+6.708203932499369*facDiff[1]*bFacFskin[19]); 
  vol_incr[22] = 0.07453559924999298*(6.708203932499369*facDiff[1]*bFacFskin[23]+6.708203932499369*(facDiff[0]*bFacFskin[22]+facDiff[3]*bFacFskin[21])+6.708203932499369*facDiff[2]*bFacFskin[19]); 
  vol_incr[23] = 0.4999999999999999*(facDiff[0]*bFacFskin[23]+facDiff[1]*bFacFskin[22]+facDiff[2]*bFacFskin[21]+facDiff[3]*bFacFskin[19]); 

  double edgeSurf_incr[24] = {0.0}; 
  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[24] = {0.0}; 
  edgeSurf[0] = (-0.5412658773652741*fskin[4]*surfVar_r)-0.5412658773652741*fedge[4]*surfVar_r-0.5625*fskin[0]*surfVar_r+0.5625*fedge[0]*surfVar_r; 
  edgeSurf[1] = (-0.5412658773652741*fskin[8]*surfVar_r)-0.5412658773652741*fedge[8]*surfVar_r-0.5625*fskin[1]*surfVar_r+0.5625*fedge[1]*surfVar_r; 
  edgeSurf[2] = (-0.5412658773652741*fskin[9]*surfVar_r)-0.5412658773652741*fedge[9]*surfVar_r-0.5625*fskin[2]*surfVar_r+0.5625*fedge[2]*surfVar_r; 
  edgeSurf[3] = (-0.5412658773652741*fskin[10]*surfVar_r)-0.5412658773652741*fedge[10]*surfVar_r-0.5625*fskin[3]*surfVar_r+0.5625*fedge[3]*surfVar_r; 
  edgeSurf[4] = (-1.4375*fskin[4]*surfVar_r)-0.4375*fedge[4]*surfVar_r-1.407291281149713*fskin[0]*surfVar_r+0.5412658773652741*fedge[0]*surfVar_r; 
  edgeSurf[5] = (-0.5412658773652741*fskin[12]*surfVar_r)-0.5412658773652741*fedge[12]*surfVar_r-0.5625*fskin[5]*surfVar_r+0.5625*fedge[5]*surfVar_r; 
  edgeSurf[6] = (-0.5412658773652741*fskin[13]*surfVar_r)-0.5412658773652741*fedge[13]*surfVar_r-0.5625*fskin[6]*surfVar_r+0.5625*fedge[6]*surfVar_r; 
  edgeSurf[7] = (-0.5412658773652741*fskin[14]*surfVar_r)-0.5412658773652741*fedge[14]*surfVar_r-0.5625*fskin[7]*surfVar_r+0.5625*fedge[7]*surfVar_r; 
  edgeSurf[8] = (-1.4375*fskin[8]*surfVar_r)-0.4375*fedge[8]*surfVar_r-1.407291281149713*fskin[1]*surfVar_r+0.5412658773652741*fedge[1]*surfVar_r; 
  edgeSurf[9] = (-1.4375*fskin[9]*surfVar_r)-0.4375*fedge[9]*surfVar_r-1.407291281149713*fskin[2]*surfVar_r+0.5412658773652741*fedge[2]*surfVar_r; 
  edgeSurf[10] = (-1.4375*fskin[10]*surfVar_r)-0.4375*fedge[10]*surfVar_r-1.407291281149713*fskin[3]*surfVar_r+0.5412658773652741*fedge[3]*surfVar_r; 
  edgeSurf[11] = (-0.5412658773652741*fskin[15]*surfVar_r)-0.5412658773652741*fedge[15]*surfVar_r-0.5625*fskin[11]*surfVar_r+0.5625*fedge[11]*surfVar_r; 
  edgeSurf[12] = (-1.4375*fskin[12]*surfVar_r)-0.4375*fedge[12]*surfVar_r-1.407291281149713*fskin[5]*surfVar_r+0.5412658773652741*fedge[5]*surfVar_r; 
  edgeSurf[13] = (-1.4375*fskin[13]*surfVar_r)-0.4375*fedge[13]*surfVar_r-1.407291281149713*fskin[6]*surfVar_r+0.5412658773652741*fedge[6]*surfVar_r; 
  edgeSurf[14] = (-1.4375*fskin[14]*surfVar_r)-0.4375*fedge[14]*surfVar_r-1.407291281149713*fskin[7]*surfVar_r+0.5412658773652741*fedge[7]*surfVar_r; 
  edgeSurf[15] = (-1.4375*fskin[15]*surfVar_r)-0.4375*fedge[15]*surfVar_r-1.407291281149713*fskin[11]*surfVar_r+0.5412658773652741*fedge[11]*surfVar_r; 
  edgeSurf[16] = (-0.5412658773652742*fskin[19]*surfVar_r)-0.5412658773652742*fedge[19]*surfVar_r-0.5625*fskin[16]*surfVar_r+0.5625*fedge[16]*surfVar_r; 
  edgeSurf[17] = (-0.5412658773652742*fskin[21]*surfVar_r)-0.5412658773652742*fedge[21]*surfVar_r-0.5625*fskin[17]*surfVar_r+0.5625*fedge[17]*surfVar_r; 
  edgeSurf[18] = (-0.5412658773652742*fskin[22]*surfVar_r)-0.5412658773652742*fedge[22]*surfVar_r-0.5625*fskin[18]*surfVar_r+0.5625*fedge[18]*surfVar_r; 
  edgeSurf[19] = (-1.4375*fskin[19]*surfVar_r)-0.4375*fedge[19]*surfVar_r-1.407291281149713*fskin[16]*surfVar_r+0.5412658773652742*fedge[16]*surfVar_r; 
  edgeSurf[20] = (-0.5412658773652742*fskin[23]*surfVar_r)-0.5412658773652742*fedge[23]*surfVar_r-0.5625*fskin[20]*surfVar_r+0.5625*fedge[20]*surfVar_r; 
  edgeSurf[21] = (-1.4375*fskin[21]*surfVar_r)-0.4375*fedge[21]*surfVar_r-1.407291281149713*fskin[17]*surfVar_r+0.5412658773652742*fedge[17]*surfVar_r; 
  edgeSurf[22] = (-1.4375*fskin[22]*surfVar_r)-0.4375*fedge[22]*surfVar_r-1.407291281149713*fskin[18]*surfVar_r+0.5412658773652742*fedge[18]*surfVar_r; 
  edgeSurf[23] = (-1.4375*fskin[23]*surfVar_r)-0.4375*fedge[23]*surfVar_r-1.407291281149713*fskin[20]*surfVar_r+0.5412658773652742*fedge[20]*surfVar_r; 

  double boundSurf[24] = {0.0}; 
  boundSurf[4] = 0.8660254037844386*fskin[0]*surfVar_l-1.5*fskin[4]*surfVar_l; 
  boundSurf[8] = 0.8660254037844386*fskin[1]*surfVar_l-1.5*fskin[8]*surfVar_l; 
  boundSurf[9] = 0.8660254037844386*fskin[2]*surfVar_l-1.5*fskin[9]*surfVar_l; 
  boundSurf[10] = 0.8660254037844386*fskin[3]*surfVar_l-1.5*fskin[10]*surfVar_l; 
  boundSurf[12] = 0.8660254037844386*fskin[5]*surfVar_l-1.5*fskin[12]*surfVar_l; 
  boundSurf[13] = 0.8660254037844386*fskin[6]*surfVar_l-1.5*fskin[13]*surfVar_l; 
  boundSurf[14] = 0.8660254037844386*fskin[7]*surfVar_l-1.5*fskin[14]*surfVar_l; 
  boundSurf[15] = 0.8660254037844386*fskin[11]*surfVar_l-1.5*fskin[15]*surfVar_l; 
  boundSurf[19] = 0.8660254037844387*fskin[16]*surfVar_l-1.5*fskin[19]*surfVar_l; 
  boundSurf[21] = 0.8660254037844387*fskin[17]*surfVar_l-1.5*fskin[21]*surfVar_l; 
  boundSurf[22] = 0.8660254037844387*fskin[18]*surfVar_l-1.5*fskin[22]*surfVar_l; 
  boundSurf[23] = 0.8660254037844387*fskin[20]*surfVar_l-1.5*fskin[23]*surfVar_l; 

  edgeSurf_incr[0] = 0.5*(facDiff[3]*edgeSurf[5]+edgeSurf[2]*facDiff[2]+edgeSurf[1]*facDiff[1]+edgeSurf[0]*facDiff[0]); 
  edgeSurf_incr[1] = 0.5*(facDiff[2]*edgeSurf[5]+edgeSurf[2]*facDiff[3]+edgeSurf[0]*facDiff[1]+facDiff[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.5*(facDiff[1]*edgeSurf[5]+edgeSurf[1]*facDiff[3]+edgeSurf[0]*facDiff[2]+facDiff[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.5*(facDiff[3]*edgeSurf[11]+facDiff[2]*edgeSurf[7]+facDiff[1]*edgeSurf[6]+facDiff[0]*edgeSurf[3]); 
  edgeSurf_incr[4] = 0.5*(facDiff[3]*edgeSurf[12]+facDiff[2]*edgeSurf[9]+facDiff[1]*edgeSurf[8]+facDiff[0]*edgeSurf[4]); 
  edgeSurf_incr[5] = 0.5*(facDiff[0]*edgeSurf[5]+edgeSurf[0]*facDiff[3]+edgeSurf[1]*facDiff[2]+facDiff[1]*edgeSurf[2]); 
  edgeSurf_incr[6] = 0.5*(facDiff[2]*edgeSurf[11]+facDiff[3]*edgeSurf[7]+facDiff[0]*edgeSurf[6]+facDiff[1]*edgeSurf[3]); 
  edgeSurf_incr[7] = 0.5*(facDiff[1]*edgeSurf[11]+facDiff[0]*edgeSurf[7]+facDiff[3]*edgeSurf[6]+facDiff[2]*edgeSurf[3]); 
  edgeSurf_incr[8] = 0.5*(facDiff[2]*edgeSurf[12]+facDiff[3]*edgeSurf[9]+facDiff[0]*edgeSurf[8]+facDiff[1]*edgeSurf[4]); 
  edgeSurf_incr[9] = 0.5*(facDiff[1]*edgeSurf[12]+facDiff[0]*edgeSurf[9]+facDiff[3]*edgeSurf[8]+facDiff[2]*edgeSurf[4]); 
  edgeSurf_incr[10] = 0.5*(facDiff[3]*edgeSurf[15]+facDiff[2]*edgeSurf[14]+facDiff[1]*edgeSurf[13]+facDiff[0]*edgeSurf[10]); 
  edgeSurf_incr[11] = 0.5*(facDiff[0]*edgeSurf[11]+facDiff[1]*edgeSurf[7]+facDiff[2]*edgeSurf[6]+edgeSurf[3]*facDiff[3]); 
  edgeSurf_incr[12] = 0.5*(facDiff[0]*edgeSurf[12]+facDiff[1]*edgeSurf[9]+facDiff[2]*edgeSurf[8]+facDiff[3]*edgeSurf[4]); 
  edgeSurf_incr[13] = 0.5*(facDiff[2]*edgeSurf[15]+facDiff[3]*edgeSurf[14]+facDiff[0]*edgeSurf[13]+facDiff[1]*edgeSurf[10]); 
  edgeSurf_incr[14] = 0.5*(facDiff[1]*edgeSurf[15]+facDiff[0]*edgeSurf[14]+facDiff[3]*edgeSurf[13]+facDiff[2]*edgeSurf[10]); 
  edgeSurf_incr[15] = 0.5*(facDiff[0]*edgeSurf[15]+facDiff[1]*edgeSurf[14]+facDiff[2]*edgeSurf[13]+facDiff[3]*edgeSurf[10]); 
  edgeSurf_incr[16] = 0.03333333333333333*(15.0*facDiff[3]*edgeSurf[20]+15.0*(facDiff[2]*edgeSurf[18]+facDiff[1]*edgeSurf[17])+15.0*facDiff[0]*edgeSurf[16]); 
  edgeSurf_incr[17] = 0.03333333333333333*(15.0*facDiff[2]*edgeSurf[20]+15.0*(facDiff[3]*edgeSurf[18]+facDiff[0]*edgeSurf[17])+15.0*facDiff[1]*edgeSurf[16]); 
  edgeSurf_incr[18] = 0.03333333333333333*(15.0*facDiff[1]*edgeSurf[20]+15.0*(facDiff[0]*edgeSurf[18]+facDiff[3]*edgeSurf[17])+15.0*facDiff[2]*edgeSurf[16]); 
  edgeSurf_incr[19] = 0.03333333333333333*(15.0*facDiff[3]*edgeSurf[23]+15.0*(facDiff[2]*edgeSurf[22]+facDiff[1]*edgeSurf[21])+15.0*facDiff[0]*edgeSurf[19]); 
  edgeSurf_incr[20] = 0.03333333333333333*(15.0*facDiff[0]*edgeSurf[20]+15.0*(facDiff[1]*edgeSurf[18]+facDiff[2]*edgeSurf[17])+15.0*facDiff[3]*edgeSurf[16]); 
  edgeSurf_incr[21] = 0.03333333333333333*(15.0*facDiff[2]*edgeSurf[23]+15.0*(facDiff[3]*edgeSurf[22]+facDiff[0]*edgeSurf[21])+15.0*facDiff[1]*edgeSurf[19]); 
  edgeSurf_incr[22] = 0.03333333333333333*(15.0*facDiff[1]*edgeSurf[23]+15.0*(facDiff[0]*edgeSurf[22]+facDiff[3]*edgeSurf[21])+15.0*facDiff[2]*edgeSurf[19]); 
  edgeSurf_incr[23] = 0.03333333333333333*(15.0*facDiff[0]*edgeSurf[23]+15.0*(facDiff[1]*edgeSurf[22]+facDiff[2]*edgeSurf[21])+15.0*facDiff[3]*edgeSurf[19]); 

  boundSurf_incr[0] = 0.5*(facDiff[3]*boundSurf[5]+boundSurf[2]*facDiff[2]+boundSurf[1]*facDiff[1]+boundSurf[0]*facDiff[0]); 
  boundSurf_incr[1] = 0.5*(facDiff[2]*boundSurf[5]+boundSurf[2]*facDiff[3]+boundSurf[0]*facDiff[1]+facDiff[0]*boundSurf[1]); 
  boundSurf_incr[2] = 0.5*(facDiff[1]*boundSurf[5]+boundSurf[1]*facDiff[3]+boundSurf[0]*facDiff[2]+facDiff[0]*boundSurf[2]); 
  boundSurf_incr[3] = 0.5*(facDiff[3]*boundSurf[11]+facDiff[2]*boundSurf[7]+facDiff[1]*boundSurf[6]+facDiff[0]*boundSurf[3]); 
  boundSurf_incr[4] = 0.5*(facDiff[3]*boundSurf[12]+facDiff[2]*boundSurf[9]+facDiff[1]*boundSurf[8]+facDiff[0]*boundSurf[4]); 
  boundSurf_incr[5] = 0.5*(facDiff[0]*boundSurf[5]+boundSurf[0]*facDiff[3]+boundSurf[1]*facDiff[2]+facDiff[1]*boundSurf[2]); 
  boundSurf_incr[6] = 0.5*(facDiff[2]*boundSurf[11]+facDiff[3]*boundSurf[7]+facDiff[0]*boundSurf[6]+facDiff[1]*boundSurf[3]); 
  boundSurf_incr[7] = 0.5*(facDiff[1]*boundSurf[11]+facDiff[0]*boundSurf[7]+facDiff[3]*boundSurf[6]+facDiff[2]*boundSurf[3]); 
  boundSurf_incr[8] = 0.5*(facDiff[2]*boundSurf[12]+facDiff[3]*boundSurf[9]+facDiff[0]*boundSurf[8]+facDiff[1]*boundSurf[4]); 
  boundSurf_incr[9] = 0.5*(facDiff[1]*boundSurf[12]+facDiff[0]*boundSurf[9]+facDiff[3]*boundSurf[8]+facDiff[2]*boundSurf[4]); 
  boundSurf_incr[10] = 0.5*(facDiff[3]*boundSurf[15]+facDiff[2]*boundSurf[14]+facDiff[1]*boundSurf[13]+facDiff[0]*boundSurf[10]); 
  boundSurf_incr[11] = 0.5*(facDiff[0]*boundSurf[11]+facDiff[1]*boundSurf[7]+facDiff[2]*boundSurf[6]+boundSurf[3]*facDiff[3]); 
  boundSurf_incr[12] = 0.5*(facDiff[0]*boundSurf[12]+facDiff[1]*boundSurf[9]+facDiff[2]*boundSurf[8]+facDiff[3]*boundSurf[4]); 
  boundSurf_incr[13] = 0.5*(facDiff[2]*boundSurf[15]+facDiff[3]*boundSurf[14]+facDiff[0]*boundSurf[13]+facDiff[1]*boundSurf[10]); 
  boundSurf_incr[14] = 0.5*(facDiff[1]*boundSurf[15]+facDiff[0]*boundSurf[14]+facDiff[3]*boundSurf[13]+facDiff[2]*boundSurf[10]); 
  boundSurf_incr[15] = 0.5*(facDiff[0]*boundSurf[15]+facDiff[1]*boundSurf[14]+facDiff[2]*boundSurf[13]+facDiff[3]*boundSurf[10]); 
  boundSurf_incr[16] = 0.03333333333333333*(15.0*facDiff[3]*boundSurf[20]+15.0*(facDiff[2]*boundSurf[18]+facDiff[1]*boundSurf[17])+15.0*facDiff[0]*boundSurf[16]); 
  boundSurf_incr[17] = 0.03333333333333333*(15.0*facDiff[2]*boundSurf[20]+15.0*(facDiff[3]*boundSurf[18]+facDiff[0]*boundSurf[17])+15.0*facDiff[1]*boundSurf[16]); 
  boundSurf_incr[18] = 0.03333333333333333*(15.0*facDiff[1]*boundSurf[20]+15.0*(facDiff[0]*boundSurf[18]+facDiff[3]*boundSurf[17])+15.0*facDiff[2]*boundSurf[16]); 
  boundSurf_incr[19] = 0.03333333333333333*(15.0*facDiff[3]*boundSurf[23]+15.0*(facDiff[2]*boundSurf[22]+facDiff[1]*boundSurf[21])+15.0*facDiff[0]*boundSurf[19]); 
  boundSurf_incr[20] = 0.03333333333333333*(15.0*facDiff[0]*boundSurf[20]+15.0*(facDiff[1]*boundSurf[18]+facDiff[2]*boundSurf[17])+15.0*facDiff[3]*boundSurf[16]); 
  boundSurf_incr[21] = 0.03333333333333333*(15.0*facDiff[2]*boundSurf[23]+15.0*(facDiff[3]*boundSurf[22]+facDiff[0]*boundSurf[21])+15.0*facDiff[1]*boundSurf[19]); 
  boundSurf_incr[22] = 0.03333333333333333*(15.0*facDiff[1]*boundSurf[23]+15.0*(facDiff[0]*boundSurf[22]+facDiff[3]*boundSurf[21])+15.0*facDiff[2]*boundSurf[19]); 
  boundSurf_incr[23] = 0.03333333333333333*(15.0*facDiff[0]*boundSurf[23]+15.0*(facDiff[1]*boundSurf[22]+facDiff[2]*boundSurf[21])+15.0*facDiff[3]*boundSurf[19]); 


  } else { 

  double edgeSurf[24] = {0.0}; 
  edgeSurf[0] = 0.0625*(8.660254037844386*(fskin[4]+fedge[4])-9.0*fskin[0]+9.0*fedge[0])*surfVar_l; 
  edgeSurf[1] = 0.0625*(8.660254037844386*(fskin[8]+fedge[8])-9.0*fskin[1]+9.0*fedge[1])*surfVar_l; 
  edgeSurf[2] = 0.0625*(8.660254037844386*(fskin[9]+fedge[9])-9.0*fskin[2]+9.0*fedge[2])*surfVar_l; 
  edgeSurf[3] = 0.0625*(8.660254037844386*(fskin[10]+fedge[10])-9.0*fskin[3]+9.0*fedge[3])*surfVar_l; 
  edgeSurf[4] = -0.0625*(23.0*fskin[4]+7.0*fedge[4]-22.5166604983954*fskin[0]+8.660254037844386*fedge[0])*surfVar_l; 
  edgeSurf[5] = 0.0625*(8.660254037844386*(fskin[12]+fedge[12])-9.0*fskin[5]+9.0*fedge[5])*surfVar_l; 
  edgeSurf[6] = 0.0625*(8.660254037844386*(fskin[13]+fedge[13])-9.0*fskin[6]+9.0*fedge[6])*surfVar_l; 
  edgeSurf[7] = 0.0625*(8.660254037844386*(fskin[14]+fedge[14])-9.0*fskin[7]+9.0*fedge[7])*surfVar_l; 
  edgeSurf[8] = -0.0625*(23.0*fskin[8]+7.0*fedge[8]-22.5166604983954*fskin[1]+8.660254037844386*fedge[1])*surfVar_l; 
  edgeSurf[9] = -0.0625*(23.0*fskin[9]+7.0*fedge[9]-22.5166604983954*fskin[2]+8.660254037844386*fedge[2])*surfVar_l; 
  edgeSurf[10] = -0.0625*(23.0*fskin[10]+7.0*fedge[10]-22.5166604983954*fskin[3]+8.660254037844386*fedge[3])*surfVar_l; 
  edgeSurf[11] = 0.0625*(8.660254037844386*(fskin[15]+fedge[15])-9.0*fskin[11]+9.0*fedge[11])*surfVar_l; 
  edgeSurf[12] = -0.0625*(23.0*fskin[12]+7.0*fedge[12]-22.5166604983954*fskin[5]+8.660254037844386*fedge[5])*surfVar_l; 
  edgeSurf[13] = -0.0625*(23.0*fskin[13]+7.0*fedge[13]-22.5166604983954*fskin[6]+8.660254037844386*fedge[6])*surfVar_l; 
  edgeSurf[14] = -0.0625*(23.0*fskin[14]+7.0*fedge[14]-22.5166604983954*fskin[7]+8.660254037844386*fedge[7])*surfVar_l; 
  edgeSurf[15] = -0.0625*(23.0*fskin[15]+7.0*fedge[15]-22.5166604983954*fskin[11]+8.660254037844386*fedge[11])*surfVar_l; 
  edgeSurf[16] = 0.0625*(8.660254037844387*(fskin[19]+fedge[19])-9.0*fskin[16]+9.0*fedge[16])*surfVar_l; 
  edgeSurf[17] = 0.0625*(8.660254037844387*(fskin[21]+fedge[21])-9.0*fskin[17]+9.0*fedge[17])*surfVar_l; 
  edgeSurf[18] = 0.0625*(8.660254037844387*(fskin[22]+fedge[22])-9.0*fskin[18]+9.0*fedge[18])*surfVar_l; 
  edgeSurf[19] = -0.0125*(115.0*fskin[19]+35.0*fedge[19]-112.583302491977*fskin[16]+43.30127018922195*fedge[16])*surfVar_l; 
  edgeSurf[20] = 0.0625*(8.660254037844387*(fskin[23]+fedge[23])-9.0*fskin[20]+9.0*fedge[20])*surfVar_l; 
  edgeSurf[21] = -0.0125*(115.0*fskin[21]+35.0*fedge[21]-112.583302491977*fskin[17]+43.30127018922195*fedge[17])*surfVar_l; 
  edgeSurf[22] = -0.0125*(115.0*fskin[22]+35.0*fedge[22]-112.583302491977*fskin[18]+43.30127018922195*fedge[18])*surfVar_l; 
  edgeSurf[23] = -0.0125*(115.0*fskin[23]+35.0*fedge[23]-112.583302491977*fskin[20]+43.30127018922195*fedge[20])*surfVar_l; 

  double boundSurf[24] = {0.0}; 
  boundSurf[4] = -0.5*(3.0*fskin[4]+1.732050807568877*fskin[0])*surfVar_r; 
  boundSurf[8] = -0.5*(3.0*fskin[8]+1.732050807568877*fskin[1])*surfVar_r; 
  boundSurf[9] = -0.5*(3.0*fskin[9]+1.732050807568877*fskin[2])*surfVar_r; 
  boundSurf[10] = -0.5*(3.0*fskin[10]+1.732050807568877*fskin[3])*surfVar_r; 
  boundSurf[12] = -0.5*(3.0*fskin[12]+1.732050807568877*fskin[5])*surfVar_r; 
  boundSurf[13] = -0.5*(3.0*fskin[13]+1.732050807568877*fskin[6])*surfVar_r; 
  boundSurf[14] = -0.5*(3.0*fskin[14]+1.732050807568877*fskin[7])*surfVar_r; 
  boundSurf[15] = -0.5*(3.0*fskin[15]+1.732050807568877*fskin[11])*surfVar_r; 
  boundSurf[19] = -0.1*(15.0*fskin[19]+8.660254037844387*fskin[16])*surfVar_r; 
  boundSurf[21] = -0.1*(15.0*fskin[21]+8.660254037844387*fskin[17])*surfVar_r; 
  boundSurf[22] = -0.1*(15.0*fskin[22]+8.660254037844387*fskin[18])*surfVar_r; 
  boundSurf[23] = -0.1*(15.0*fskin[23]+8.660254037844387*fskin[20])*surfVar_r; 

  edgeSurf_incr[0] = 0.5*(facDiff[3]*edgeSurf[5]+edgeSurf[2]*facDiff[2]+edgeSurf[1]*facDiff[1]+edgeSurf[0]*facDiff[0]); 
  edgeSurf_incr[1] = 0.5*(facDiff[2]*edgeSurf[5]+edgeSurf[2]*facDiff[3]+edgeSurf[0]*facDiff[1]+facDiff[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.5*(facDiff[1]*edgeSurf[5]+edgeSurf[1]*facDiff[3]+edgeSurf[0]*facDiff[2]+facDiff[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.5*(facDiff[3]*edgeSurf[11]+facDiff[2]*edgeSurf[7]+facDiff[1]*edgeSurf[6]+facDiff[0]*edgeSurf[3]); 
  edgeSurf_incr[4] = 0.5*(facDiff[3]*edgeSurf[12]+facDiff[2]*edgeSurf[9]+facDiff[1]*edgeSurf[8]+facDiff[0]*edgeSurf[4]); 
  edgeSurf_incr[5] = 0.5*(facDiff[0]*edgeSurf[5]+edgeSurf[0]*facDiff[3]+edgeSurf[1]*facDiff[2]+facDiff[1]*edgeSurf[2]); 
  edgeSurf_incr[6] = 0.5*(facDiff[2]*edgeSurf[11]+facDiff[3]*edgeSurf[7]+facDiff[0]*edgeSurf[6]+facDiff[1]*edgeSurf[3]); 
  edgeSurf_incr[7] = 0.5*(facDiff[1]*edgeSurf[11]+facDiff[0]*edgeSurf[7]+facDiff[3]*edgeSurf[6]+facDiff[2]*edgeSurf[3]); 
  edgeSurf_incr[8] = 0.5*(facDiff[2]*edgeSurf[12]+facDiff[3]*edgeSurf[9]+facDiff[0]*edgeSurf[8]+facDiff[1]*edgeSurf[4]); 
  edgeSurf_incr[9] = 0.5*(facDiff[1]*edgeSurf[12]+facDiff[0]*edgeSurf[9]+facDiff[3]*edgeSurf[8]+facDiff[2]*edgeSurf[4]); 
  edgeSurf_incr[10] = 0.5*(facDiff[3]*edgeSurf[15]+facDiff[2]*edgeSurf[14]+facDiff[1]*edgeSurf[13]+facDiff[0]*edgeSurf[10]); 
  edgeSurf_incr[11] = 0.5*(facDiff[0]*edgeSurf[11]+facDiff[1]*edgeSurf[7]+facDiff[2]*edgeSurf[6]+edgeSurf[3]*facDiff[3]); 
  edgeSurf_incr[12] = 0.5*(facDiff[0]*edgeSurf[12]+facDiff[1]*edgeSurf[9]+facDiff[2]*edgeSurf[8]+facDiff[3]*edgeSurf[4]); 
  edgeSurf_incr[13] = 0.5*(facDiff[2]*edgeSurf[15]+facDiff[3]*edgeSurf[14]+facDiff[0]*edgeSurf[13]+facDiff[1]*edgeSurf[10]); 
  edgeSurf_incr[14] = 0.5*(facDiff[1]*edgeSurf[15]+facDiff[0]*edgeSurf[14]+facDiff[3]*edgeSurf[13]+facDiff[2]*edgeSurf[10]); 
  edgeSurf_incr[15] = 0.5*(facDiff[0]*edgeSurf[15]+facDiff[1]*edgeSurf[14]+facDiff[2]*edgeSurf[13]+facDiff[3]*edgeSurf[10]); 
  edgeSurf_incr[16] = 0.03333333333333333*(15.0*facDiff[3]*edgeSurf[20]+15.0*(facDiff[2]*edgeSurf[18]+facDiff[1]*edgeSurf[17])+15.0*facDiff[0]*edgeSurf[16]); 
  edgeSurf_incr[17] = 0.03333333333333333*(15.0*facDiff[2]*edgeSurf[20]+15.0*(facDiff[3]*edgeSurf[18]+facDiff[0]*edgeSurf[17])+15.0*facDiff[1]*edgeSurf[16]); 
  edgeSurf_incr[18] = 0.03333333333333333*(15.0*facDiff[1]*edgeSurf[20]+15.0*(facDiff[0]*edgeSurf[18]+facDiff[3]*edgeSurf[17])+15.0*facDiff[2]*edgeSurf[16]); 
  edgeSurf_incr[19] = 0.03333333333333333*(15.0*facDiff[3]*edgeSurf[23]+15.0*(facDiff[2]*edgeSurf[22]+facDiff[1]*edgeSurf[21])+15.0*facDiff[0]*edgeSurf[19]); 
  edgeSurf_incr[20] = 0.03333333333333333*(15.0*facDiff[0]*edgeSurf[20]+15.0*(facDiff[1]*edgeSurf[18]+facDiff[2]*edgeSurf[17])+15.0*facDiff[3]*edgeSurf[16]); 
  edgeSurf_incr[21] = 0.03333333333333333*(15.0*facDiff[2]*edgeSurf[23]+15.0*(facDiff[3]*edgeSurf[22]+facDiff[0]*edgeSurf[21])+15.0*facDiff[1]*edgeSurf[19]); 
  edgeSurf_incr[22] = 0.03333333333333333*(15.0*facDiff[1]*edgeSurf[23]+15.0*(facDiff[0]*edgeSurf[22]+facDiff[3]*edgeSurf[21])+15.0*facDiff[2]*edgeSurf[19]); 
  edgeSurf_incr[23] = 0.03333333333333333*(15.0*facDiff[0]*edgeSurf[23]+15.0*(facDiff[1]*edgeSurf[22]+facDiff[2]*edgeSurf[21])+15.0*facDiff[3]*edgeSurf[19]); 

  boundSurf_incr[0] = 0.5*(facDiff[3]*boundSurf[5]+boundSurf[2]*facDiff[2]+boundSurf[1]*facDiff[1]+boundSurf[0]*facDiff[0]); 
  boundSurf_incr[1] = 0.5*(facDiff[2]*boundSurf[5]+boundSurf[2]*facDiff[3]+boundSurf[0]*facDiff[1]+facDiff[0]*boundSurf[1]); 
  boundSurf_incr[2] = 0.5*(facDiff[1]*boundSurf[5]+boundSurf[1]*facDiff[3]+boundSurf[0]*facDiff[2]+facDiff[0]*boundSurf[2]); 
  boundSurf_incr[3] = 0.5*(facDiff[3]*boundSurf[11]+facDiff[2]*boundSurf[7]+facDiff[1]*boundSurf[6]+facDiff[0]*boundSurf[3]); 
  boundSurf_incr[4] = 0.5*(facDiff[3]*boundSurf[12]+facDiff[2]*boundSurf[9]+facDiff[1]*boundSurf[8]+facDiff[0]*boundSurf[4]); 
  boundSurf_incr[5] = 0.5*(facDiff[0]*boundSurf[5]+boundSurf[0]*facDiff[3]+boundSurf[1]*facDiff[2]+facDiff[1]*boundSurf[2]); 
  boundSurf_incr[6] = 0.5*(facDiff[2]*boundSurf[11]+facDiff[3]*boundSurf[7]+facDiff[0]*boundSurf[6]+facDiff[1]*boundSurf[3]); 
  boundSurf_incr[7] = 0.5*(facDiff[1]*boundSurf[11]+facDiff[0]*boundSurf[7]+facDiff[3]*boundSurf[6]+facDiff[2]*boundSurf[3]); 
  boundSurf_incr[8] = 0.5*(facDiff[2]*boundSurf[12]+facDiff[3]*boundSurf[9]+facDiff[0]*boundSurf[8]+facDiff[1]*boundSurf[4]); 
  boundSurf_incr[9] = 0.5*(facDiff[1]*boundSurf[12]+facDiff[0]*boundSurf[9]+facDiff[3]*boundSurf[8]+facDiff[2]*boundSurf[4]); 
  boundSurf_incr[10] = 0.5*(facDiff[3]*boundSurf[15]+facDiff[2]*boundSurf[14]+facDiff[1]*boundSurf[13]+facDiff[0]*boundSurf[10]); 
  boundSurf_incr[11] = 0.5*(facDiff[0]*boundSurf[11]+facDiff[1]*boundSurf[7]+facDiff[2]*boundSurf[6]+boundSurf[3]*facDiff[3]); 
  boundSurf_incr[12] = 0.5*(facDiff[0]*boundSurf[12]+facDiff[1]*boundSurf[9]+facDiff[2]*boundSurf[8]+facDiff[3]*boundSurf[4]); 
  boundSurf_incr[13] = 0.5*(facDiff[2]*boundSurf[15]+facDiff[3]*boundSurf[14]+facDiff[0]*boundSurf[13]+facDiff[1]*boundSurf[10]); 
  boundSurf_incr[14] = 0.5*(facDiff[1]*boundSurf[15]+facDiff[0]*boundSurf[14]+facDiff[3]*boundSurf[13]+facDiff[2]*boundSurf[10]); 
  boundSurf_incr[15] = 0.5*(facDiff[0]*boundSurf[15]+facDiff[1]*boundSurf[14]+facDiff[2]*boundSurf[13]+facDiff[3]*boundSurf[10]); 
  boundSurf_incr[16] = 0.03333333333333333*(15.0*facDiff[3]*boundSurf[20]+15.0*(facDiff[2]*boundSurf[18]+facDiff[1]*boundSurf[17])+15.0*facDiff[0]*boundSurf[16]); 
  boundSurf_incr[17] = 0.03333333333333333*(15.0*facDiff[2]*boundSurf[20]+15.0*(facDiff[3]*boundSurf[18]+facDiff[0]*boundSurf[17])+15.0*facDiff[1]*boundSurf[16]); 
  boundSurf_incr[18] = 0.03333333333333333*(15.0*facDiff[1]*boundSurf[20]+15.0*(facDiff[0]*boundSurf[18]+facDiff[3]*boundSurf[17])+15.0*facDiff[2]*boundSurf[16]); 
  boundSurf_incr[19] = 0.03333333333333333*(15.0*facDiff[3]*boundSurf[23]+15.0*(facDiff[2]*boundSurf[22]+facDiff[1]*boundSurf[21])+15.0*facDiff[0]*boundSurf[19]); 
  boundSurf_incr[20] = 0.03333333333333333*(15.0*facDiff[0]*boundSurf[20]+15.0*(facDiff[1]*boundSurf[18]+facDiff[2]*boundSurf[17])+15.0*facDiff[3]*boundSurf[16]); 
  boundSurf_incr[21] = 0.03333333333333333*(15.0*facDiff[2]*boundSurf[23]+15.0*(facDiff[3]*boundSurf[22]+facDiff[0]*boundSurf[21])+15.0*facDiff[1]*boundSurf[19]); 
  boundSurf_incr[22] = 0.03333333333333333*(15.0*facDiff[1]*boundSurf[23]+15.0*(facDiff[0]*boundSurf[22]+facDiff[3]*boundSurf[21])+15.0*facDiff[2]*boundSurf[19]); 
  boundSurf_incr[23] = 0.03333333333333333*(15.0*facDiff[0]*boundSurf[23]+15.0*(facDiff[1]*boundSurf[22]+facDiff[2]*boundSurf[21])+15.0*facDiff[3]*boundSurf[19]); 

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
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*rdvSq4; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*rdvSq4; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*rdvSq4; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*rdvSq4; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*rdvSq4; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*rdvSq4; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*rdvSq4; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*rdvSq4; 
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*rdvSq4; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*rdvSq4; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*rdvSq4; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*rdvSq4; 
  return 0.;

} 
