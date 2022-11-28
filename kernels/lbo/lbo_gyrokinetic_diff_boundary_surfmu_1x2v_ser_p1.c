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

  double facDiff[2]; 
  // Expand diffusion coefficient in conf basis.
  facDiff[0] = 1.414213562373095*(bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0])*m_; 
  facDiff[1] = 1.414213562373095*(bmag_inv[0]*nuVtSqSum[1]+nuVtSqSum[0]*bmag_inv[1])*m_; 

  double vol_incr[12] = {0.0};
  vol_incr[3] = 0.6123724356957944*facDiff[1]*fskin[1]*dxv[2]+0.6123724356957944*facDiff[0]*fskin[0]*dxv[2]; 
  vol_incr[5] = 0.6123724356957944*facDiff[0]*fskin[1]*dxv[2]+0.6123724356957944*fskin[0]*facDiff[1]*dxv[2]; 
  vol_incr[6] = 0.6123724356957944*facDiff[1]*dxv[2]*fskin[4]+0.6123724356957944*facDiff[0]*dxv[2]*fskin[2]; 
  vol_incr[7] = 0.6123724356957944*facDiff[0]*dxv[2]*fskin[4]+0.6123724356957944*facDiff[1]*dxv[2]*fskin[2]; 
  vol_incr[10] = 0.6123724356957944*facDiff[1]*dxv[2]*fskin[9]+0.6123724356957944*facDiff[0]*dxv[2]*fskin[8]; 
  vol_incr[11] = 0.6123724356957944*facDiff[0]*dxv[2]*fskin[9]+0.6123724356957944*facDiff[1]*dxv[2]*fskin[8]; 

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[12] = {0.0}; 
  edgeSurf[0] = (-0.5412658773652741*w[2]*fskin[3])-0.270632938682637*dxv[2]*fskin[3]-0.5412658773652741*w[2]*fedge[3]-0.270632938682637*dxv[2]*fedge[3]-0.5625*fskin[0]*w[2]+0.5625*fedge[0]*w[2]-0.28125*fskin[0]*dxv[2]+0.28125*fedge[0]*dxv[2]; 
  edgeSurf[1] = (-0.5412658773652741*w[2]*fskin[5])-0.270632938682637*dxv[2]*fskin[5]-0.5412658773652741*w[2]*fedge[5]-0.270632938682637*dxv[2]*fedge[5]-0.5625*fskin[1]*w[2]+0.5625*fedge[1]*w[2]-0.28125*fskin[1]*dxv[2]+0.28125*fedge[1]*dxv[2]; 
  edgeSurf[2] = (-0.5412658773652741*w[2]*fskin[6])-0.270632938682637*dxv[2]*fskin[6]-0.5412658773652741*w[2]*fedge[6]-0.270632938682637*dxv[2]*fedge[6]-0.5625*fskin[2]*w[2]+0.5625*fedge[2]*w[2]-0.28125*dxv[2]*fskin[2]+0.28125*dxv[2]*fedge[2]; 
  edgeSurf[3] = (-1.4375*w[2]*fskin[3])-0.71875*dxv[2]*fskin[3]-0.4375*w[2]*fedge[3]-0.21875*dxv[2]*fedge[3]-1.407291281149713*fskin[0]*w[2]+0.5412658773652741*fedge[0]*w[2]-0.7036456405748563*fskin[0]*dxv[2]+0.270632938682637*fedge[0]*dxv[2]; 
  edgeSurf[4] = (-0.5412658773652741*w[2]*fskin[7])-0.270632938682637*dxv[2]*fskin[7]-0.5412658773652741*w[2]*fedge[7]-0.270632938682637*dxv[2]*fedge[7]-0.5625*w[2]*fskin[4]-0.28125*dxv[2]*fskin[4]+0.5625*w[2]*fedge[4]+0.28125*dxv[2]*fedge[4]; 
  edgeSurf[5] = (-1.4375*w[2]*fskin[5])-0.71875*dxv[2]*fskin[5]-0.4375*w[2]*fedge[5]-0.21875*dxv[2]*fedge[5]-1.407291281149713*fskin[1]*w[2]+0.5412658773652741*fedge[1]*w[2]-0.7036456405748563*fskin[1]*dxv[2]+0.270632938682637*fedge[1]*dxv[2]; 
  edgeSurf[6] = (-1.4375*w[2]*fskin[6])-0.71875*dxv[2]*fskin[6]-0.4375*w[2]*fedge[6]-0.21875*dxv[2]*fedge[6]-1.407291281149713*fskin[2]*w[2]+0.5412658773652741*fedge[2]*w[2]-0.7036456405748563*dxv[2]*fskin[2]+0.270632938682637*dxv[2]*fedge[2]; 
  edgeSurf[7] = (-1.4375*w[2]*fskin[7])-0.71875*dxv[2]*fskin[7]-0.4375*w[2]*fedge[7]-0.21875*dxv[2]*fedge[7]-1.407291281149713*w[2]*fskin[4]-0.7036456405748563*dxv[2]*fskin[4]+0.5412658773652741*w[2]*fedge[4]+0.270632938682637*dxv[2]*fedge[4]; 
  edgeSurf[8] = (-0.5412658773652742*w[2]*fskin[10])-0.2706329386826371*dxv[2]*fskin[10]-0.5412658773652742*w[2]*fedge[10]-0.2706329386826371*dxv[2]*fedge[10]-0.5625*w[2]*fskin[8]-0.28125*dxv[2]*fskin[8]+0.5625*w[2]*fedge[8]+0.28125*dxv[2]*fedge[8]; 
  edgeSurf[9] = (-0.5412658773652742*w[2]*fskin[11])-0.2706329386826371*dxv[2]*fskin[11]-0.5412658773652742*w[2]*fedge[11]-0.2706329386826371*dxv[2]*fedge[11]-0.5625*w[2]*fskin[9]-0.28125*dxv[2]*fskin[9]+0.5625*w[2]*fedge[9]+0.28125*dxv[2]*fedge[9]; 
  edgeSurf[10] = (-1.4375*w[2]*fskin[10])-0.71875*dxv[2]*fskin[10]-0.4375*w[2]*fedge[10]-0.21875*dxv[2]*fedge[10]-1.407291281149713*w[2]*fskin[8]-0.7036456405748563*dxv[2]*fskin[8]+0.5412658773652742*w[2]*fedge[8]+0.2706329386826371*dxv[2]*fedge[8]; 
  edgeSurf[11] = (-1.4375*w[2]*fskin[11])-0.71875*dxv[2]*fskin[11]-0.4375*w[2]*fedge[11]-0.21875*dxv[2]*fedge[11]-1.407291281149713*w[2]*fskin[9]-0.7036456405748563*dxv[2]*fskin[9]+0.5412658773652742*w[2]*fedge[9]+0.2706329386826371*dxv[2]*fedge[9]; 

  double boundSurf[12] = {0.0}; 
  boundSurf[3] = (-1.5*w[2]*fskin[3])+0.75*dxv[2]*fskin[3]+0.8660254037844386*fskin[0]*w[2]-0.4330127018922193*fskin[0]*dxv[2]; 
  boundSurf[5] = (-1.5*w[2]*fskin[5])+0.75*dxv[2]*fskin[5]+0.8660254037844386*fskin[1]*w[2]-0.4330127018922193*fskin[1]*dxv[2]; 
  boundSurf[6] = (-1.5*w[2]*fskin[6])+0.75*dxv[2]*fskin[6]+0.8660254037844386*fskin[2]*w[2]-0.4330127018922193*dxv[2]*fskin[2]; 
  boundSurf[7] = (-1.5*w[2]*fskin[7])+0.75*dxv[2]*fskin[7]+0.8660254037844386*w[2]*fskin[4]-0.4330127018922193*dxv[2]*fskin[4]; 
  boundSurf[10] = (-1.5*w[2]*fskin[10])+0.75*dxv[2]*fskin[10]+0.8660254037844387*w[2]*fskin[8]-0.4330127018922194*dxv[2]*fskin[8]; 
  boundSurf[11] = (-1.5*w[2]*fskin[11])+0.75*dxv[2]*fskin[11]+0.8660254037844387*w[2]*fskin[9]-0.4330127018922194*dxv[2]*fskin[9]; 

  edgeSurf_incr[0] = 0.7071067811865475*edgeSurf[1]*facDiff[1]+0.7071067811865475*edgeSurf[0]*facDiff[0]; 
  edgeSurf_incr[1] = 0.7071067811865475*edgeSurf[0]*facDiff[1]+0.7071067811865475*facDiff[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865475*facDiff[1]*edgeSurf[4]+0.7071067811865475*facDiff[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.7071067811865475*facDiff[1]*edgeSurf[5]+0.7071067811865475*facDiff[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.7071067811865475*facDiff[0]*edgeSurf[4]+0.7071067811865475*facDiff[1]*edgeSurf[2]; 
  edgeSurf_incr[5] = 0.7071067811865475*facDiff[0]*edgeSurf[5]+0.7071067811865475*facDiff[1]*edgeSurf[3]; 
  edgeSurf_incr[6] = 0.7071067811865475*facDiff[1]*edgeSurf[7]+0.7071067811865475*facDiff[0]*edgeSurf[6]; 
  edgeSurf_incr[7] = 0.7071067811865475*facDiff[0]*edgeSurf[7]+0.7071067811865475*facDiff[1]*edgeSurf[6]; 
  edgeSurf_incr[8] = 0.7071067811865475*facDiff[1]*edgeSurf[9]+0.7071067811865475*facDiff[0]*edgeSurf[8]; 
  edgeSurf_incr[9] = 0.7071067811865475*facDiff[0]*edgeSurf[9]+0.7071067811865475*facDiff[1]*edgeSurf[8]; 
  edgeSurf_incr[10] = 0.7071067811865475*facDiff[1]*edgeSurf[11]+0.7071067811865475*facDiff[0]*edgeSurf[10]; 
  edgeSurf_incr[11] = 0.7071067811865475*facDiff[0]*edgeSurf[11]+0.7071067811865475*facDiff[1]*edgeSurf[10]; 

  boundSurf_incr[0] = 0.7071067811865475*boundSurf[1]*facDiff[1]+0.7071067811865475*boundSurf[0]*facDiff[0]; 
  boundSurf_incr[1] = 0.7071067811865475*boundSurf[0]*facDiff[1]+0.7071067811865475*facDiff[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865475*facDiff[1]*boundSurf[4]+0.7071067811865475*facDiff[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.7071067811865475*facDiff[1]*boundSurf[5]+0.7071067811865475*facDiff[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.7071067811865475*facDiff[0]*boundSurf[4]+0.7071067811865475*facDiff[1]*boundSurf[2]; 
  boundSurf_incr[5] = 0.7071067811865475*facDiff[0]*boundSurf[5]+0.7071067811865475*facDiff[1]*boundSurf[3]; 
  boundSurf_incr[6] = 0.7071067811865475*facDiff[1]*boundSurf[7]+0.7071067811865475*facDiff[0]*boundSurf[6]; 
  boundSurf_incr[7] = 0.7071067811865475*facDiff[0]*boundSurf[7]+0.7071067811865475*facDiff[1]*boundSurf[6]; 
  boundSurf_incr[8] = 0.7071067811865475*facDiff[1]*boundSurf[9]+0.7071067811865475*facDiff[0]*boundSurf[8]; 
  boundSurf_incr[9] = 0.7071067811865475*facDiff[0]*boundSurf[9]+0.7071067811865475*facDiff[1]*boundSurf[8]; 
  boundSurf_incr[10] = 0.7071067811865475*facDiff[1]*boundSurf[11]+0.7071067811865475*facDiff[0]*boundSurf[10]; 
  boundSurf_incr[11] = 0.7071067811865475*facDiff[0]*boundSurf[11]+0.7071067811865475*facDiff[1]*boundSurf[10]; 


  } else { 

  double edgeSurf[12] = {0.0}; 
  edgeSurf[0] = 0.5412658773652741*w[2]*fskin[3]-0.270632938682637*dxv[2]*fskin[3]+0.5412658773652741*w[2]*fedge[3]-0.270632938682637*dxv[2]*fedge[3]-0.5625*fskin[0]*w[2]+0.5625*fedge[0]*w[2]+0.28125*fskin[0]*dxv[2]-0.28125*fedge[0]*dxv[2]; 
  edgeSurf[1] = 0.5412658773652741*w[2]*fskin[5]-0.270632938682637*dxv[2]*fskin[5]+0.5412658773652741*w[2]*fedge[5]-0.270632938682637*dxv[2]*fedge[5]-0.5625*fskin[1]*w[2]+0.5625*fedge[1]*w[2]+0.28125*fskin[1]*dxv[2]-0.28125*fedge[1]*dxv[2]; 
  edgeSurf[2] = 0.5412658773652741*w[2]*fskin[6]-0.270632938682637*dxv[2]*fskin[6]+0.5412658773652741*w[2]*fedge[6]-0.270632938682637*dxv[2]*fedge[6]-0.5625*fskin[2]*w[2]+0.5625*fedge[2]*w[2]+0.28125*dxv[2]*fskin[2]-0.28125*dxv[2]*fedge[2]; 
  edgeSurf[3] = (-1.4375*w[2]*fskin[3])+0.71875*dxv[2]*fskin[3]-0.4375*w[2]*fedge[3]+0.21875*dxv[2]*fedge[3]+1.407291281149713*fskin[0]*w[2]-0.5412658773652741*fedge[0]*w[2]-0.7036456405748563*fskin[0]*dxv[2]+0.270632938682637*fedge[0]*dxv[2]; 
  edgeSurf[4] = 0.5412658773652741*w[2]*fskin[7]-0.270632938682637*dxv[2]*fskin[7]+0.5412658773652741*w[2]*fedge[7]-0.270632938682637*dxv[2]*fedge[7]-0.5625*w[2]*fskin[4]+0.28125*dxv[2]*fskin[4]+0.5625*w[2]*fedge[4]-0.28125*dxv[2]*fedge[4]; 
  edgeSurf[5] = (-1.4375*w[2]*fskin[5])+0.71875*dxv[2]*fskin[5]-0.4375*w[2]*fedge[5]+0.21875*dxv[2]*fedge[5]+1.407291281149713*fskin[1]*w[2]-0.5412658773652741*fedge[1]*w[2]-0.7036456405748563*fskin[1]*dxv[2]+0.270632938682637*fedge[1]*dxv[2]; 
  edgeSurf[6] = (-1.4375*w[2]*fskin[6])+0.71875*dxv[2]*fskin[6]-0.4375*w[2]*fedge[6]+0.21875*dxv[2]*fedge[6]+1.407291281149713*fskin[2]*w[2]-0.5412658773652741*fedge[2]*w[2]-0.7036456405748563*dxv[2]*fskin[2]+0.270632938682637*dxv[2]*fedge[2]; 
  edgeSurf[7] = (-1.4375*w[2]*fskin[7])+0.71875*dxv[2]*fskin[7]-0.4375*w[2]*fedge[7]+0.21875*dxv[2]*fedge[7]+1.407291281149713*w[2]*fskin[4]-0.7036456405748563*dxv[2]*fskin[4]-0.5412658773652741*w[2]*fedge[4]+0.270632938682637*dxv[2]*fedge[4]; 
  edgeSurf[8] = 0.5412658773652742*w[2]*fskin[10]-0.2706329386826371*dxv[2]*fskin[10]+0.5412658773652742*w[2]*fedge[10]-0.2706329386826371*dxv[2]*fedge[10]-0.5625*w[2]*fskin[8]+0.28125*dxv[2]*fskin[8]+0.5625*w[2]*fedge[8]-0.28125*dxv[2]*fedge[8]; 
  edgeSurf[9] = 0.5412658773652742*w[2]*fskin[11]-0.2706329386826371*dxv[2]*fskin[11]+0.5412658773652742*w[2]*fedge[11]-0.2706329386826371*dxv[2]*fedge[11]-0.5625*w[2]*fskin[9]+0.28125*dxv[2]*fskin[9]+0.5625*w[2]*fedge[9]-0.28125*dxv[2]*fedge[9]; 
  edgeSurf[10] = (-1.4375*w[2]*fskin[10])+0.71875*dxv[2]*fskin[10]-0.4375*w[2]*fedge[10]+0.21875*dxv[2]*fedge[10]+1.407291281149713*w[2]*fskin[8]-0.7036456405748563*dxv[2]*fskin[8]-0.5412658773652742*w[2]*fedge[8]+0.2706329386826371*dxv[2]*fedge[8]; 
  edgeSurf[11] = (-1.4375*w[2]*fskin[11])+0.71875*dxv[2]*fskin[11]-0.4375*w[2]*fedge[11]+0.21875*dxv[2]*fedge[11]+1.407291281149713*w[2]*fskin[9]-0.7036456405748563*dxv[2]*fskin[9]-0.5412658773652742*w[2]*fedge[9]+0.2706329386826371*dxv[2]*fedge[9]; 

  double boundSurf[12] = {0.0}; 
  boundSurf[3] = (-1.5*w[2]*fskin[3])-0.75*dxv[2]*fskin[3]-0.8660254037844386*fskin[0]*w[2]-0.4330127018922193*fskin[0]*dxv[2]; 
  boundSurf[5] = (-1.5*w[2]*fskin[5])-0.75*dxv[2]*fskin[5]-0.8660254037844386*fskin[1]*w[2]-0.4330127018922193*fskin[1]*dxv[2]; 
  boundSurf[6] = (-1.5*w[2]*fskin[6])-0.75*dxv[2]*fskin[6]-0.8660254037844386*fskin[2]*w[2]-0.4330127018922193*dxv[2]*fskin[2]; 
  boundSurf[7] = (-1.5*w[2]*fskin[7])-0.75*dxv[2]*fskin[7]-0.8660254037844386*w[2]*fskin[4]-0.4330127018922193*dxv[2]*fskin[4]; 
  boundSurf[10] = (-1.5*w[2]*fskin[10])-0.75*dxv[2]*fskin[10]-0.8660254037844387*w[2]*fskin[8]-0.4330127018922194*dxv[2]*fskin[8]; 
  boundSurf[11] = (-1.5*w[2]*fskin[11])-0.75*dxv[2]*fskin[11]-0.8660254037844387*w[2]*fskin[9]-0.4330127018922194*dxv[2]*fskin[9]; 

  edgeSurf_incr[0] = 0.7071067811865475*edgeSurf[1]*facDiff[1]+0.7071067811865475*edgeSurf[0]*facDiff[0]; 
  edgeSurf_incr[1] = 0.7071067811865475*edgeSurf[0]*facDiff[1]+0.7071067811865475*facDiff[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865475*facDiff[1]*edgeSurf[4]+0.7071067811865475*facDiff[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.7071067811865475*facDiff[1]*edgeSurf[5]+0.7071067811865475*facDiff[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.7071067811865475*facDiff[0]*edgeSurf[4]+0.7071067811865475*facDiff[1]*edgeSurf[2]; 
  edgeSurf_incr[5] = 0.7071067811865475*facDiff[0]*edgeSurf[5]+0.7071067811865475*facDiff[1]*edgeSurf[3]; 
  edgeSurf_incr[6] = 0.7071067811865475*facDiff[1]*edgeSurf[7]+0.7071067811865475*facDiff[0]*edgeSurf[6]; 
  edgeSurf_incr[7] = 0.7071067811865475*facDiff[0]*edgeSurf[7]+0.7071067811865475*facDiff[1]*edgeSurf[6]; 
  edgeSurf_incr[8] = 0.7071067811865475*facDiff[1]*edgeSurf[9]+0.7071067811865475*facDiff[0]*edgeSurf[8]; 
  edgeSurf_incr[9] = 0.7071067811865475*facDiff[0]*edgeSurf[9]+0.7071067811865475*facDiff[1]*edgeSurf[8]; 
  edgeSurf_incr[10] = 0.7071067811865475*facDiff[1]*edgeSurf[11]+0.7071067811865475*facDiff[0]*edgeSurf[10]; 
  edgeSurf_incr[11] = 0.7071067811865475*facDiff[0]*edgeSurf[11]+0.7071067811865475*facDiff[1]*edgeSurf[10]; 

  boundSurf_incr[0] = 0.7071067811865475*boundSurf[1]*facDiff[1]+0.7071067811865475*boundSurf[0]*facDiff[0]; 
  boundSurf_incr[1] = 0.7071067811865475*boundSurf[0]*facDiff[1]+0.7071067811865475*facDiff[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865475*facDiff[1]*boundSurf[4]+0.7071067811865475*facDiff[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.7071067811865475*facDiff[1]*boundSurf[5]+0.7071067811865475*facDiff[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.7071067811865475*facDiff[0]*boundSurf[4]+0.7071067811865475*facDiff[1]*boundSurf[2]; 
  boundSurf_incr[5] = 0.7071067811865475*facDiff[0]*boundSurf[5]+0.7071067811865475*facDiff[1]*boundSurf[3]; 
  boundSurf_incr[6] = 0.7071067811865475*facDiff[1]*boundSurf[7]+0.7071067811865475*facDiff[0]*boundSurf[6]; 
  boundSurf_incr[7] = 0.7071067811865475*facDiff[0]*boundSurf[7]+0.7071067811865475*facDiff[1]*boundSurf[6]; 
  boundSurf_incr[8] = 0.7071067811865475*facDiff[1]*boundSurf[9]+0.7071067811865475*facDiff[0]*boundSurf[8]; 
  boundSurf_incr[9] = 0.7071067811865475*facDiff[0]*boundSurf[9]+0.7071067811865475*facDiff[1]*boundSurf[8]; 
  boundSurf_incr[10] = 0.7071067811865475*facDiff[1]*boundSurf[11]+0.7071067811865475*facDiff[0]*boundSurf[10]; 
  boundSurf_incr[11] = 0.7071067811865475*facDiff[0]*boundSurf[11]+0.7071067811865475*facDiff[1]*boundSurf[10]; 

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
