#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_boundary_surfmu_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
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

  double facDiff[4]; 
  // Expand diffusion coefficient in conf basis.
  facDiff[0] = (bmag_inv[3]*nuVtSqSum[3]+bmag_inv[2]*nuVtSqSum[2]+bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0])*m_; 
  facDiff[1] = (bmag_inv[2]*nuVtSqSum[3]+nuVtSqSum[2]*bmag_inv[3]+bmag_inv[0]*nuVtSqSum[1]+nuVtSqSum[0]*bmag_inv[1])*m_; 
  facDiff[2] = (bmag_inv[1]*nuVtSqSum[3]+nuVtSqSum[1]*bmag_inv[3]+bmag_inv[0]*nuVtSqSum[2]+nuVtSqSum[0]*bmag_inv[2])*m_; 
  facDiff[3] = (bmag_inv[0]*nuVtSqSum[3]+nuVtSqSum[0]*bmag_inv[3]+bmag_inv[1]*nuVtSqSum[2]+nuVtSqSum[1]*bmag_inv[2])*m_; 

  double vol_incr[24] = {0.0};
  vol_incr[4] = 0.4330127018922193*dxv[3]*facDiff[3]*fskin[5]+0.4330127018922193*facDiff[2]*fskin[2]*dxv[3]+0.4330127018922193*facDiff[1]*fskin[1]*dxv[3]+0.4330127018922193*facDiff[0]*fskin[0]*dxv[3]; 
  vol_incr[8] = 0.4330127018922193*facDiff[2]*dxv[3]*fskin[5]+0.4330127018922193*fskin[2]*dxv[3]*facDiff[3]+0.4330127018922193*facDiff[0]*fskin[1]*dxv[3]+0.4330127018922193*fskin[0]*facDiff[1]*dxv[3]; 
  vol_incr[9] = 0.4330127018922193*facDiff[1]*dxv[3]*fskin[5]+0.4330127018922193*fskin[1]*dxv[3]*facDiff[3]+0.4330127018922193*facDiff[0]*fskin[2]*dxv[3]+0.4330127018922193*fskin[0]*facDiff[2]*dxv[3]; 
  vol_incr[10] = 0.4330127018922193*dxv[3]*facDiff[3]*fskin[11]+0.4330127018922193*facDiff[2]*dxv[3]*fskin[7]+0.4330127018922193*facDiff[1]*dxv[3]*fskin[6]+0.4330127018922193*facDiff[0]*dxv[3]*fskin[3]; 
  vol_incr[12] = 0.4330127018922193*facDiff[0]*dxv[3]*fskin[5]+0.4330127018922193*fskin[0]*dxv[3]*facDiff[3]+0.4330127018922193*facDiff[1]*fskin[2]*dxv[3]+0.4330127018922193*fskin[1]*facDiff[2]*dxv[3]; 
  vol_incr[13] = 0.4330127018922193*facDiff[2]*dxv[3]*fskin[11]+0.4330127018922193*dxv[3]*facDiff[3]*fskin[7]+0.4330127018922193*facDiff[0]*dxv[3]*fskin[6]+0.4330127018922193*facDiff[1]*dxv[3]*fskin[3]; 
  vol_incr[14] = 0.4330127018922193*facDiff[1]*dxv[3]*fskin[11]+0.4330127018922193*facDiff[0]*dxv[3]*fskin[7]+0.4330127018922193*dxv[3]*facDiff[3]*fskin[6]+0.4330127018922193*facDiff[2]*dxv[3]*fskin[3]; 
  vol_incr[15] = 0.4330127018922193*facDiff[0]*dxv[3]*fskin[11]+0.4330127018922193*facDiff[1]*dxv[3]*fskin[7]+0.4330127018922193*facDiff[2]*dxv[3]*fskin[6]+0.4330127018922193*dxv[3]*facDiff[3]*fskin[3]; 
  vol_incr[19] = 0.4330127018922194*dxv[3]*facDiff[3]*fskin[20]+0.4330127018922193*facDiff[2]*dxv[3]*fskin[18]+0.4330127018922193*facDiff[1]*dxv[3]*fskin[17]+0.4330127018922194*facDiff[0]*dxv[3]*fskin[16]; 
  vol_incr[21] = 0.4330127018922193*facDiff[2]*dxv[3]*fskin[20]+0.4330127018922194*dxv[3]*facDiff[3]*fskin[18]+0.4330127018922194*facDiff[0]*dxv[3]*fskin[17]+0.4330127018922193*facDiff[1]*dxv[3]*fskin[16]; 
  vol_incr[22] = 0.4330127018922193*facDiff[1]*dxv[3]*fskin[20]+0.4330127018922194*facDiff[0]*dxv[3]*fskin[18]+0.4330127018922194*dxv[3]*facDiff[3]*fskin[17]+0.4330127018922193*facDiff[2]*dxv[3]*fskin[16]; 
  vol_incr[23] = 0.4330127018922194*facDiff[0]*dxv[3]*fskin[20]+0.4330127018922193*facDiff[1]*dxv[3]*fskin[18]+0.4330127018922193*facDiff[2]*dxv[3]*fskin[17]+0.4330127018922194*dxv[3]*facDiff[3]*fskin[16]; 

  double edgeSurf_incr[24] = {0.0}; 
  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[24] = {0.0}; 
  edgeSurf[0] = (-0.5412658773652741*w[3]*fskin[4])-0.270632938682637*dxv[3]*fskin[4]-0.5412658773652741*w[3]*fedge[4]-0.270632938682637*dxv[3]*fedge[4]-0.5625*fskin[0]*w[3]+0.5625*fedge[0]*w[3]-0.28125*fskin[0]*dxv[3]+0.28125*fedge[0]*dxv[3]; 
  edgeSurf[1] = (-0.5412658773652741*w[3]*fskin[8])-0.270632938682637*dxv[3]*fskin[8]-0.5412658773652741*w[3]*fedge[8]-0.270632938682637*dxv[3]*fedge[8]-0.5625*fskin[1]*w[3]+0.5625*fedge[1]*w[3]-0.28125*fskin[1]*dxv[3]+0.28125*fedge[1]*dxv[3]; 
  edgeSurf[2] = (-0.5412658773652741*w[3]*fskin[9])-0.270632938682637*dxv[3]*fskin[9]-0.5412658773652741*w[3]*fedge[9]-0.270632938682637*dxv[3]*fedge[9]-0.5625*fskin[2]*w[3]+0.5625*fedge[2]*w[3]-0.28125*fskin[2]*dxv[3]+0.28125*fedge[2]*dxv[3]; 
  edgeSurf[3] = (-0.5412658773652741*w[3]*fskin[10])-0.270632938682637*dxv[3]*fskin[10]-0.5412658773652741*w[3]*fedge[10]-0.270632938682637*dxv[3]*fedge[10]-0.5625*fskin[3]*w[3]+0.5625*fedge[3]*w[3]-0.28125*dxv[3]*fskin[3]+0.28125*dxv[3]*fedge[3]; 
  edgeSurf[4] = (-1.4375*w[3]*fskin[4])-0.71875*dxv[3]*fskin[4]-0.4375*w[3]*fedge[4]-0.21875*dxv[3]*fedge[4]-1.407291281149713*fskin[0]*w[3]+0.5412658773652741*fedge[0]*w[3]-0.7036456405748563*fskin[0]*dxv[3]+0.270632938682637*fedge[0]*dxv[3]; 
  edgeSurf[5] = (-0.5412658773652741*w[3]*fskin[12])-0.270632938682637*dxv[3]*fskin[12]-0.5412658773652741*w[3]*fedge[12]-0.270632938682637*dxv[3]*fedge[12]-0.5625*w[3]*fskin[5]-0.28125*dxv[3]*fskin[5]+0.5625*w[3]*fedge[5]+0.28125*dxv[3]*fedge[5]; 
  edgeSurf[6] = (-0.5412658773652741*w[3]*fskin[13])-0.270632938682637*dxv[3]*fskin[13]-0.5412658773652741*w[3]*fedge[13]-0.270632938682637*dxv[3]*fedge[13]-0.5625*w[3]*fskin[6]-0.28125*dxv[3]*fskin[6]+0.5625*w[3]*fedge[6]+0.28125*dxv[3]*fedge[6]; 
  edgeSurf[7] = (-0.5412658773652741*w[3]*fskin[14])-0.270632938682637*dxv[3]*fskin[14]-0.5412658773652741*w[3]*fedge[14]-0.270632938682637*dxv[3]*fedge[14]-0.5625*w[3]*fskin[7]-0.28125*dxv[3]*fskin[7]+0.5625*w[3]*fedge[7]+0.28125*dxv[3]*fedge[7]; 
  edgeSurf[8] = (-1.4375*w[3]*fskin[8])-0.71875*dxv[3]*fskin[8]-0.4375*w[3]*fedge[8]-0.21875*dxv[3]*fedge[8]-1.407291281149713*fskin[1]*w[3]+0.5412658773652741*fedge[1]*w[3]-0.7036456405748563*fskin[1]*dxv[3]+0.270632938682637*fedge[1]*dxv[3]; 
  edgeSurf[9] = (-1.4375*w[3]*fskin[9])-0.71875*dxv[3]*fskin[9]-0.4375*w[3]*fedge[9]-0.21875*dxv[3]*fedge[9]-1.407291281149713*fskin[2]*w[3]+0.5412658773652741*fedge[2]*w[3]-0.7036456405748563*fskin[2]*dxv[3]+0.270632938682637*fedge[2]*dxv[3]; 
  edgeSurf[10] = (-1.4375*w[3]*fskin[10])-0.71875*dxv[3]*fskin[10]-0.4375*w[3]*fedge[10]-0.21875*dxv[3]*fedge[10]-1.407291281149713*fskin[3]*w[3]+0.5412658773652741*fedge[3]*w[3]-0.7036456405748563*dxv[3]*fskin[3]+0.270632938682637*dxv[3]*fedge[3]; 
  edgeSurf[11] = (-0.5412658773652741*w[3]*fskin[15])-0.270632938682637*dxv[3]*fskin[15]-0.5412658773652741*w[3]*fedge[15]-0.270632938682637*dxv[3]*fedge[15]-0.5625*w[3]*fskin[11]-0.28125*dxv[3]*fskin[11]+0.5625*w[3]*fedge[11]+0.28125*dxv[3]*fedge[11]; 
  edgeSurf[12] = (-1.4375*w[3]*fskin[12])-0.71875*dxv[3]*fskin[12]-0.4375*w[3]*fedge[12]-0.21875*dxv[3]*fedge[12]-1.407291281149713*w[3]*fskin[5]-0.7036456405748563*dxv[3]*fskin[5]+0.5412658773652741*w[3]*fedge[5]+0.270632938682637*dxv[3]*fedge[5]; 
  edgeSurf[13] = (-1.4375*w[3]*fskin[13])-0.71875*dxv[3]*fskin[13]-0.4375*w[3]*fedge[13]-0.21875*dxv[3]*fedge[13]-1.407291281149713*w[3]*fskin[6]-0.7036456405748563*dxv[3]*fskin[6]+0.5412658773652741*w[3]*fedge[6]+0.270632938682637*dxv[3]*fedge[6]; 
  edgeSurf[14] = (-1.4375*w[3]*fskin[14])-0.71875*dxv[3]*fskin[14]-0.4375*w[3]*fedge[14]-0.21875*dxv[3]*fedge[14]-1.407291281149713*w[3]*fskin[7]-0.7036456405748563*dxv[3]*fskin[7]+0.5412658773652741*w[3]*fedge[7]+0.270632938682637*dxv[3]*fedge[7]; 
  edgeSurf[15] = (-1.4375*w[3]*fskin[15])-0.71875*dxv[3]*fskin[15]-0.4375*w[3]*fedge[15]-0.21875*dxv[3]*fedge[15]-1.407291281149713*w[3]*fskin[11]-0.7036456405748563*dxv[3]*fskin[11]+0.5412658773652741*w[3]*fedge[11]+0.270632938682637*dxv[3]*fedge[11]; 
  edgeSurf[16] = (-0.5412658773652742*w[3]*fskin[19])-0.2706329386826371*dxv[3]*fskin[19]-0.5412658773652742*w[3]*fedge[19]-0.2706329386826371*dxv[3]*fedge[19]-0.5625*w[3]*fskin[16]-0.28125*dxv[3]*fskin[16]+0.5625*w[3]*fedge[16]+0.28125*dxv[3]*fedge[16]; 
  edgeSurf[17] = (-0.5412658773652742*w[3]*fskin[21])-0.2706329386826371*dxv[3]*fskin[21]-0.5412658773652742*w[3]*fedge[21]-0.2706329386826371*dxv[3]*fedge[21]-0.5625*w[3]*fskin[17]-0.28125*dxv[3]*fskin[17]+0.5625*w[3]*fedge[17]+0.28125*dxv[3]*fedge[17]; 
  edgeSurf[18] = (-0.5412658773652742*w[3]*fskin[22])-0.2706329386826371*dxv[3]*fskin[22]-0.5412658773652742*w[3]*fedge[22]-0.2706329386826371*dxv[3]*fedge[22]-0.5625*w[3]*fskin[18]-0.28125*dxv[3]*fskin[18]+0.5625*w[3]*fedge[18]+0.28125*dxv[3]*fedge[18]; 
  edgeSurf[19] = (-1.4375*w[3]*fskin[19])-0.71875*dxv[3]*fskin[19]-0.4375*w[3]*fedge[19]-0.21875*dxv[3]*fedge[19]-1.407291281149713*w[3]*fskin[16]-0.7036456405748563*dxv[3]*fskin[16]+0.5412658773652742*w[3]*fedge[16]+0.2706329386826371*dxv[3]*fedge[16]; 
  edgeSurf[20] = (-0.5412658773652742*w[3]*fskin[23])-0.2706329386826371*dxv[3]*fskin[23]-0.5412658773652742*w[3]*fedge[23]-0.2706329386826371*dxv[3]*fedge[23]-0.5625*w[3]*fskin[20]-0.28125*dxv[3]*fskin[20]+0.5625*w[3]*fedge[20]+0.28125*dxv[3]*fedge[20]; 
  edgeSurf[21] = (-1.4375*w[3]*fskin[21])-0.71875*dxv[3]*fskin[21]-0.4375*w[3]*fedge[21]-0.21875*dxv[3]*fedge[21]-1.407291281149713*w[3]*fskin[17]-0.7036456405748563*dxv[3]*fskin[17]+0.5412658773652742*w[3]*fedge[17]+0.2706329386826371*dxv[3]*fedge[17]; 
  edgeSurf[22] = (-1.4375*w[3]*fskin[22])-0.71875*dxv[3]*fskin[22]-0.4375*w[3]*fedge[22]-0.21875*dxv[3]*fedge[22]-1.407291281149713*w[3]*fskin[18]-0.7036456405748563*dxv[3]*fskin[18]+0.5412658773652742*w[3]*fedge[18]+0.2706329386826371*dxv[3]*fedge[18]; 
  edgeSurf[23] = (-1.4375*w[3]*fskin[23])-0.71875*dxv[3]*fskin[23]-0.4375*w[3]*fedge[23]-0.21875*dxv[3]*fedge[23]-1.407291281149713*w[3]*fskin[20]-0.7036456405748563*dxv[3]*fskin[20]+0.5412658773652742*w[3]*fedge[20]+0.2706329386826371*dxv[3]*fedge[20]; 

  double boundSurf[24] = {0.0}; 
  boundSurf[4] = (-1.5*w[3]*fskin[4])+0.75*dxv[3]*fskin[4]+0.8660254037844386*fskin[0]*w[3]-0.4330127018922193*fskin[0]*dxv[3]; 
  boundSurf[8] = (-1.5*w[3]*fskin[8])+0.75*dxv[3]*fskin[8]+0.8660254037844386*fskin[1]*w[3]-0.4330127018922193*fskin[1]*dxv[3]; 
  boundSurf[9] = (-1.5*w[3]*fskin[9])+0.75*dxv[3]*fskin[9]+0.8660254037844386*fskin[2]*w[3]-0.4330127018922193*fskin[2]*dxv[3]; 
  boundSurf[10] = (-1.5*w[3]*fskin[10])+0.75*dxv[3]*fskin[10]+0.8660254037844386*fskin[3]*w[3]-0.4330127018922193*dxv[3]*fskin[3]; 
  boundSurf[12] = (-1.5*w[3]*fskin[12])+0.75*dxv[3]*fskin[12]+0.8660254037844386*w[3]*fskin[5]-0.4330127018922193*dxv[3]*fskin[5]; 
  boundSurf[13] = (-1.5*w[3]*fskin[13])+0.75*dxv[3]*fskin[13]+0.8660254037844386*w[3]*fskin[6]-0.4330127018922193*dxv[3]*fskin[6]; 
  boundSurf[14] = (-1.5*w[3]*fskin[14])+0.75*dxv[3]*fskin[14]+0.8660254037844386*w[3]*fskin[7]-0.4330127018922193*dxv[3]*fskin[7]; 
  boundSurf[15] = (-1.5*w[3]*fskin[15])+0.75*dxv[3]*fskin[15]+0.8660254037844386*w[3]*fskin[11]-0.4330127018922193*dxv[3]*fskin[11]; 
  boundSurf[19] = (-1.5*w[3]*fskin[19])+0.75*dxv[3]*fskin[19]+0.8660254037844387*w[3]*fskin[16]-0.4330127018922194*dxv[3]*fskin[16]; 
  boundSurf[21] = (-1.5*w[3]*fskin[21])+0.75*dxv[3]*fskin[21]+0.8660254037844387*w[3]*fskin[17]-0.4330127018922194*dxv[3]*fskin[17]; 
  boundSurf[22] = (-1.5*w[3]*fskin[22])+0.75*dxv[3]*fskin[22]+0.8660254037844387*w[3]*fskin[18]-0.4330127018922194*dxv[3]*fskin[18]; 
  boundSurf[23] = (-1.5*w[3]*fskin[23])+0.75*dxv[3]*fskin[23]+0.8660254037844387*w[3]*fskin[20]-0.4330127018922194*dxv[3]*fskin[20]; 

  edgeSurf_incr[0] = 0.5*facDiff[3]*edgeSurf[5]+0.5*edgeSurf[2]*facDiff[2]+0.5*edgeSurf[1]*facDiff[1]+0.5*edgeSurf[0]*facDiff[0]; 
  edgeSurf_incr[1] = 0.5*facDiff[2]*edgeSurf[5]+0.5*edgeSurf[2]*facDiff[3]+0.5*edgeSurf[0]*facDiff[1]+0.5*facDiff[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.5*facDiff[1]*edgeSurf[5]+0.5*edgeSurf[1]*facDiff[3]+0.5*edgeSurf[0]*facDiff[2]+0.5*facDiff[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.5*facDiff[3]*edgeSurf[11]+0.5*facDiff[2]*edgeSurf[7]+0.5*facDiff[1]*edgeSurf[6]+0.5*facDiff[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.5*facDiff[3]*edgeSurf[12]+0.5*facDiff[2]*edgeSurf[9]+0.5*facDiff[1]*edgeSurf[8]+0.5*facDiff[0]*edgeSurf[4]; 
  edgeSurf_incr[5] = 0.5*facDiff[0]*edgeSurf[5]+0.5*edgeSurf[0]*facDiff[3]+0.5*edgeSurf[1]*facDiff[2]+0.5*facDiff[1]*edgeSurf[2]; 
  edgeSurf_incr[6] = 0.5*facDiff[2]*edgeSurf[11]+0.5*facDiff[3]*edgeSurf[7]+0.5*facDiff[0]*edgeSurf[6]+0.5*facDiff[1]*edgeSurf[3]; 
  edgeSurf_incr[7] = 0.5*facDiff[1]*edgeSurf[11]+0.5*facDiff[0]*edgeSurf[7]+0.5*facDiff[3]*edgeSurf[6]+0.5*facDiff[2]*edgeSurf[3]; 
  edgeSurf_incr[8] = 0.5*facDiff[2]*edgeSurf[12]+0.5*facDiff[3]*edgeSurf[9]+0.5*facDiff[0]*edgeSurf[8]+0.5*facDiff[1]*edgeSurf[4]; 
  edgeSurf_incr[9] = 0.5*facDiff[1]*edgeSurf[12]+0.5*facDiff[0]*edgeSurf[9]+0.5*facDiff[3]*edgeSurf[8]+0.5*facDiff[2]*edgeSurf[4]; 
  edgeSurf_incr[10] = 0.5*facDiff[3]*edgeSurf[15]+0.5*facDiff[2]*edgeSurf[14]+0.5*facDiff[1]*edgeSurf[13]+0.5*facDiff[0]*edgeSurf[10]; 
  edgeSurf_incr[11] = 0.5*facDiff[0]*edgeSurf[11]+0.5*facDiff[1]*edgeSurf[7]+0.5*facDiff[2]*edgeSurf[6]+0.5*edgeSurf[3]*facDiff[3]; 
  edgeSurf_incr[12] = 0.5*facDiff[0]*edgeSurf[12]+0.5*facDiff[1]*edgeSurf[9]+0.5*facDiff[2]*edgeSurf[8]+0.5*facDiff[3]*edgeSurf[4]; 
  edgeSurf_incr[13] = 0.5*facDiff[2]*edgeSurf[15]+0.5*facDiff[3]*edgeSurf[14]+0.5*facDiff[0]*edgeSurf[13]+0.5*facDiff[1]*edgeSurf[10]; 
  edgeSurf_incr[14] = 0.5*facDiff[1]*edgeSurf[15]+0.5*facDiff[0]*edgeSurf[14]+0.5*facDiff[3]*edgeSurf[13]+0.5*facDiff[2]*edgeSurf[10]; 
  edgeSurf_incr[15] = 0.5*facDiff[0]*edgeSurf[15]+0.5*facDiff[1]*edgeSurf[14]+0.5*facDiff[2]*edgeSurf[13]+0.5*facDiff[3]*edgeSurf[10]; 
  edgeSurf_incr[16] = 0.5*facDiff[3]*edgeSurf[20]+0.5000000000000001*facDiff[2]*edgeSurf[18]+0.5000000000000001*facDiff[1]*edgeSurf[17]+0.5*facDiff[0]*edgeSurf[16]; 
  edgeSurf_incr[17] = 0.5000000000000001*facDiff[2]*edgeSurf[20]+0.5*facDiff[3]*edgeSurf[18]+0.5*facDiff[0]*edgeSurf[17]+0.5000000000000001*facDiff[1]*edgeSurf[16]; 
  edgeSurf_incr[18] = 0.5000000000000001*facDiff[1]*edgeSurf[20]+0.5*facDiff[0]*edgeSurf[18]+0.5*facDiff[3]*edgeSurf[17]+0.5000000000000001*facDiff[2]*edgeSurf[16]; 
  edgeSurf_incr[19] = 0.5*facDiff[3]*edgeSurf[23]+0.5000000000000001*facDiff[2]*edgeSurf[22]+0.5000000000000001*facDiff[1]*edgeSurf[21]+0.5*facDiff[0]*edgeSurf[19]; 
  edgeSurf_incr[20] = 0.5*facDiff[0]*edgeSurf[20]+0.5000000000000001*facDiff[1]*edgeSurf[18]+0.5000000000000001*facDiff[2]*edgeSurf[17]+0.5*facDiff[3]*edgeSurf[16]; 
  edgeSurf_incr[21] = 0.5000000000000001*facDiff[2]*edgeSurf[23]+0.5*facDiff[3]*edgeSurf[22]+0.5*facDiff[0]*edgeSurf[21]+0.5000000000000001*facDiff[1]*edgeSurf[19]; 
  edgeSurf_incr[22] = 0.5000000000000001*facDiff[1]*edgeSurf[23]+0.5*facDiff[0]*edgeSurf[22]+0.5*facDiff[3]*edgeSurf[21]+0.5000000000000001*facDiff[2]*edgeSurf[19]; 
  edgeSurf_incr[23] = 0.5*facDiff[0]*edgeSurf[23]+0.5000000000000001*facDiff[1]*edgeSurf[22]+0.5000000000000001*facDiff[2]*edgeSurf[21]+0.5*facDiff[3]*edgeSurf[19]; 

  boundSurf_incr[0] = 0.5*facDiff[3]*boundSurf[5]+0.5*boundSurf[2]*facDiff[2]+0.5*boundSurf[1]*facDiff[1]+0.5*boundSurf[0]*facDiff[0]; 
  boundSurf_incr[1] = 0.5*facDiff[2]*boundSurf[5]+0.5*boundSurf[2]*facDiff[3]+0.5*boundSurf[0]*facDiff[1]+0.5*facDiff[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.5*facDiff[1]*boundSurf[5]+0.5*boundSurf[1]*facDiff[3]+0.5*boundSurf[0]*facDiff[2]+0.5*facDiff[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.5*facDiff[3]*boundSurf[11]+0.5*facDiff[2]*boundSurf[7]+0.5*facDiff[1]*boundSurf[6]+0.5*facDiff[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.5*facDiff[3]*boundSurf[12]+0.5*facDiff[2]*boundSurf[9]+0.5*facDiff[1]*boundSurf[8]+0.5*facDiff[0]*boundSurf[4]; 
  boundSurf_incr[5] = 0.5*facDiff[0]*boundSurf[5]+0.5*boundSurf[0]*facDiff[3]+0.5*boundSurf[1]*facDiff[2]+0.5*facDiff[1]*boundSurf[2]; 
  boundSurf_incr[6] = 0.5*facDiff[2]*boundSurf[11]+0.5*facDiff[3]*boundSurf[7]+0.5*facDiff[0]*boundSurf[6]+0.5*facDiff[1]*boundSurf[3]; 
  boundSurf_incr[7] = 0.5*facDiff[1]*boundSurf[11]+0.5*facDiff[0]*boundSurf[7]+0.5*facDiff[3]*boundSurf[6]+0.5*facDiff[2]*boundSurf[3]; 
  boundSurf_incr[8] = 0.5*facDiff[2]*boundSurf[12]+0.5*facDiff[3]*boundSurf[9]+0.5*facDiff[0]*boundSurf[8]+0.5*facDiff[1]*boundSurf[4]; 
  boundSurf_incr[9] = 0.5*facDiff[1]*boundSurf[12]+0.5*facDiff[0]*boundSurf[9]+0.5*facDiff[3]*boundSurf[8]+0.5*facDiff[2]*boundSurf[4]; 
  boundSurf_incr[10] = 0.5*facDiff[3]*boundSurf[15]+0.5*facDiff[2]*boundSurf[14]+0.5*facDiff[1]*boundSurf[13]+0.5*facDiff[0]*boundSurf[10]; 
  boundSurf_incr[11] = 0.5*facDiff[0]*boundSurf[11]+0.5*facDiff[1]*boundSurf[7]+0.5*facDiff[2]*boundSurf[6]+0.5*boundSurf[3]*facDiff[3]; 
  boundSurf_incr[12] = 0.5*facDiff[0]*boundSurf[12]+0.5*facDiff[1]*boundSurf[9]+0.5*facDiff[2]*boundSurf[8]+0.5*facDiff[3]*boundSurf[4]; 
  boundSurf_incr[13] = 0.5*facDiff[2]*boundSurf[15]+0.5*facDiff[3]*boundSurf[14]+0.5*facDiff[0]*boundSurf[13]+0.5*facDiff[1]*boundSurf[10]; 
  boundSurf_incr[14] = 0.5*facDiff[1]*boundSurf[15]+0.5*facDiff[0]*boundSurf[14]+0.5*facDiff[3]*boundSurf[13]+0.5*facDiff[2]*boundSurf[10]; 
  boundSurf_incr[15] = 0.5*facDiff[0]*boundSurf[15]+0.5*facDiff[1]*boundSurf[14]+0.5*facDiff[2]*boundSurf[13]+0.5*facDiff[3]*boundSurf[10]; 
  boundSurf_incr[16] = 0.5*facDiff[3]*boundSurf[20]+0.5000000000000001*facDiff[2]*boundSurf[18]+0.5000000000000001*facDiff[1]*boundSurf[17]+0.5*facDiff[0]*boundSurf[16]; 
  boundSurf_incr[17] = 0.5000000000000001*facDiff[2]*boundSurf[20]+0.5*facDiff[3]*boundSurf[18]+0.5*facDiff[0]*boundSurf[17]+0.5000000000000001*facDiff[1]*boundSurf[16]; 
  boundSurf_incr[18] = 0.5000000000000001*facDiff[1]*boundSurf[20]+0.5*facDiff[0]*boundSurf[18]+0.5*facDiff[3]*boundSurf[17]+0.5000000000000001*facDiff[2]*boundSurf[16]; 
  boundSurf_incr[19] = 0.5*facDiff[3]*boundSurf[23]+0.5000000000000001*facDiff[2]*boundSurf[22]+0.5000000000000001*facDiff[1]*boundSurf[21]+0.5*facDiff[0]*boundSurf[19]; 
  boundSurf_incr[20] = 0.5*facDiff[0]*boundSurf[20]+0.5000000000000001*facDiff[1]*boundSurf[18]+0.5000000000000001*facDiff[2]*boundSurf[17]+0.5*facDiff[3]*boundSurf[16]; 
  boundSurf_incr[21] = 0.5000000000000001*facDiff[2]*boundSurf[23]+0.5*facDiff[3]*boundSurf[22]+0.5*facDiff[0]*boundSurf[21]+0.5000000000000001*facDiff[1]*boundSurf[19]; 
  boundSurf_incr[22] = 0.5000000000000001*facDiff[1]*boundSurf[23]+0.5*facDiff[0]*boundSurf[22]+0.5*facDiff[3]*boundSurf[21]+0.5000000000000001*facDiff[2]*boundSurf[19]; 
  boundSurf_incr[23] = 0.5*facDiff[0]*boundSurf[23]+0.5000000000000001*facDiff[1]*boundSurf[22]+0.5000000000000001*facDiff[2]*boundSurf[21]+0.5*facDiff[3]*boundSurf[19]; 


  } else { 

  double edgeSurf[24] = {0.0}; 
  edgeSurf[0] = 0.5412658773652741*w[3]*fskin[4]-0.270632938682637*dxv[3]*fskin[4]+0.5412658773652741*w[3]*fedge[4]-0.270632938682637*dxv[3]*fedge[4]-0.5625*fskin[0]*w[3]+0.5625*fedge[0]*w[3]+0.28125*fskin[0]*dxv[3]-0.28125*fedge[0]*dxv[3]; 
  edgeSurf[1] = 0.5412658773652741*w[3]*fskin[8]-0.270632938682637*dxv[3]*fskin[8]+0.5412658773652741*w[3]*fedge[8]-0.270632938682637*dxv[3]*fedge[8]-0.5625*fskin[1]*w[3]+0.5625*fedge[1]*w[3]+0.28125*fskin[1]*dxv[3]-0.28125*fedge[1]*dxv[3]; 
  edgeSurf[2] = 0.5412658773652741*w[3]*fskin[9]-0.270632938682637*dxv[3]*fskin[9]+0.5412658773652741*w[3]*fedge[9]-0.270632938682637*dxv[3]*fedge[9]-0.5625*fskin[2]*w[3]+0.5625*fedge[2]*w[3]+0.28125*fskin[2]*dxv[3]-0.28125*fedge[2]*dxv[3]; 
  edgeSurf[3] = 0.5412658773652741*w[3]*fskin[10]-0.270632938682637*dxv[3]*fskin[10]+0.5412658773652741*w[3]*fedge[10]-0.270632938682637*dxv[3]*fedge[10]-0.5625*fskin[3]*w[3]+0.5625*fedge[3]*w[3]+0.28125*dxv[3]*fskin[3]-0.28125*dxv[3]*fedge[3]; 
  edgeSurf[4] = (-1.4375*w[3]*fskin[4])+0.71875*dxv[3]*fskin[4]-0.4375*w[3]*fedge[4]+0.21875*dxv[3]*fedge[4]+1.407291281149713*fskin[0]*w[3]-0.5412658773652741*fedge[0]*w[3]-0.7036456405748563*fskin[0]*dxv[3]+0.270632938682637*fedge[0]*dxv[3]; 
  edgeSurf[5] = 0.5412658773652741*w[3]*fskin[12]-0.270632938682637*dxv[3]*fskin[12]+0.5412658773652741*w[3]*fedge[12]-0.270632938682637*dxv[3]*fedge[12]-0.5625*w[3]*fskin[5]+0.28125*dxv[3]*fskin[5]+0.5625*w[3]*fedge[5]-0.28125*dxv[3]*fedge[5]; 
  edgeSurf[6] = 0.5412658773652741*w[3]*fskin[13]-0.270632938682637*dxv[3]*fskin[13]+0.5412658773652741*w[3]*fedge[13]-0.270632938682637*dxv[3]*fedge[13]-0.5625*w[3]*fskin[6]+0.28125*dxv[3]*fskin[6]+0.5625*w[3]*fedge[6]-0.28125*dxv[3]*fedge[6]; 
  edgeSurf[7] = 0.5412658773652741*w[3]*fskin[14]-0.270632938682637*dxv[3]*fskin[14]+0.5412658773652741*w[3]*fedge[14]-0.270632938682637*dxv[3]*fedge[14]-0.5625*w[3]*fskin[7]+0.28125*dxv[3]*fskin[7]+0.5625*w[3]*fedge[7]-0.28125*dxv[3]*fedge[7]; 
  edgeSurf[8] = (-1.4375*w[3]*fskin[8])+0.71875*dxv[3]*fskin[8]-0.4375*w[3]*fedge[8]+0.21875*dxv[3]*fedge[8]+1.407291281149713*fskin[1]*w[3]-0.5412658773652741*fedge[1]*w[3]-0.7036456405748563*fskin[1]*dxv[3]+0.270632938682637*fedge[1]*dxv[3]; 
  edgeSurf[9] = (-1.4375*w[3]*fskin[9])+0.71875*dxv[3]*fskin[9]-0.4375*w[3]*fedge[9]+0.21875*dxv[3]*fedge[9]+1.407291281149713*fskin[2]*w[3]-0.5412658773652741*fedge[2]*w[3]-0.7036456405748563*fskin[2]*dxv[3]+0.270632938682637*fedge[2]*dxv[3]; 
  edgeSurf[10] = (-1.4375*w[3]*fskin[10])+0.71875*dxv[3]*fskin[10]-0.4375*w[3]*fedge[10]+0.21875*dxv[3]*fedge[10]+1.407291281149713*fskin[3]*w[3]-0.5412658773652741*fedge[3]*w[3]-0.7036456405748563*dxv[3]*fskin[3]+0.270632938682637*dxv[3]*fedge[3]; 
  edgeSurf[11] = 0.5412658773652741*w[3]*fskin[15]-0.270632938682637*dxv[3]*fskin[15]+0.5412658773652741*w[3]*fedge[15]-0.270632938682637*dxv[3]*fedge[15]-0.5625*w[3]*fskin[11]+0.28125*dxv[3]*fskin[11]+0.5625*w[3]*fedge[11]-0.28125*dxv[3]*fedge[11]; 
  edgeSurf[12] = (-1.4375*w[3]*fskin[12])+0.71875*dxv[3]*fskin[12]-0.4375*w[3]*fedge[12]+0.21875*dxv[3]*fedge[12]+1.407291281149713*w[3]*fskin[5]-0.7036456405748563*dxv[3]*fskin[5]-0.5412658773652741*w[3]*fedge[5]+0.270632938682637*dxv[3]*fedge[5]; 
  edgeSurf[13] = (-1.4375*w[3]*fskin[13])+0.71875*dxv[3]*fskin[13]-0.4375*w[3]*fedge[13]+0.21875*dxv[3]*fedge[13]+1.407291281149713*w[3]*fskin[6]-0.7036456405748563*dxv[3]*fskin[6]-0.5412658773652741*w[3]*fedge[6]+0.270632938682637*dxv[3]*fedge[6]; 
  edgeSurf[14] = (-1.4375*w[3]*fskin[14])+0.71875*dxv[3]*fskin[14]-0.4375*w[3]*fedge[14]+0.21875*dxv[3]*fedge[14]+1.407291281149713*w[3]*fskin[7]-0.7036456405748563*dxv[3]*fskin[7]-0.5412658773652741*w[3]*fedge[7]+0.270632938682637*dxv[3]*fedge[7]; 
  edgeSurf[15] = (-1.4375*w[3]*fskin[15])+0.71875*dxv[3]*fskin[15]-0.4375*w[3]*fedge[15]+0.21875*dxv[3]*fedge[15]+1.407291281149713*w[3]*fskin[11]-0.7036456405748563*dxv[3]*fskin[11]-0.5412658773652741*w[3]*fedge[11]+0.270632938682637*dxv[3]*fedge[11]; 
  edgeSurf[16] = 0.5412658773652742*w[3]*fskin[19]-0.2706329386826371*dxv[3]*fskin[19]+0.5412658773652742*w[3]*fedge[19]-0.2706329386826371*dxv[3]*fedge[19]-0.5625*w[3]*fskin[16]+0.28125*dxv[3]*fskin[16]+0.5625*w[3]*fedge[16]-0.28125*dxv[3]*fedge[16]; 
  edgeSurf[17] = 0.5412658773652742*w[3]*fskin[21]-0.2706329386826371*dxv[3]*fskin[21]+0.5412658773652742*w[3]*fedge[21]-0.2706329386826371*dxv[3]*fedge[21]-0.5625*w[3]*fskin[17]+0.28125*dxv[3]*fskin[17]+0.5625*w[3]*fedge[17]-0.28125*dxv[3]*fedge[17]; 
  edgeSurf[18] = 0.5412658773652742*w[3]*fskin[22]-0.2706329386826371*dxv[3]*fskin[22]+0.5412658773652742*w[3]*fedge[22]-0.2706329386826371*dxv[3]*fedge[22]-0.5625*w[3]*fskin[18]+0.28125*dxv[3]*fskin[18]+0.5625*w[3]*fedge[18]-0.28125*dxv[3]*fedge[18]; 
  edgeSurf[19] = (-1.4375*w[3]*fskin[19])+0.71875*dxv[3]*fskin[19]-0.4375*w[3]*fedge[19]+0.21875*dxv[3]*fedge[19]+1.407291281149713*w[3]*fskin[16]-0.7036456405748563*dxv[3]*fskin[16]-0.5412658773652742*w[3]*fedge[16]+0.2706329386826371*dxv[3]*fedge[16]; 
  edgeSurf[20] = 0.5412658773652742*w[3]*fskin[23]-0.2706329386826371*dxv[3]*fskin[23]+0.5412658773652742*w[3]*fedge[23]-0.2706329386826371*dxv[3]*fedge[23]-0.5625*w[3]*fskin[20]+0.28125*dxv[3]*fskin[20]+0.5625*w[3]*fedge[20]-0.28125*dxv[3]*fedge[20]; 
  edgeSurf[21] = (-1.4375*w[3]*fskin[21])+0.71875*dxv[3]*fskin[21]-0.4375*w[3]*fedge[21]+0.21875*dxv[3]*fedge[21]+1.407291281149713*w[3]*fskin[17]-0.7036456405748563*dxv[3]*fskin[17]-0.5412658773652742*w[3]*fedge[17]+0.2706329386826371*dxv[3]*fedge[17]; 
  edgeSurf[22] = (-1.4375*w[3]*fskin[22])+0.71875*dxv[3]*fskin[22]-0.4375*w[3]*fedge[22]+0.21875*dxv[3]*fedge[22]+1.407291281149713*w[3]*fskin[18]-0.7036456405748563*dxv[3]*fskin[18]-0.5412658773652742*w[3]*fedge[18]+0.2706329386826371*dxv[3]*fedge[18]; 
  edgeSurf[23] = (-1.4375*w[3]*fskin[23])+0.71875*dxv[3]*fskin[23]-0.4375*w[3]*fedge[23]+0.21875*dxv[3]*fedge[23]+1.407291281149713*w[3]*fskin[20]-0.7036456405748563*dxv[3]*fskin[20]-0.5412658773652742*w[3]*fedge[20]+0.2706329386826371*dxv[3]*fedge[20]; 

  double boundSurf[24] = {0.0}; 
  boundSurf[4] = (-1.5*w[3]*fskin[4])-0.75*dxv[3]*fskin[4]-0.8660254037844386*fskin[0]*w[3]-0.4330127018922193*fskin[0]*dxv[3]; 
  boundSurf[8] = (-1.5*w[3]*fskin[8])-0.75*dxv[3]*fskin[8]-0.8660254037844386*fskin[1]*w[3]-0.4330127018922193*fskin[1]*dxv[3]; 
  boundSurf[9] = (-1.5*w[3]*fskin[9])-0.75*dxv[3]*fskin[9]-0.8660254037844386*fskin[2]*w[3]-0.4330127018922193*fskin[2]*dxv[3]; 
  boundSurf[10] = (-1.5*w[3]*fskin[10])-0.75*dxv[3]*fskin[10]-0.8660254037844386*fskin[3]*w[3]-0.4330127018922193*dxv[3]*fskin[3]; 
  boundSurf[12] = (-1.5*w[3]*fskin[12])-0.75*dxv[3]*fskin[12]-0.8660254037844386*w[3]*fskin[5]-0.4330127018922193*dxv[3]*fskin[5]; 
  boundSurf[13] = (-1.5*w[3]*fskin[13])-0.75*dxv[3]*fskin[13]-0.8660254037844386*w[3]*fskin[6]-0.4330127018922193*dxv[3]*fskin[6]; 
  boundSurf[14] = (-1.5*w[3]*fskin[14])-0.75*dxv[3]*fskin[14]-0.8660254037844386*w[3]*fskin[7]-0.4330127018922193*dxv[3]*fskin[7]; 
  boundSurf[15] = (-1.5*w[3]*fskin[15])-0.75*dxv[3]*fskin[15]-0.8660254037844386*w[3]*fskin[11]-0.4330127018922193*dxv[3]*fskin[11]; 
  boundSurf[19] = (-1.5*w[3]*fskin[19])-0.75*dxv[3]*fskin[19]-0.8660254037844387*w[3]*fskin[16]-0.4330127018922194*dxv[3]*fskin[16]; 
  boundSurf[21] = (-1.5*w[3]*fskin[21])-0.75*dxv[3]*fskin[21]-0.8660254037844387*w[3]*fskin[17]-0.4330127018922194*dxv[3]*fskin[17]; 
  boundSurf[22] = (-1.5*w[3]*fskin[22])-0.75*dxv[3]*fskin[22]-0.8660254037844387*w[3]*fskin[18]-0.4330127018922194*dxv[3]*fskin[18]; 
  boundSurf[23] = (-1.5*w[3]*fskin[23])-0.75*dxv[3]*fskin[23]-0.8660254037844387*w[3]*fskin[20]-0.4330127018922194*dxv[3]*fskin[20]; 

  edgeSurf_incr[0] = 0.5*facDiff[3]*edgeSurf[5]+0.5*edgeSurf[2]*facDiff[2]+0.5*edgeSurf[1]*facDiff[1]+0.5*edgeSurf[0]*facDiff[0]; 
  edgeSurf_incr[1] = 0.5*facDiff[2]*edgeSurf[5]+0.5*edgeSurf[2]*facDiff[3]+0.5*edgeSurf[0]*facDiff[1]+0.5*facDiff[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.5*facDiff[1]*edgeSurf[5]+0.5*edgeSurf[1]*facDiff[3]+0.5*edgeSurf[0]*facDiff[2]+0.5*facDiff[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.5*facDiff[3]*edgeSurf[11]+0.5*facDiff[2]*edgeSurf[7]+0.5*facDiff[1]*edgeSurf[6]+0.5*facDiff[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.5*facDiff[3]*edgeSurf[12]+0.5*facDiff[2]*edgeSurf[9]+0.5*facDiff[1]*edgeSurf[8]+0.5*facDiff[0]*edgeSurf[4]; 
  edgeSurf_incr[5] = 0.5*facDiff[0]*edgeSurf[5]+0.5*edgeSurf[0]*facDiff[3]+0.5*edgeSurf[1]*facDiff[2]+0.5*facDiff[1]*edgeSurf[2]; 
  edgeSurf_incr[6] = 0.5*facDiff[2]*edgeSurf[11]+0.5*facDiff[3]*edgeSurf[7]+0.5*facDiff[0]*edgeSurf[6]+0.5*facDiff[1]*edgeSurf[3]; 
  edgeSurf_incr[7] = 0.5*facDiff[1]*edgeSurf[11]+0.5*facDiff[0]*edgeSurf[7]+0.5*facDiff[3]*edgeSurf[6]+0.5*facDiff[2]*edgeSurf[3]; 
  edgeSurf_incr[8] = 0.5*facDiff[2]*edgeSurf[12]+0.5*facDiff[3]*edgeSurf[9]+0.5*facDiff[0]*edgeSurf[8]+0.5*facDiff[1]*edgeSurf[4]; 
  edgeSurf_incr[9] = 0.5*facDiff[1]*edgeSurf[12]+0.5*facDiff[0]*edgeSurf[9]+0.5*facDiff[3]*edgeSurf[8]+0.5*facDiff[2]*edgeSurf[4]; 
  edgeSurf_incr[10] = 0.5*facDiff[3]*edgeSurf[15]+0.5*facDiff[2]*edgeSurf[14]+0.5*facDiff[1]*edgeSurf[13]+0.5*facDiff[0]*edgeSurf[10]; 
  edgeSurf_incr[11] = 0.5*facDiff[0]*edgeSurf[11]+0.5*facDiff[1]*edgeSurf[7]+0.5*facDiff[2]*edgeSurf[6]+0.5*edgeSurf[3]*facDiff[3]; 
  edgeSurf_incr[12] = 0.5*facDiff[0]*edgeSurf[12]+0.5*facDiff[1]*edgeSurf[9]+0.5*facDiff[2]*edgeSurf[8]+0.5*facDiff[3]*edgeSurf[4]; 
  edgeSurf_incr[13] = 0.5*facDiff[2]*edgeSurf[15]+0.5*facDiff[3]*edgeSurf[14]+0.5*facDiff[0]*edgeSurf[13]+0.5*facDiff[1]*edgeSurf[10]; 
  edgeSurf_incr[14] = 0.5*facDiff[1]*edgeSurf[15]+0.5*facDiff[0]*edgeSurf[14]+0.5*facDiff[3]*edgeSurf[13]+0.5*facDiff[2]*edgeSurf[10]; 
  edgeSurf_incr[15] = 0.5*facDiff[0]*edgeSurf[15]+0.5*facDiff[1]*edgeSurf[14]+0.5*facDiff[2]*edgeSurf[13]+0.5*facDiff[3]*edgeSurf[10]; 
  edgeSurf_incr[16] = 0.5*facDiff[3]*edgeSurf[20]+0.5000000000000001*facDiff[2]*edgeSurf[18]+0.5000000000000001*facDiff[1]*edgeSurf[17]+0.5*facDiff[0]*edgeSurf[16]; 
  edgeSurf_incr[17] = 0.5000000000000001*facDiff[2]*edgeSurf[20]+0.5*facDiff[3]*edgeSurf[18]+0.5*facDiff[0]*edgeSurf[17]+0.5000000000000001*facDiff[1]*edgeSurf[16]; 
  edgeSurf_incr[18] = 0.5000000000000001*facDiff[1]*edgeSurf[20]+0.5*facDiff[0]*edgeSurf[18]+0.5*facDiff[3]*edgeSurf[17]+0.5000000000000001*facDiff[2]*edgeSurf[16]; 
  edgeSurf_incr[19] = 0.5*facDiff[3]*edgeSurf[23]+0.5000000000000001*facDiff[2]*edgeSurf[22]+0.5000000000000001*facDiff[1]*edgeSurf[21]+0.5*facDiff[0]*edgeSurf[19]; 
  edgeSurf_incr[20] = 0.5*facDiff[0]*edgeSurf[20]+0.5000000000000001*facDiff[1]*edgeSurf[18]+0.5000000000000001*facDiff[2]*edgeSurf[17]+0.5*facDiff[3]*edgeSurf[16]; 
  edgeSurf_incr[21] = 0.5000000000000001*facDiff[2]*edgeSurf[23]+0.5*facDiff[3]*edgeSurf[22]+0.5*facDiff[0]*edgeSurf[21]+0.5000000000000001*facDiff[1]*edgeSurf[19]; 
  edgeSurf_incr[22] = 0.5000000000000001*facDiff[1]*edgeSurf[23]+0.5*facDiff[0]*edgeSurf[22]+0.5*facDiff[3]*edgeSurf[21]+0.5000000000000001*facDiff[2]*edgeSurf[19]; 
  edgeSurf_incr[23] = 0.5*facDiff[0]*edgeSurf[23]+0.5000000000000001*facDiff[1]*edgeSurf[22]+0.5000000000000001*facDiff[2]*edgeSurf[21]+0.5*facDiff[3]*edgeSurf[19]; 

  boundSurf_incr[0] = 0.5*facDiff[3]*boundSurf[5]+0.5*boundSurf[2]*facDiff[2]+0.5*boundSurf[1]*facDiff[1]+0.5*boundSurf[0]*facDiff[0]; 
  boundSurf_incr[1] = 0.5*facDiff[2]*boundSurf[5]+0.5*boundSurf[2]*facDiff[3]+0.5*boundSurf[0]*facDiff[1]+0.5*facDiff[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.5*facDiff[1]*boundSurf[5]+0.5*boundSurf[1]*facDiff[3]+0.5*boundSurf[0]*facDiff[2]+0.5*facDiff[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.5*facDiff[3]*boundSurf[11]+0.5*facDiff[2]*boundSurf[7]+0.5*facDiff[1]*boundSurf[6]+0.5*facDiff[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.5*facDiff[3]*boundSurf[12]+0.5*facDiff[2]*boundSurf[9]+0.5*facDiff[1]*boundSurf[8]+0.5*facDiff[0]*boundSurf[4]; 
  boundSurf_incr[5] = 0.5*facDiff[0]*boundSurf[5]+0.5*boundSurf[0]*facDiff[3]+0.5*boundSurf[1]*facDiff[2]+0.5*facDiff[1]*boundSurf[2]; 
  boundSurf_incr[6] = 0.5*facDiff[2]*boundSurf[11]+0.5*facDiff[3]*boundSurf[7]+0.5*facDiff[0]*boundSurf[6]+0.5*facDiff[1]*boundSurf[3]; 
  boundSurf_incr[7] = 0.5*facDiff[1]*boundSurf[11]+0.5*facDiff[0]*boundSurf[7]+0.5*facDiff[3]*boundSurf[6]+0.5*facDiff[2]*boundSurf[3]; 
  boundSurf_incr[8] = 0.5*facDiff[2]*boundSurf[12]+0.5*facDiff[3]*boundSurf[9]+0.5*facDiff[0]*boundSurf[8]+0.5*facDiff[1]*boundSurf[4]; 
  boundSurf_incr[9] = 0.5*facDiff[1]*boundSurf[12]+0.5*facDiff[0]*boundSurf[9]+0.5*facDiff[3]*boundSurf[8]+0.5*facDiff[2]*boundSurf[4]; 
  boundSurf_incr[10] = 0.5*facDiff[3]*boundSurf[15]+0.5*facDiff[2]*boundSurf[14]+0.5*facDiff[1]*boundSurf[13]+0.5*facDiff[0]*boundSurf[10]; 
  boundSurf_incr[11] = 0.5*facDiff[0]*boundSurf[11]+0.5*facDiff[1]*boundSurf[7]+0.5*facDiff[2]*boundSurf[6]+0.5*boundSurf[3]*facDiff[3]; 
  boundSurf_incr[12] = 0.5*facDiff[0]*boundSurf[12]+0.5*facDiff[1]*boundSurf[9]+0.5*facDiff[2]*boundSurf[8]+0.5*facDiff[3]*boundSurf[4]; 
  boundSurf_incr[13] = 0.5*facDiff[2]*boundSurf[15]+0.5*facDiff[3]*boundSurf[14]+0.5*facDiff[0]*boundSurf[13]+0.5*facDiff[1]*boundSurf[10]; 
  boundSurf_incr[14] = 0.5*facDiff[1]*boundSurf[15]+0.5*facDiff[0]*boundSurf[14]+0.5*facDiff[3]*boundSurf[13]+0.5*facDiff[2]*boundSurf[10]; 
  boundSurf_incr[15] = 0.5*facDiff[0]*boundSurf[15]+0.5*facDiff[1]*boundSurf[14]+0.5*facDiff[2]*boundSurf[13]+0.5*facDiff[3]*boundSurf[10]; 
  boundSurf_incr[16] = 0.5*facDiff[3]*boundSurf[20]+0.5000000000000001*facDiff[2]*boundSurf[18]+0.5000000000000001*facDiff[1]*boundSurf[17]+0.5*facDiff[0]*boundSurf[16]; 
  boundSurf_incr[17] = 0.5000000000000001*facDiff[2]*boundSurf[20]+0.5*facDiff[3]*boundSurf[18]+0.5*facDiff[0]*boundSurf[17]+0.5000000000000001*facDiff[1]*boundSurf[16]; 
  boundSurf_incr[18] = 0.5000000000000001*facDiff[1]*boundSurf[20]+0.5*facDiff[0]*boundSurf[18]+0.5*facDiff[3]*boundSurf[17]+0.5000000000000001*facDiff[2]*boundSurf[16]; 
  boundSurf_incr[19] = 0.5*facDiff[3]*boundSurf[23]+0.5000000000000001*facDiff[2]*boundSurf[22]+0.5000000000000001*facDiff[1]*boundSurf[21]+0.5*facDiff[0]*boundSurf[19]; 
  boundSurf_incr[20] = 0.5*facDiff[0]*boundSurf[20]+0.5000000000000001*facDiff[1]*boundSurf[18]+0.5000000000000001*facDiff[2]*boundSurf[17]+0.5*facDiff[3]*boundSurf[16]; 
  boundSurf_incr[21] = 0.5000000000000001*facDiff[2]*boundSurf[23]+0.5*facDiff[3]*boundSurf[22]+0.5*facDiff[0]*boundSurf[21]+0.5000000000000001*facDiff[1]*boundSurf[19]; 
  boundSurf_incr[22] = 0.5000000000000001*facDiff[1]*boundSurf[23]+0.5*facDiff[0]*boundSurf[22]+0.5*facDiff[3]*boundSurf[21]+0.5000000000000001*facDiff[2]*boundSurf[19]; 
  boundSurf_incr[23] = 0.5*facDiff[0]*boundSurf[23]+0.5000000000000001*facDiff[1]*boundSurf[22]+0.5000000000000001*facDiff[2]*boundSurf[21]+0.5*facDiff[3]*boundSurf[19]; 

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
} 
