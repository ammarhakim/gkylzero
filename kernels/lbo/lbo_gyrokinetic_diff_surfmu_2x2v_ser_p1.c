#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfmu_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // w[4]: cell-center coordinates. 
  // dxv[4]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 
  double incr[24] = {0.0}; 

  double diffFac[4] = {0.}; 
  diffFac[0] = bmag_inv[3]*nuVtSqSum[3]*m_+bmag_inv[2]*nuVtSqSum[2]*m_+bmag_inv[1]*nuVtSqSum[1]*m_+bmag_inv[0]*nuVtSqSum[0]*m_; 
  diffFac[1] = bmag_inv[2]*nuVtSqSum[3]*m_+nuVtSqSum[2]*bmag_inv[3]*m_+bmag_inv[0]*nuVtSqSum[1]*m_+nuVtSqSum[0]*bmag_inv[1]*m_; 
  diffFac[2] = bmag_inv[1]*nuVtSqSum[3]*m_+nuVtSqSum[1]*bmag_inv[3]*m_+bmag_inv[0]*nuVtSqSum[2]*m_+nuVtSqSum[0]*bmag_inv[2]*m_; 
  diffFac[3] = bmag_inv[0]*nuVtSqSum[3]*m_+nuVtSqSum[0]*bmag_inv[3]*m_+bmag_inv[1]*nuVtSqSum[2]*m_+nuVtSqSum[1]*bmag_inv[2]*m_; 

  double f_xx[24] = {0.0}; 
  f_xx[0] = (-0.5412658773652741*w[3]*fr[4])-0.270632938682637*dxv[3]*fr[4]+0.5412658773652741*w[3]*fl[4]-0.270632938682637*dxv[3]*fl[4]-0.5412658773652741*dxv[3]*fc[4]+0.5625*fr[0]*w[3]+0.5625*fl[0]*w[3]-1.125*fc[0]*w[3]+0.28125*fr[0]*dxv[3]-0.28125*fl[0]*dxv[3]; 
  f_xx[1] = (-0.5412658773652741*w[3]*fr[8])-0.270632938682637*dxv[3]*fr[8]+0.5412658773652741*w[3]*fl[8]-0.270632938682637*dxv[3]*fl[8]-0.5412658773652741*dxv[3]*fc[8]+0.5625*fr[1]*w[3]+0.5625*fl[1]*w[3]-1.125*fc[1]*w[3]+0.28125*fr[1]*dxv[3]-0.28125*fl[1]*dxv[3]; 
  f_xx[2] = (-0.5412658773652741*w[3]*fr[9])-0.270632938682637*dxv[3]*fr[9]+0.5412658773652741*w[3]*fl[9]-0.270632938682637*dxv[3]*fl[9]-0.5412658773652741*dxv[3]*fc[9]+0.5625*fr[2]*w[3]+0.5625*fl[2]*w[3]-1.125*fc[2]*w[3]+0.28125*fr[2]*dxv[3]-0.28125*fl[2]*dxv[3]; 
  f_xx[3] = (-0.5412658773652741*w[3]*fr[10])-0.270632938682637*dxv[3]*fr[10]+0.5412658773652741*w[3]*fl[10]-0.270632938682637*dxv[3]*fl[10]-0.5412658773652741*dxv[3]*fc[10]+0.5625*fr[3]*w[3]+0.5625*fl[3]*w[3]-1.125*fc[3]*w[3]+0.28125*dxv[3]*fr[3]-0.28125*dxv[3]*fl[3]; 
  f_xx[4] = (-0.4375*w[3]*fr[4])-0.21875*dxv[3]*fr[4]-0.4375*w[3]*fl[4]+0.21875*dxv[3]*fl[4]-2.875*w[3]*fc[4]+0.5412658773652741*fr[0]*w[3]-0.5412658773652741*fl[0]*w[3]+0.270632938682637*fr[0]*dxv[3]+0.270632938682637*fl[0]*dxv[3]-0.5412658773652741*fc[0]*dxv[3]; 
  f_xx[5] = (-0.5412658773652741*w[3]*fr[12])-0.270632938682637*dxv[3]*fr[12]+0.5412658773652741*w[3]*fl[12]-0.270632938682637*dxv[3]*fl[12]-0.5412658773652741*dxv[3]*fc[12]+0.5625*w[3]*fr[5]+0.28125*dxv[3]*fr[5]+0.5625*w[3]*fl[5]-0.28125*dxv[3]*fl[5]-1.125*w[3]*fc[5]; 
  f_xx[6] = (-0.5412658773652741*w[3]*fr[13])-0.270632938682637*dxv[3]*fr[13]+0.5412658773652741*w[3]*fl[13]-0.270632938682637*dxv[3]*fl[13]-0.5412658773652741*dxv[3]*fc[13]+0.5625*w[3]*fr[6]+0.28125*dxv[3]*fr[6]+0.5625*w[3]*fl[6]-0.28125*dxv[3]*fl[6]-1.125*w[3]*fc[6]; 
  f_xx[7] = (-0.5412658773652741*w[3]*fr[14])-0.270632938682637*dxv[3]*fr[14]+0.5412658773652741*w[3]*fl[14]-0.270632938682637*dxv[3]*fl[14]-0.5412658773652741*dxv[3]*fc[14]+0.5625*w[3]*fr[7]+0.28125*dxv[3]*fr[7]+0.5625*w[3]*fl[7]-0.28125*dxv[3]*fl[7]-1.125*w[3]*fc[7]; 
  f_xx[8] = (-0.4375*w[3]*fr[8])-0.21875*dxv[3]*fr[8]-0.4375*w[3]*fl[8]+0.21875*dxv[3]*fl[8]-2.875*w[3]*fc[8]+0.5412658773652741*fr[1]*w[3]-0.5412658773652741*fl[1]*w[3]+0.270632938682637*fr[1]*dxv[3]+0.270632938682637*fl[1]*dxv[3]-0.5412658773652741*fc[1]*dxv[3]; 
  f_xx[9] = (-0.4375*w[3]*fr[9])-0.21875*dxv[3]*fr[9]-0.4375*w[3]*fl[9]+0.21875*dxv[3]*fl[9]-2.875*w[3]*fc[9]+0.5412658773652741*fr[2]*w[3]-0.5412658773652741*fl[2]*w[3]+0.270632938682637*fr[2]*dxv[3]+0.270632938682637*fl[2]*dxv[3]-0.5412658773652741*fc[2]*dxv[3]; 
  f_xx[10] = (-0.4375*w[3]*fr[10])-0.21875*dxv[3]*fr[10]-0.4375*w[3]*fl[10]+0.21875*dxv[3]*fl[10]-2.875*w[3]*fc[10]+0.5412658773652741*fr[3]*w[3]-0.5412658773652741*fl[3]*w[3]+0.270632938682637*dxv[3]*fr[3]+0.270632938682637*dxv[3]*fl[3]-0.5412658773652741*dxv[3]*fc[3]; 
  f_xx[11] = (-0.5412658773652741*w[3]*fr[15])-0.270632938682637*dxv[3]*fr[15]+0.5412658773652741*w[3]*fl[15]-0.270632938682637*dxv[3]*fl[15]-0.5412658773652741*dxv[3]*fc[15]+0.5625*w[3]*fr[11]+0.28125*dxv[3]*fr[11]+0.5625*w[3]*fl[11]-0.28125*dxv[3]*fl[11]-1.125*w[3]*fc[11]; 
  f_xx[12] = (-0.4375*w[3]*fr[12])-0.21875*dxv[3]*fr[12]-0.4375*w[3]*fl[12]+0.21875*dxv[3]*fl[12]-2.875*w[3]*fc[12]+0.5412658773652741*w[3]*fr[5]+0.270632938682637*dxv[3]*fr[5]-0.5412658773652741*w[3]*fl[5]+0.270632938682637*dxv[3]*fl[5]-0.5412658773652741*dxv[3]*fc[5]; 
  f_xx[13] = (-0.4375*w[3]*fr[13])-0.21875*dxv[3]*fr[13]-0.4375*w[3]*fl[13]+0.21875*dxv[3]*fl[13]-2.875*w[3]*fc[13]+0.5412658773652741*w[3]*fr[6]+0.270632938682637*dxv[3]*fr[6]-0.5412658773652741*w[3]*fl[6]+0.270632938682637*dxv[3]*fl[6]-0.5412658773652741*dxv[3]*fc[6]; 
  f_xx[14] = (-0.4375*w[3]*fr[14])-0.21875*dxv[3]*fr[14]-0.4375*w[3]*fl[14]+0.21875*dxv[3]*fl[14]-2.875*w[3]*fc[14]+0.5412658773652741*w[3]*fr[7]+0.270632938682637*dxv[3]*fr[7]-0.5412658773652741*w[3]*fl[7]+0.270632938682637*dxv[3]*fl[7]-0.5412658773652741*dxv[3]*fc[7]; 
  f_xx[15] = (-0.4375*w[3]*fr[15])-0.21875*dxv[3]*fr[15]-0.4375*w[3]*fl[15]+0.21875*dxv[3]*fl[15]-2.875*w[3]*fc[15]+0.5412658773652741*w[3]*fr[11]+0.270632938682637*dxv[3]*fr[11]-0.5412658773652741*w[3]*fl[11]+0.270632938682637*dxv[3]*fl[11]-0.5412658773652741*dxv[3]*fc[11]; 
  f_xx[16] = (-0.5412658773652742*w[3]*fr[19])-0.2706329386826371*dxv[3]*fr[19]+0.5412658773652742*w[3]*fl[19]-0.2706329386826371*dxv[3]*fl[19]-0.5412658773652742*dxv[3]*fc[19]+0.5625*w[3]*fr[16]+0.28125*dxv[3]*fr[16]+0.5625*w[3]*fl[16]-0.28125*dxv[3]*fl[16]-1.125*w[3]*fc[16]; 
  f_xx[17] = (-0.5412658773652742*w[3]*fr[21])-0.2706329386826371*dxv[3]*fr[21]+0.5412658773652742*w[3]*fl[21]-0.2706329386826371*dxv[3]*fl[21]-0.5412658773652742*dxv[3]*fc[21]+0.5625*w[3]*fr[17]+0.28125*dxv[3]*fr[17]+0.5625*w[3]*fl[17]-0.28125*dxv[3]*fl[17]-1.125*w[3]*fc[17]; 
  f_xx[18] = (-0.5412658773652742*w[3]*fr[22])-0.2706329386826371*dxv[3]*fr[22]+0.5412658773652742*w[3]*fl[22]-0.2706329386826371*dxv[3]*fl[22]-0.5412658773652742*dxv[3]*fc[22]+0.5625*w[3]*fr[18]+0.28125*dxv[3]*fr[18]+0.5625*w[3]*fl[18]-0.28125*dxv[3]*fl[18]-1.125*w[3]*fc[18]; 
  f_xx[19] = (-0.4375*w[3]*fr[19])-0.21875*dxv[3]*fr[19]-0.4375*w[3]*fl[19]+0.21875*dxv[3]*fl[19]-2.875*w[3]*fc[19]+0.5412658773652742*w[3]*fr[16]+0.2706329386826371*dxv[3]*fr[16]-0.5412658773652742*w[3]*fl[16]+0.2706329386826371*dxv[3]*fl[16]-0.5412658773652742*dxv[3]*fc[16]; 
  f_xx[20] = (-0.5412658773652742*w[3]*fr[23])-0.2706329386826371*dxv[3]*fr[23]+0.5412658773652742*w[3]*fl[23]-0.2706329386826371*dxv[3]*fl[23]-0.5412658773652742*dxv[3]*fc[23]+0.5625*w[3]*fr[20]+0.28125*dxv[3]*fr[20]+0.5625*w[3]*fl[20]-0.28125*dxv[3]*fl[20]-1.125*w[3]*fc[20]; 
  f_xx[21] = (-0.4375*w[3]*fr[21])-0.21875*dxv[3]*fr[21]-0.4375*w[3]*fl[21]+0.21875*dxv[3]*fl[21]-2.875*w[3]*fc[21]+0.5412658773652742*w[3]*fr[17]+0.2706329386826371*dxv[3]*fr[17]-0.5412658773652742*w[3]*fl[17]+0.2706329386826371*dxv[3]*fl[17]-0.5412658773652742*dxv[3]*fc[17]; 
  f_xx[22] = (-0.4375*w[3]*fr[22])-0.21875*dxv[3]*fr[22]-0.4375*w[3]*fl[22]+0.21875*dxv[3]*fl[22]-2.875*w[3]*fc[22]+0.5412658773652742*w[3]*fr[18]+0.2706329386826371*dxv[3]*fr[18]-0.5412658773652742*w[3]*fl[18]+0.2706329386826371*dxv[3]*fl[18]-0.5412658773652742*dxv[3]*fc[18]; 
  f_xx[23] = (-0.4375*w[3]*fr[23])-0.21875*dxv[3]*fr[23]-0.4375*w[3]*fl[23]+0.21875*dxv[3]*fl[23]-2.875*w[3]*fc[23]+0.5412658773652742*w[3]*fr[20]+0.2706329386826371*dxv[3]*fr[20]-0.5412658773652742*w[3]*fl[20]+0.2706329386826371*dxv[3]*fl[20]-0.5412658773652742*dxv[3]*fc[20]; 

  incr[0] = 0.5*diffFac[3]*f_xx[5]+0.5*diffFac[2]*f_xx[2]+0.5*diffFac[1]*f_xx[1]+0.5*diffFac[0]*f_xx[0]; 
  incr[1] = 0.5*diffFac[2]*f_xx[5]+0.5*f_xx[2]*diffFac[3]+0.5*diffFac[0]*f_xx[1]+0.5*f_xx[0]*diffFac[1]; 
  incr[2] = 0.5*diffFac[1]*f_xx[5]+0.5*f_xx[1]*diffFac[3]+0.5*diffFac[0]*f_xx[2]+0.5*f_xx[0]*diffFac[2]; 
  incr[3] = 0.5*diffFac[3]*f_xx[11]+0.5*diffFac[2]*f_xx[7]+0.5*diffFac[1]*f_xx[6]+0.5*diffFac[0]*f_xx[3]; 
  incr[4] = 0.5*diffFac[3]*f_xx[12]+0.5*diffFac[2]*f_xx[9]+0.5*diffFac[1]*f_xx[8]+0.5*diffFac[0]*f_xx[4]; 
  incr[5] = 0.5*diffFac[0]*f_xx[5]+0.5*f_xx[0]*diffFac[3]+0.5*diffFac[1]*f_xx[2]+0.5*f_xx[1]*diffFac[2]; 
  incr[6] = 0.5*diffFac[2]*f_xx[11]+0.5*diffFac[3]*f_xx[7]+0.5*diffFac[0]*f_xx[6]+0.5*diffFac[1]*f_xx[3]; 
  incr[7] = 0.5*diffFac[1]*f_xx[11]+0.5*diffFac[0]*f_xx[7]+0.5*diffFac[3]*f_xx[6]+0.5*diffFac[2]*f_xx[3]; 
  incr[8] = 0.5*diffFac[2]*f_xx[12]+0.5*diffFac[3]*f_xx[9]+0.5*diffFac[0]*f_xx[8]+0.5*diffFac[1]*f_xx[4]; 
  incr[9] = 0.5*diffFac[1]*f_xx[12]+0.5*diffFac[0]*f_xx[9]+0.5*diffFac[3]*f_xx[8]+0.5*diffFac[2]*f_xx[4]; 
  incr[10] = 0.5*diffFac[3]*f_xx[15]+0.5*diffFac[2]*f_xx[14]+0.5*diffFac[1]*f_xx[13]+0.5*diffFac[0]*f_xx[10]; 
  incr[11] = 0.5*diffFac[0]*f_xx[11]+0.5*diffFac[1]*f_xx[7]+0.5*diffFac[2]*f_xx[6]+0.5*diffFac[3]*f_xx[3]; 
  incr[12] = 0.5*diffFac[0]*f_xx[12]+0.5*diffFac[1]*f_xx[9]+0.5*diffFac[2]*f_xx[8]+0.5*diffFac[3]*f_xx[4]; 
  incr[13] = 0.5*diffFac[2]*f_xx[15]+0.5*diffFac[3]*f_xx[14]+0.5*diffFac[0]*f_xx[13]+0.5*diffFac[1]*f_xx[10]; 
  incr[14] = 0.5*diffFac[1]*f_xx[15]+0.5*diffFac[0]*f_xx[14]+0.5*diffFac[3]*f_xx[13]+0.5*diffFac[2]*f_xx[10]; 
  incr[15] = 0.5*diffFac[0]*f_xx[15]+0.5*diffFac[1]*f_xx[14]+0.5*diffFac[2]*f_xx[13]+0.5*diffFac[3]*f_xx[10]; 
  incr[16] = 0.5*diffFac[3]*f_xx[20]+0.5000000000000001*diffFac[2]*f_xx[18]+0.5000000000000001*diffFac[1]*f_xx[17]+0.5*diffFac[0]*f_xx[16]; 
  incr[17] = 0.5000000000000001*diffFac[2]*f_xx[20]+0.5*diffFac[3]*f_xx[18]+0.5*diffFac[0]*f_xx[17]+0.5000000000000001*diffFac[1]*f_xx[16]; 
  incr[18] = 0.5000000000000001*diffFac[1]*f_xx[20]+0.5*diffFac[0]*f_xx[18]+0.5*diffFac[3]*f_xx[17]+0.5000000000000001*diffFac[2]*f_xx[16]; 
  incr[19] = 0.5*diffFac[3]*f_xx[23]+0.5000000000000001*diffFac[2]*f_xx[22]+0.5000000000000001*diffFac[1]*f_xx[21]+0.5*diffFac[0]*f_xx[19]; 
  incr[20] = 0.5*diffFac[0]*f_xx[20]+0.5000000000000001*diffFac[1]*f_xx[18]+0.5000000000000001*diffFac[2]*f_xx[17]+0.5*diffFac[3]*f_xx[16]; 
  incr[21] = 0.5000000000000001*diffFac[2]*f_xx[23]+0.5*diffFac[3]*f_xx[22]+0.5*diffFac[0]*f_xx[21]+0.5000000000000001*diffFac[1]*f_xx[19]; 
  incr[22] = 0.5000000000000001*diffFac[1]*f_xx[23]+0.5*diffFac[0]*f_xx[22]+0.5*diffFac[3]*f_xx[21]+0.5000000000000001*diffFac[2]*f_xx[19]; 
  incr[23] = 0.5*diffFac[0]*f_xx[23]+0.5000000000000001*diffFac[1]*f_xx[22]+0.5000000000000001*diffFac[2]*f_xx[21]+0.5*diffFac[3]*f_xx[19]; 

  out[0] += incr[0]*rdvSq4; 
  out[1] += incr[1]*rdvSq4; 
  out[2] += incr[2]*rdvSq4; 
  out[3] += incr[3]*rdvSq4; 
  out[4] += incr[4]*rdvSq4; 
  out[5] += incr[5]*rdvSq4; 
  out[6] += incr[6]*rdvSq4; 
  out[7] += incr[7]*rdvSq4; 
  out[8] += incr[8]*rdvSq4; 
  out[9] += incr[9]*rdvSq4; 
  out[10] += incr[10]*rdvSq4; 
  out[11] += incr[11]*rdvSq4; 
  out[12] += incr[12]*rdvSq4; 
  out[13] += incr[13]*rdvSq4; 
  out[14] += incr[14]*rdvSq4; 
  out[15] += incr[15]*rdvSq4; 
  out[16] += incr[16]*rdvSq4; 
  out[17] += incr[17]*rdvSq4; 
  out[18] += incr[18]*rdvSq4; 
  out[19] += incr[19]*rdvSq4; 
  out[20] += incr[20]*rdvSq4; 
  out[21] += incr[21]*rdvSq4; 
  out[22] += incr[22]*rdvSq4; 
  out[23] += incr[23]*rdvSq4; 

  return 0.;

} 
