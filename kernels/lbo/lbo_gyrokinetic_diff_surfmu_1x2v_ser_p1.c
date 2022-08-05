#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // m_:            species mass.
  // bmag_inv:      1/(magnetic field magnitude). 
  // w[3]:         cell-center coordinates. 
  // dxv[3]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 
  double temp_diff[8] = {0.0}; 
  double diff_incr[8] = {0.0}; 

  temp_diff[0] = (-0.5412658773652741*w[2]*fr[3])-0.270632938682637*dxv[2]*fr[3]+0.5412658773652741*w[2]*fl[3]-0.270632938682637*dxv[2]*fl[3]-0.5412658773652741*dxv[2]*fc[3]+0.5625*fr[0]*w[2]+0.5625*fl[0]*w[2]-1.125*fc[0]*w[2]+0.28125*fr[0]*dxv[2]-0.28125*fl[0]*dxv[2]; 
  temp_diff[1] = (-0.5412658773652741*w[2]*fr[5])-0.270632938682637*dxv[2]*fr[5]+0.5412658773652741*w[2]*fl[5]-0.270632938682637*dxv[2]*fl[5]-0.5412658773652741*dxv[2]*fc[5]+0.5625*fr[1]*w[2]+0.5625*fl[1]*w[2]-1.125*fc[1]*w[2]+0.28125*fr[1]*dxv[2]-0.28125*fl[1]*dxv[2]; 
  temp_diff[2] = (-0.5412658773652741*w[2]*fr[6])-0.270632938682637*dxv[2]*fr[6]+0.5412658773652741*w[2]*fl[6]-0.270632938682637*dxv[2]*fl[6]-0.5412658773652741*dxv[2]*fc[6]+0.5625*fr[2]*w[2]+0.5625*fl[2]*w[2]-1.125*fc[2]*w[2]+0.28125*dxv[2]*fr[2]-0.28125*dxv[2]*fl[2]; 
  temp_diff[3] = (-0.4375*w[2]*fr[3])-0.21875*dxv[2]*fr[3]-0.4375*w[2]*fl[3]+0.21875*dxv[2]*fl[3]-2.875*w[2]*fc[3]+0.5412658773652741*fr[0]*w[2]-0.5412658773652741*fl[0]*w[2]+0.270632938682637*fr[0]*dxv[2]+0.270632938682637*fl[0]*dxv[2]-0.5412658773652741*fc[0]*dxv[2]; 
  temp_diff[4] = (-0.5412658773652741*w[2]*fr[7])-0.270632938682637*dxv[2]*fr[7]+0.5412658773652741*w[2]*fl[7]-0.270632938682637*dxv[2]*fl[7]-0.5412658773652741*dxv[2]*fc[7]+0.5625*w[2]*fr[4]+0.28125*dxv[2]*fr[4]+0.5625*w[2]*fl[4]-0.28125*dxv[2]*fl[4]-1.125*w[2]*fc[4]; 
  temp_diff[5] = (-0.4375*w[2]*fr[5])-0.21875*dxv[2]*fr[5]-0.4375*w[2]*fl[5]+0.21875*dxv[2]*fl[5]-2.875*w[2]*fc[5]+0.5412658773652741*fr[1]*w[2]-0.5412658773652741*fl[1]*w[2]+0.270632938682637*fr[1]*dxv[2]+0.270632938682637*fl[1]*dxv[2]-0.5412658773652741*fc[1]*dxv[2]; 
  temp_diff[6] = (-0.4375*w[2]*fr[6])-0.21875*dxv[2]*fr[6]-0.4375*w[2]*fl[6]+0.21875*dxv[2]*fl[6]-2.875*w[2]*fc[6]+0.5412658773652741*fr[2]*w[2]-0.5412658773652741*fl[2]*w[2]+0.270632938682637*dxv[2]*fr[2]+0.270632938682637*dxv[2]*fl[2]-0.5412658773652741*dxv[2]*fc[2]; 
  temp_diff[7] = (-0.4375*w[2]*fr[7])-0.21875*dxv[2]*fr[7]-0.4375*w[2]*fl[7]+0.21875*dxv[2]*fl[7]-2.875*w[2]*fc[7]+0.5412658773652741*w[2]*fr[4]+0.270632938682637*dxv[2]*fr[4]-0.5412658773652741*w[2]*fl[4]+0.270632938682637*dxv[2]*fl[4]-0.5412658773652741*dxv[2]*fc[4]; 

  diff_incr[0] = bmag_inv[0]*nuVtSqSum[1]*temp_diff[1]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[1]*m_+temp_diff[0]*bmag_inv[1]*nuVtSqSum[1]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[0]*m_; 
  diff_incr[1] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[1]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[1]*m_+bmag_inv[0]*temp_diff[0]*nuVtSqSum[1]*m_+nuVtSqSum[0]*temp_diff[0]*bmag_inv[1]*m_; 
  diff_incr[2] = bmag_inv[0]*nuVtSqSum[1]*temp_diff[4]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[4]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[2]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[2]*m_; 
  diff_incr[3] = bmag_inv[0]*nuVtSqSum[1]*temp_diff[5]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[5]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[3]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[3]*m_; 
  diff_incr[4] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[4]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[4]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[2]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[2]*m_; 
  diff_incr[5] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[5]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[5]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[3]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[3]*m_; 
  diff_incr[6] = bmag_inv[0]*nuVtSqSum[1]*temp_diff[7]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[7]*m_+bmag_inv[1]*nuVtSqSum[1]*temp_diff[6]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[6]*m_; 
  diff_incr[7] = 1.8*bmag_inv[1]*nuVtSqSum[1]*temp_diff[7]*m_+bmag_inv[0]*nuVtSqSum[0]*temp_diff[7]*m_+bmag_inv[0]*nuVtSqSum[1]*temp_diff[6]*m_+nuVtSqSum[0]*bmag_inv[1]*temp_diff[6]*m_; 

  out[0] += diff_incr[0]*rdvSq4; 
  out[1] += diff_incr[1]*rdvSq4; 
  out[2] += diff_incr[2]*rdvSq4; 
  out[3] += diff_incr[3]*rdvSq4; 
  out[4] += diff_incr[4]*rdvSq4; 
  out[5] += diff_incr[5]*rdvSq4; 
  out[6] += diff_incr[6]*rdvSq4; 
  out[7] += diff_incr[7]*rdvSq4; 
} 
