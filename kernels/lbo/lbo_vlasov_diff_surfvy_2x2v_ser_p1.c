#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_diff_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         cell-center coordinates. 
  // dxv[4]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[8]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 
  double temp_diff[16] = {0.0}; 
  double diff_incr[16] = {0.0}; 

  temp_diff[0] = (-0.5412658773652741*fr[4])+0.5412658773652741*fl[4]+0.5625*fr[0]+0.5625*fl[0]-1.125*fc[0]; 
  temp_diff[1] = (-0.5412658773652741*fr[8])+0.5412658773652741*fl[8]+0.5625*fr[1]+0.5625*fl[1]-1.125*fc[1]; 
  temp_diff[2] = (-0.5412658773652741*fr[9])+0.5412658773652741*fl[9]+0.5625*fr[2]+0.5625*fl[2]-1.125*fc[2]; 
  temp_diff[3] = (-0.5412658773652741*fr[10])+0.5412658773652741*fl[10]+0.5625*fr[3]+0.5625*fl[3]-1.125*fc[3]; 
  temp_diff[4] = (-0.4375*fr[4])-0.4375*fl[4]-2.875*fc[4]+0.5412658773652741*fr[0]-0.5412658773652741*fl[0]; 
  temp_diff[5] = (-0.5412658773652741*fr[12])+0.5412658773652741*fl[12]+0.5625*fr[5]+0.5625*fl[5]-1.125*fc[5]; 
  temp_diff[6] = (-0.5412658773652741*fr[13])+0.5412658773652741*fl[13]+0.5625*fr[6]+0.5625*fl[6]-1.125*fc[6]; 
  temp_diff[7] = (-0.5412658773652741*fr[14])+0.5412658773652741*fl[14]+0.5625*fr[7]+0.5625*fl[7]-1.125*fc[7]; 
  temp_diff[8] = (-0.4375*fr[8])-0.4375*fl[8]-2.875*fc[8]+0.5412658773652741*fr[1]-0.5412658773652741*fl[1]; 
  temp_diff[9] = (-0.4375*fr[9])-0.4375*fl[9]-2.875*fc[9]+0.5412658773652741*fr[2]-0.5412658773652741*fl[2]; 
  temp_diff[10] = (-0.4375*fr[10])-0.4375*fl[10]-2.875*fc[10]+0.5412658773652741*fr[3]-0.5412658773652741*fl[3]; 
  temp_diff[11] = (-0.5412658773652741*fr[15])+0.5412658773652741*fl[15]+0.5625*fr[11]+0.5625*fl[11]-1.125*fc[11]; 
  temp_diff[12] = (-0.4375*fr[12])-0.4375*fl[12]-2.875*fc[12]+0.5412658773652741*fr[5]-0.5412658773652741*fl[5]; 
  temp_diff[13] = (-0.4375*fr[13])-0.4375*fl[13]-2.875*fc[13]+0.5412658773652741*fr[6]-0.5412658773652741*fl[6]; 
  temp_diff[14] = (-0.4375*fr[14])-0.4375*fl[14]-2.875*fc[14]+0.5412658773652741*fr[7]-0.5412658773652741*fl[7]; 
  temp_diff[15] = (-0.4375*fr[15])-0.4375*fl[15]-2.875*fc[15]+0.5412658773652741*fr[11]-0.5412658773652741*fl[11]; 

  diff_incr[0] = 0.5*nuVtSqSum[3]*temp_diff[5]+0.5*nuVtSqSum[2]*temp_diff[2]+0.5*nuVtSqSum[1]*temp_diff[1]+0.5*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.5*nuVtSqSum[2]*temp_diff[5]+0.5*temp_diff[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[1]+0.5*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.5*nuVtSqSum[1]*temp_diff[5]+0.5*temp_diff[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[2]+0.5*temp_diff[0]*nuVtSqSum[2]; 
  diff_incr[3] = 0.5*nuVtSqSum[3]*temp_diff[11]+0.5*nuVtSqSum[2]*temp_diff[7]+0.5*nuVtSqSum[1]*temp_diff[6]+0.5*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.5*nuVtSqSum[3]*temp_diff[12]+0.5*nuVtSqSum[2]*temp_diff[9]+0.5*nuVtSqSum[1]*temp_diff[8]+0.5*nuVtSqSum[0]*temp_diff[4]; 
  diff_incr[5] = 0.5*nuVtSqSum[0]*temp_diff[5]+0.5*temp_diff[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_diff[2]+0.5*temp_diff[1]*nuVtSqSum[2]; 
  diff_incr[6] = 0.5*nuVtSqSum[2]*temp_diff[11]+0.5*nuVtSqSum[3]*temp_diff[7]+0.5*nuVtSqSum[0]*temp_diff[6]+0.5*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[7] = 0.5*nuVtSqSum[1]*temp_diff[11]+0.5*nuVtSqSum[0]*temp_diff[7]+0.5*nuVtSqSum[3]*temp_diff[6]+0.5*nuVtSqSum[2]*temp_diff[3]; 
  diff_incr[8] = 0.5*nuVtSqSum[2]*temp_diff[12]+0.5*nuVtSqSum[3]*temp_diff[9]+0.5*nuVtSqSum[0]*temp_diff[8]+0.5*nuVtSqSum[1]*temp_diff[4]; 
  diff_incr[9] = 0.5*nuVtSqSum[1]*temp_diff[12]+0.5*nuVtSqSum[0]*temp_diff[9]+0.5*nuVtSqSum[3]*temp_diff[8]+0.5*nuVtSqSum[2]*temp_diff[4]; 
  diff_incr[10] = 0.5*nuVtSqSum[3]*temp_diff[15]+0.5*nuVtSqSum[2]*temp_diff[14]+0.5*nuVtSqSum[1]*temp_diff[13]+0.5*nuVtSqSum[0]*temp_diff[10]; 
  diff_incr[11] = 0.5*nuVtSqSum[0]*temp_diff[11]+0.5*nuVtSqSum[1]*temp_diff[7]+0.5*nuVtSqSum[2]*temp_diff[6]+0.5*nuVtSqSum[3]*temp_diff[3]; 
  diff_incr[12] = 0.5*nuVtSqSum[0]*temp_diff[12]+0.5*nuVtSqSum[1]*temp_diff[9]+0.5*nuVtSqSum[2]*temp_diff[8]+0.5*nuVtSqSum[3]*temp_diff[4]; 
  diff_incr[13] = 0.5*nuVtSqSum[2]*temp_diff[15]+0.5*nuVtSqSum[3]*temp_diff[14]+0.5*nuVtSqSum[0]*temp_diff[13]+0.5*nuVtSqSum[1]*temp_diff[10]; 
  diff_incr[14] = 0.5*nuVtSqSum[1]*temp_diff[15]+0.5*nuVtSqSum[0]*temp_diff[14]+0.5*nuVtSqSum[3]*temp_diff[13]+0.5*nuVtSqSum[2]*temp_diff[10]; 
  diff_incr[15] = 0.5*nuVtSqSum[0]*temp_diff[15]+0.5*nuVtSqSum[1]*temp_diff[14]+0.5*nuVtSqSum[2]*temp_diff[13]+0.5*nuVtSqSum[3]*temp_diff[10]; 

  out[0] += diff_incr[0]*rdvSq4; 
  out[1] += diff_incr[1]*rdvSq4; 
  out[2] += diff_incr[2]*rdvSq4; 
  out[3] += diff_incr[3]*rdvSq4; 
  out[4] += diff_incr[4]*rdvSq4; 
  out[5] += diff_incr[5]*rdvSq4; 
  out[6] += diff_incr[6]*rdvSq4; 
  out[7] += diff_incr[7]*rdvSq4; 
  out[8] += diff_incr[8]*rdvSq4; 
  out[9] += diff_incr[9]*rdvSq4; 
  out[10] += diff_incr[10]*rdvSq4; 
  out[11] += diff_incr[11]*rdvSq4; 
  out[12] += diff_incr[12]*rdvSq4; 
  out[13] += diff_incr[13]*rdvSq4; 
  out[14] += diff_incr[14]*rdvSq4; 
  out[15] += diff_incr[15]*rdvSq4; 
} 
