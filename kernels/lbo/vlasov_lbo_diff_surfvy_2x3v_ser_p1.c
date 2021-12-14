#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_diff_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[5]:         cell-center coordinates. 
  // dxv[5]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[12]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 
  double temp_diff[32] = {0.0}; 
  double diff_incr[32] = {0.0}; 

  temp_diff[0] = (-0.5412658773652741*fr[4])+0.5412658773652741*fl[4]+0.5625*fr[0]+0.5625*fl[0]-1.125*fc[0]; 
  temp_diff[1] = (-0.5412658773652741*fr[9])+0.5412658773652741*fl[9]+0.5625*fr[1]+0.5625*fl[1]-1.125*fc[1]; 
  temp_diff[2] = (-0.5412658773652741*fr[10])+0.5412658773652741*fl[10]+0.5625*fr[2]+0.5625*fl[2]-1.125*fc[2]; 
  temp_diff[3] = (-0.5412658773652741*fr[11])+0.5412658773652741*fl[11]+0.5625*fr[3]+0.5625*fl[3]-1.125*fc[3]; 
  temp_diff[4] = (-0.4375*fr[4])-0.4375*fl[4]-2.875*fc[4]+0.5412658773652741*fr[0]-0.5412658773652741*fl[0]; 
  temp_diff[5] = (-0.5412658773652741*fr[15])+0.5412658773652741*fl[15]+0.5625*fr[5]+0.5625*fl[5]-1.125*fc[5]; 
  temp_diff[6] = (-0.5412658773652741*fr[17])+0.5412658773652741*fl[17]+0.5625*fr[6]+0.5625*fl[6]-1.125*fc[6]; 
  temp_diff[7] = (-0.5412658773652741*fr[18])+0.5412658773652741*fl[18]+0.5625*fr[7]+0.5625*fl[7]-1.125*fc[7]; 
  temp_diff[8] = (-0.5412658773652741*fr[19])+0.5412658773652741*fl[19]+0.5625*fr[8]+0.5625*fl[8]-1.125*fc[8]; 
  temp_diff[9] = (-0.4375*fr[9])-0.4375*fl[9]-2.875*fc[9]+0.5412658773652741*fr[1]-0.5412658773652741*fl[1]; 
  temp_diff[10] = (-0.4375*fr[10])-0.4375*fl[10]-2.875*fc[10]+0.5412658773652741*fr[2]-0.5412658773652741*fl[2]; 
  temp_diff[11] = (-0.4375*fr[11])-0.4375*fl[11]-2.875*fc[11]+0.5412658773652741*fr[3]-0.5412658773652741*fl[3]; 
  temp_diff[12] = (-0.5412658773652741*fr[23])+0.5412658773652741*fl[23]+0.5625*fr[12]+0.5625*fl[12]-1.125*fc[12]; 
  temp_diff[13] = (-0.5412658773652741*fr[24])+0.5412658773652741*fl[24]+0.5625*fr[13]+0.5625*fl[13]-1.125*fc[13]; 
  temp_diff[14] = (-0.5412658773652741*fr[25])+0.5412658773652741*fl[25]+0.5625*fr[14]+0.5625*fl[14]-1.125*fc[14]; 
  temp_diff[15] = (-0.4375*fr[15])-0.4375*fl[15]-2.875*fc[15]+0.5412658773652741*fr[5]-0.5412658773652741*fl[5]; 
  temp_diff[16] = (-0.5412658773652741*fr[26])+0.5412658773652741*fl[26]+0.5625*fr[16]+0.5625*fl[16]-1.125*fc[16]; 
  temp_diff[17] = (-0.4375*fr[17])-0.4375*fl[17]-2.875*fc[17]+0.5412658773652741*fr[6]-0.5412658773652741*fl[6]; 
  temp_diff[18] = (-0.4375*fr[18])-0.4375*fl[18]-2.875*fc[18]+0.5412658773652741*fr[7]-0.5412658773652741*fl[7]; 
  temp_diff[19] = (-0.4375*fr[19])-0.4375*fl[19]-2.875*fc[19]+0.5412658773652741*fr[8]-0.5412658773652741*fl[8]; 
  temp_diff[20] = (-0.5412658773652741*fr[28])+0.5412658773652741*fl[28]+0.5625*fr[20]+0.5625*fl[20]-1.125*fc[20]; 
  temp_diff[21] = (-0.5412658773652741*fr[29])+0.5412658773652741*fl[29]+0.5625*fr[21]+0.5625*fl[21]-1.125*fc[21]; 
  temp_diff[22] = (-0.5412658773652741*fr[30])+0.5412658773652741*fl[30]+0.5625*fr[22]+0.5625*fl[22]-1.125*fc[22]; 
  temp_diff[23] = (-0.4375*fr[23])-0.4375*fl[23]-2.875*fc[23]+0.5412658773652741*fr[12]-0.5412658773652741*fl[12]; 
  temp_diff[24] = (-0.4375*fr[24])-0.4375*fl[24]-2.875*fc[24]+0.5412658773652741*fr[13]-0.5412658773652741*fl[13]; 
  temp_diff[25] = (-0.4375*fr[25])-0.4375*fl[25]-2.875*fc[25]+0.5412658773652741*fr[14]-0.5412658773652741*fl[14]; 
  temp_diff[26] = (-0.4375*fr[26])-0.4375*fl[26]-2.875*fc[26]+0.5412658773652741*fr[16]-0.5412658773652741*fl[16]; 
  temp_diff[27] = (-0.5412658773652741*fr[31])+0.5412658773652741*fl[31]+0.5625*fr[27]+0.5625*fl[27]-1.125*fc[27]; 
  temp_diff[28] = (-0.4375*fr[28])-0.4375*fl[28]-2.875*fc[28]+0.5412658773652741*fr[20]-0.5412658773652741*fl[20]; 
  temp_diff[29] = (-0.4375*fr[29])-0.4375*fl[29]-2.875*fc[29]+0.5412658773652741*fr[21]-0.5412658773652741*fl[21]; 
  temp_diff[30] = (-0.4375*fr[30])-0.4375*fl[30]-2.875*fc[30]+0.5412658773652741*fr[22]-0.5412658773652741*fl[22]; 
  temp_diff[31] = (-0.4375*fr[31])-0.4375*fl[31]-2.875*fc[31]+0.5412658773652741*fr[27]-0.5412658773652741*fl[27]; 

  diff_incr[0] = 0.5*nuVtSqSum[3]*temp_diff[6]+0.5*nuVtSqSum[2]*temp_diff[2]+0.5*nuVtSqSum[1]*temp_diff[1]+0.5*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.5*nuVtSqSum[2]*temp_diff[6]+0.5*temp_diff[2]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[1]+0.5*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.5*nuVtSqSum[1]*temp_diff[6]+0.5*temp_diff[1]*nuVtSqSum[3]+0.5*nuVtSqSum[0]*temp_diff[2]+0.5*temp_diff[0]*nuVtSqSum[2]; 
  diff_incr[3] = 0.5*nuVtSqSum[3]*temp_diff[16]+0.5*nuVtSqSum[2]*temp_diff[8]+0.5*nuVtSqSum[1]*temp_diff[7]+0.5*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.5*nuVtSqSum[3]*temp_diff[17]+0.5*nuVtSqSum[2]*temp_diff[10]+0.5*nuVtSqSum[1]*temp_diff[9]+0.5*nuVtSqSum[0]*temp_diff[4]; 
  diff_incr[5] = 0.5*nuVtSqSum[3]*temp_diff[20]+0.5*nuVtSqSum[2]*temp_diff[13]+0.5*nuVtSqSum[1]*temp_diff[12]+0.5*nuVtSqSum[0]*temp_diff[5]; 
  diff_incr[6] = 0.5*nuVtSqSum[0]*temp_diff[6]+0.5*temp_diff[0]*nuVtSqSum[3]+0.5*nuVtSqSum[1]*temp_diff[2]+0.5*temp_diff[1]*nuVtSqSum[2]; 
  diff_incr[7] = 0.5*nuVtSqSum[2]*temp_diff[16]+0.5*nuVtSqSum[3]*temp_diff[8]+0.5*nuVtSqSum[0]*temp_diff[7]+0.5*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[8] = 0.5*nuVtSqSum[1]*temp_diff[16]+0.5*nuVtSqSum[0]*temp_diff[8]+0.5*nuVtSqSum[3]*temp_diff[7]+0.5*nuVtSqSum[2]*temp_diff[3]; 
  diff_incr[9] = 0.5*nuVtSqSum[2]*temp_diff[17]+0.5*nuVtSqSum[3]*temp_diff[10]+0.5*nuVtSqSum[0]*temp_diff[9]+0.5*nuVtSqSum[1]*temp_diff[4]; 
  diff_incr[10] = 0.5*nuVtSqSum[1]*temp_diff[17]+0.5*nuVtSqSum[0]*temp_diff[10]+0.5*nuVtSqSum[3]*temp_diff[9]+0.5*nuVtSqSum[2]*temp_diff[4]; 
  diff_incr[11] = 0.5*nuVtSqSum[3]*temp_diff[26]+0.5*nuVtSqSum[2]*temp_diff[19]+0.5*nuVtSqSum[1]*temp_diff[18]+0.5*nuVtSqSum[0]*temp_diff[11]; 
  diff_incr[12] = 0.5*nuVtSqSum[2]*temp_diff[20]+0.5*nuVtSqSum[3]*temp_diff[13]+0.5*nuVtSqSum[0]*temp_diff[12]+0.5*nuVtSqSum[1]*temp_diff[5]; 
  diff_incr[13] = 0.5*nuVtSqSum[1]*temp_diff[20]+0.5*nuVtSqSum[0]*temp_diff[13]+0.5*nuVtSqSum[3]*temp_diff[12]+0.5*nuVtSqSum[2]*temp_diff[5]; 
  diff_incr[14] = 0.5*nuVtSqSum[3]*temp_diff[27]+0.5*nuVtSqSum[2]*temp_diff[22]+0.5*nuVtSqSum[1]*temp_diff[21]+0.5*nuVtSqSum[0]*temp_diff[14]; 
  diff_incr[15] = 0.5*nuVtSqSum[3]*temp_diff[28]+0.5*nuVtSqSum[2]*temp_diff[24]+0.5*nuVtSqSum[1]*temp_diff[23]+0.5*nuVtSqSum[0]*temp_diff[15]; 
  diff_incr[16] = 0.5*nuVtSqSum[0]*temp_diff[16]+0.5*nuVtSqSum[1]*temp_diff[8]+0.5*nuVtSqSum[2]*temp_diff[7]+0.5*nuVtSqSum[3]*temp_diff[3]; 
  diff_incr[17] = 0.5*nuVtSqSum[0]*temp_diff[17]+0.5*nuVtSqSum[1]*temp_diff[10]+0.5*nuVtSqSum[2]*temp_diff[9]+0.5*nuVtSqSum[3]*temp_diff[4]; 
  diff_incr[18] = 0.5*nuVtSqSum[2]*temp_diff[26]+0.5*nuVtSqSum[3]*temp_diff[19]+0.5*nuVtSqSum[0]*temp_diff[18]+0.5*nuVtSqSum[1]*temp_diff[11]; 
  diff_incr[19] = 0.5*nuVtSqSum[1]*temp_diff[26]+0.5*nuVtSqSum[0]*temp_diff[19]+0.5*nuVtSqSum[3]*temp_diff[18]+0.5*nuVtSqSum[2]*temp_diff[11]; 
  diff_incr[20] = 0.5*nuVtSqSum[0]*temp_diff[20]+0.5*nuVtSqSum[1]*temp_diff[13]+0.5*nuVtSqSum[2]*temp_diff[12]+0.5*nuVtSqSum[3]*temp_diff[5]; 
  diff_incr[21] = 0.5*nuVtSqSum[2]*temp_diff[27]+0.5*nuVtSqSum[3]*temp_diff[22]+0.5*nuVtSqSum[0]*temp_diff[21]+0.5*nuVtSqSum[1]*temp_diff[14]; 
  diff_incr[22] = 0.5*nuVtSqSum[1]*temp_diff[27]+0.5*nuVtSqSum[0]*temp_diff[22]+0.5*nuVtSqSum[3]*temp_diff[21]+0.5*nuVtSqSum[2]*temp_diff[14]; 
  diff_incr[23] = 0.5*nuVtSqSum[2]*temp_diff[28]+0.5*nuVtSqSum[3]*temp_diff[24]+0.5*nuVtSqSum[0]*temp_diff[23]+0.5*nuVtSqSum[1]*temp_diff[15]; 
  diff_incr[24] = 0.5*nuVtSqSum[1]*temp_diff[28]+0.5*nuVtSqSum[0]*temp_diff[24]+0.5*nuVtSqSum[3]*temp_diff[23]+0.5*nuVtSqSum[2]*temp_diff[15]; 
  diff_incr[25] = 0.5*nuVtSqSum[3]*temp_diff[31]+0.5*nuVtSqSum[2]*temp_diff[30]+0.5*nuVtSqSum[1]*temp_diff[29]+0.5*nuVtSqSum[0]*temp_diff[25]; 
  diff_incr[26] = 0.5*nuVtSqSum[0]*temp_diff[26]+0.5*nuVtSqSum[1]*temp_diff[19]+0.5*nuVtSqSum[2]*temp_diff[18]+0.5*nuVtSqSum[3]*temp_diff[11]; 
  diff_incr[27] = 0.5*nuVtSqSum[0]*temp_diff[27]+0.5*nuVtSqSum[1]*temp_diff[22]+0.5*nuVtSqSum[2]*temp_diff[21]+0.5*nuVtSqSum[3]*temp_diff[14]; 
  diff_incr[28] = 0.5*nuVtSqSum[0]*temp_diff[28]+0.5*nuVtSqSum[1]*temp_diff[24]+0.5*nuVtSqSum[2]*temp_diff[23]+0.5*nuVtSqSum[3]*temp_diff[15]; 
  diff_incr[29] = 0.5*nuVtSqSum[2]*temp_diff[31]+0.5*nuVtSqSum[3]*temp_diff[30]+0.5*nuVtSqSum[0]*temp_diff[29]+0.5*nuVtSqSum[1]*temp_diff[25]; 
  diff_incr[30] = 0.5*nuVtSqSum[1]*temp_diff[31]+0.5*nuVtSqSum[0]*temp_diff[30]+0.5*nuVtSqSum[3]*temp_diff[29]+0.5*nuVtSqSum[2]*temp_diff[25]; 
  diff_incr[31] = 0.5*nuVtSqSum[0]*temp_diff[31]+0.5*nuVtSqSum[1]*temp_diff[30]+0.5*nuVtSqSum[2]*temp_diff[29]+0.5*nuVtSqSum[3]*temp_diff[25]; 

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
  out[16] += diff_incr[16]*rdvSq4; 
  out[17] += diff_incr[17]*rdvSq4; 
  out[18] += diff_incr[18]*rdvSq4; 
  out[19] += diff_incr[19]*rdvSq4; 
  out[20] += diff_incr[20]*rdvSq4; 
  out[21] += diff_incr[21]*rdvSq4; 
  out[22] += diff_incr[22]*rdvSq4; 
  out[23] += diff_incr[23]*rdvSq4; 
  out[24] += diff_incr[24]*rdvSq4; 
  out[25] += diff_incr[25]*rdvSq4; 
  out[26] += diff_incr[26]*rdvSq4; 
  out[27] += diff_incr[27]*rdvSq4; 
  out[28] += diff_incr[28]*rdvSq4; 
  out[29] += diff_incr[29]*rdvSq4; 
  out[30] += diff_incr[30]*rdvSq4; 
  out[31] += diff_incr[31]*rdvSq4; 
} 
