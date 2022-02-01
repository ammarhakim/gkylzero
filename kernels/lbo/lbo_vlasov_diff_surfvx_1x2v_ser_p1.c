#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_diff_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         cell-center coordinates. 
  // dxv[3]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 
  double temp_diff[8] = {0.0}; 
  double diff_incr[8] = {0.0}; 

  temp_diff[0] = (-0.5412658773652741*fr[2])+0.5412658773652741*fl[2]+0.5625*fr[0]+0.5625*fl[0]-1.125*fc[0]; 
  temp_diff[1] = (-0.5412658773652741*fr[4])+0.5412658773652741*fl[4]+0.5625*fr[1]+0.5625*fl[1]-1.125*fc[1]; 
  temp_diff[2] = (-0.4375*fr[2])-0.4375*fl[2]-2.875*fc[2]+0.5412658773652741*fr[0]-0.5412658773652741*fl[0]; 
  temp_diff[3] = (-0.5412658773652741*fr[6])+0.5412658773652741*fl[6]+0.5625*fr[3]+0.5625*fl[3]-1.125*fc[3]; 
  temp_diff[4] = (-0.4375*fr[4])-0.4375*fl[4]-2.875*fc[4]+0.5412658773652741*fr[1]-0.5412658773652741*fl[1]; 
  temp_diff[5] = (-0.5412658773652741*fr[7])+0.5412658773652741*fl[7]+0.5625*fr[5]+0.5625*fl[5]-1.125*fc[5]; 
  temp_diff[6] = (-0.4375*fr[6])-0.4375*fl[6]-2.875*fc[6]+0.5412658773652741*fr[3]-0.5412658773652741*fl[3]; 
  temp_diff[7] = (-0.4375*fr[7])-0.4375*fl[7]-2.875*fc[7]+0.5412658773652741*fr[5]-0.5412658773652741*fl[5]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[4]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[5] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[5]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[6] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[6]; 
  diff_incr[7] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[7]+0.7071067811865475*nuVtSqSum[1]*temp_diff[6]; 

  out[0] += diff_incr[0]*rdvSq4; 
  out[1] += diff_incr[1]*rdvSq4; 
  out[2] += diff_incr[2]*rdvSq4; 
  out[3] += diff_incr[3]*rdvSq4; 
  out[4] += diff_incr[4]*rdvSq4; 
  out[5] += diff_incr[5]*rdvSq4; 
  out[6] += diff_incr[6]*rdvSq4; 
  out[7] += diff_incr[7]*rdvSq4; 
} 
