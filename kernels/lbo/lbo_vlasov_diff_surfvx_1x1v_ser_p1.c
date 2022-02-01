#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_diff_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         cell-center coordinates. 
  // dxv[2]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 
  double temp_diff[4] = {0.0}; 
  double diff_incr[4] = {0.0}; 

  temp_diff[0] = (-0.5412658773652741*fr[2])+0.5412658773652741*fl[2]+0.5625*fr[0]+0.5625*fl[0]-1.125*fc[0]; 
  temp_diff[1] = (-0.5412658773652741*fr[3])+0.5412658773652741*fl[3]+0.5625*fr[1]+0.5625*fl[1]-1.125*fc[1]; 
  temp_diff[2] = (-0.4375*fr[2])-0.4375*fl[2]-2.875*fc[2]+0.5412658773652741*fr[0]-0.5412658773652741*fl[0]; 
  temp_diff[3] = (-0.4375*fr[3])-0.4375*fl[3]-2.875*fc[3]+0.5412658773652741*fr[1]-0.5412658773652741*fl[1]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[3]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[3]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 

  out[0] += diff_incr[0]*rdvSq4; 
  out[1] += diff_incr[1]*rdvSq4; 
  out[2] += diff_incr[2]*rdvSq4; 
  out[3] += diff_incr[3]*rdvSq4; 
} 
