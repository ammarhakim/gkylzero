#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_diff_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         cell-center coordinates. 
  // dxv[3]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 
  double temp_diff[16] = {0.0}; 
  double diff_incr[16] = {0.0}; 

  temp_diff[0] = 0.6708203932499369*fr[12]+0.6708203932499369*fl[12]-1.341640786499874*fc[12]-1.190784930203603*fr[3]+1.190784930203603*fl[3]+0.9375*fr[0]+0.9375*fl[0]-1.875*fc[0]; 
  temp_diff[1] = 0.6708203932499369*fr[13]+0.6708203932499369*fl[13]-1.341640786499874*fc[13]-1.190784930203603*fr[5]+1.190784930203603*fl[5]+0.9375*fr[1]+0.9375*fl[1]-1.875*fc[1]; 
  temp_diff[2] = 0.6708203932499369*fr[14]+0.6708203932499369*fl[14]-1.341640786499874*fc[14]-1.190784930203603*fr[6]+1.190784930203603*fl[6]+0.9375*fr[2]+0.9375*fl[2]-1.875*fc[2]; 
  temp_diff[3] = 0.7382874503707888*fr[12]-0.7382874503707888*fl[12]-1.453125*fr[3]-1.453125*fl[3]-5.34375*fc[3]+1.190784930203603*fr[0]-1.190784930203603*fl[0]; 
  temp_diff[4] = 0.6708203932499369*fr[15]+0.6708203932499369*fl[15]-1.341640786499874*fc[15]-1.190784930203603*fr[7]+1.190784930203603*fl[7]+0.9375*fr[4]+0.9375*fl[4]-1.875*fc[4]; 
  temp_diff[5] = 0.7382874503707888*fr[13]-0.7382874503707888*fl[13]-1.453125*fr[5]-1.453125*fl[5]-5.34375*fc[5]+1.190784930203603*fr[1]-1.190784930203603*fl[1]; 
  temp_diff[6] = 0.7382874503707888*fr[14]-0.7382874503707888*fl[14]-1.453125*fr[6]-1.453125*fl[6]-5.34375*fc[6]+1.190784930203603*fr[2]-1.190784930203603*fl[2]; 
  temp_diff[7] = 0.7382874503707888*fr[15]-0.7382874503707888*fl[15]-1.453125*fr[7]-1.453125*fl[7]-5.34375*fc[7]+1.190784930203603*fr[4]-1.190784930203603*fl[4]; 
  temp_diff[8] = (-1.190784930203603*fr[10])+1.190784930203603*fl[10]+0.9375*fr[8]+0.9375*fl[8]-1.875*fc[8]; 
  temp_diff[9] = (-1.190784930203603*fr[11])+1.190784930203603*fl[11]+0.9375*fr[9]+0.9375*fl[9]-1.875*fc[9]; 
  temp_diff[10] = (-1.453125*fr[10])-1.453125*fl[10]-5.34375*fc[10]+1.190784930203603*fr[8]-1.190784930203603*fl[8]; 
  temp_diff[11] = (-1.453125*fr[11])-1.453125*fl[11]-5.34375*fc[11]+1.190784930203603*fr[9]-1.190784930203603*fl[9]; 
  temp_diff[12] = (-0.140625*fr[12])-0.140625*fl[12]-6.28125*fc[12]-0.3025768239224545*fr[3]+0.3025768239224545*fl[3]+0.4192627457812106*fr[0]+0.4192627457812106*fl[0]-0.8385254915624212*fc[0]; 
  temp_diff[13] = (-0.140625*fr[13])-0.140625*fl[13]-6.28125*fc[13]-0.3025768239224544*fr[5]+0.3025768239224544*fl[5]+0.4192627457812105*fr[1]+0.4192627457812105*fl[1]-0.8385254915624211*fc[1]; 
  temp_diff[14] = (-0.140625*fr[14])-0.140625*fl[14]-6.28125*fc[14]-0.3025768239224544*fr[6]+0.3025768239224544*fl[6]+0.4192627457812105*fr[2]+0.4192627457812105*fl[2]-0.8385254915624211*fc[2]; 
  temp_diff[15] = (-0.140625*fr[15])-0.140625*fl[15]-6.28125*fc[15]-0.3025768239224545*fr[7]+0.3025768239224545*fl[7]+0.4192627457812106*fr[4]+0.4192627457812106*fl[4]-0.8385254915624212*fc[4]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[4]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[5] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[5]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[6] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[6]; 
  diff_incr[7] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[7]+0.7071067811865475*nuVtSqSum[1]*temp_diff[6]; 
  diff_incr[8] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[9]+0.7071067811865475*nuVtSqSum[0]*temp_diff[8]; 
  diff_incr[9] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[9]+0.7071067811865475*nuVtSqSum[1]*temp_diff[8]; 
  diff_incr[10] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[11]+0.7071067811865475*nuVtSqSum[0]*temp_diff[10]; 
  diff_incr[11] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[11]+0.7071067811865475*nuVtSqSum[1]*temp_diff[10]; 
  diff_incr[12] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[13]+0.7071067811865475*nuVtSqSum[0]*temp_diff[12]; 
  diff_incr[13] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[13]+0.7071067811865475*nuVtSqSum[1]*temp_diff[12]; 
  diff_incr[14] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[15]+0.7071067811865475*nuVtSqSum[0]*temp_diff[14]; 
  diff_incr[15] = 0.7071067811865475*nuVtSqSum[0]*temp_diff[15]+0.7071067811865475*nuVtSqSum[1]*temp_diff[14]; 

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
