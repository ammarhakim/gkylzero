#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_diff_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]: cell-center coordinates. 
  // dxv[3]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[9]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  const double *nuVtSqSum = &nuPrimMomsSum[6];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 
  double temp_diff[20] = {0.0}; 
  double diff_incr[20] = {0.0}; 

  temp_diff[0] = 0.6708203932499369*fr[9]+0.6708203932499369*fl[9]-1.341640786499874*fc[9]-1.190784930203603*fr[3]+1.190784930203603*fl[3]+0.9375*fr[0]+0.9375*fl[0]-1.875*fc[0]; 
  temp_diff[1] = 0.6708203932499369*fr[15]+0.6708203932499369*fl[15]-1.341640786499874*fc[15]-1.190784930203603*fr[5]+1.190784930203603*fl[5]+0.9375*fr[1]+0.9375*fl[1]-1.875*fc[1]; 
  temp_diff[2] = 0.6708203932499369*fr[16]+0.6708203932499369*fl[16]-1.341640786499874*fc[16]-1.190784930203603*fr[6]+1.190784930203603*fl[6]+0.9375*fr[2]+0.9375*fl[2]-1.875*fc[2]; 
  temp_diff[3] = 0.7382874503707888*fr[9]-0.7382874503707888*fl[9]-1.453125*fr[3]-1.453125*fl[3]-5.34375*fc[3]+1.190784930203603*fr[0]-1.190784930203603*fl[0]; 
  temp_diff[4] = 0.6708203932499369*fr[19]+0.6708203932499369*fl[19]-1.341640786499874*fc[19]-1.190784930203603*fr[10]+1.190784930203603*fl[10]+0.9375*fr[4]+0.9375*fl[4]-1.875*fc[4]; 
  temp_diff[5] = 0.7382874503707888*fr[15]-0.7382874503707888*fl[15]-1.453125*fr[5]-1.453125*fl[5]-5.34375*fc[5]+1.190784930203603*fr[1]-1.190784930203603*fl[1]; 
  temp_diff[6] = 0.7382874503707888*fr[16]-0.7382874503707888*fl[16]-1.453125*fr[6]-1.453125*fl[6]-5.34375*fc[6]+1.190784930203603*fr[2]-1.190784930203603*fl[2]; 
  temp_diff[7] = (-1.190784930203603*fr[13])+1.190784930203603*fl[13]+0.9375*fr[7]+0.9375*fl[7]-1.875*fc[7]; 
  temp_diff[8] = (-1.190784930203603*fr[14])+1.190784930203603*fl[14]+0.9375*fr[8]+0.9375*fl[8]-1.875*fc[8]; 
  temp_diff[9] = (-0.140625*fr[9])-0.140625*fl[9]-6.28125*fc[9]-0.3025768239224545*fr[3]+0.3025768239224545*fl[3]+0.4192627457812106*fr[0]+0.4192627457812106*fl[0]-0.8385254915624212*fc[0]; 
  temp_diff[10] = 0.7382874503707888*fr[19]-0.7382874503707888*fl[19]-1.453125*fr[10]-1.453125*fl[10]-5.34375*fc[10]+1.190784930203603*fr[4]-1.190784930203603*fl[4]; 
  temp_diff[11] = (-1.190784930203603*fr[17])+1.190784930203603*fl[17]+0.9375*fr[11]+0.9375*fl[11]-1.875*fc[11]; 
  temp_diff[12] = (-1.190784930203603*fr[18])+1.190784930203603*fl[18]+0.9375*fr[12]+0.9375*fl[12]-1.875*fc[12]; 
  temp_diff[13] = (-1.453125*fr[13])-1.453125*fl[13]-5.34375*fc[13]+1.190784930203603*fr[7]-1.190784930203603*fl[7]; 
  temp_diff[14] = (-1.453125*fr[14])-1.453125*fl[14]-5.34375*fc[14]+1.190784930203603*fr[8]-1.190784930203603*fl[8]; 
  temp_diff[15] = (-0.140625*fr[15])-0.140625*fl[15]-6.28125*fc[15]-0.3025768239224544*fr[5]+0.3025768239224544*fl[5]+0.4192627457812105*fr[1]+0.4192627457812105*fl[1]-0.8385254915624211*fc[1]; 
  temp_diff[16] = (-0.140625*fr[16])-0.140625*fl[16]-6.28125*fc[16]-0.3025768239224544*fr[6]+0.3025768239224544*fl[6]+0.4192627457812105*fr[2]+0.4192627457812105*fl[2]-0.8385254915624211*fc[2]; 
  temp_diff[17] = (-1.453125*fr[17])-1.453125*fl[17]-5.34375*fc[17]+1.190784930203603*fr[11]-1.190784930203603*fl[11]; 
  temp_diff[18] = (-1.453125*fr[18])-1.453125*fl[18]-5.34375*fc[18]+1.190784930203603*fr[12]-1.190784930203603*fl[12]; 
  temp_diff[19] = (-0.140625*fr[19])-0.140625*fl[19]-6.28125*fc[19]-0.3025768239224545*fr[10]+0.3025768239224545*fl[10]+0.4192627457812106*fr[4]+0.4192627457812106*fl[4]-0.8385254915624212*fc[4]; 

  diff_incr[0] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[7]+0.7071067811865475*nuVtSqSum[1]*temp_diff[1]+0.7071067811865475*nuVtSqSum[0]*temp_diff[0]; 
  diff_incr[1] = 0.6324555320336759*nuVtSqSum[1]*temp_diff[7]+0.6324555320336759*temp_diff[1]*nuVtSqSum[2]+0.7071067811865475*nuVtSqSum[0]*temp_diff[1]+0.7071067811865475*temp_diff[0]*nuVtSqSum[1]; 
  diff_incr[2] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[11]+0.7071067811865475*nuVtSqSum[1]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[2]; 
  diff_incr[3] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[13]+0.7071067811865475*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[3]; 
  diff_incr[4] = 0.632455532033676*nuVtSqSum[1]*temp_diff[11]+0.6324555320336759*nuVtSqSum[2]*temp_diff[4]+0.7071067811865475*nuVtSqSum[0]*temp_diff[4]+0.7071067811865475*nuVtSqSum[1]*temp_diff[2]; 
  diff_incr[5] = 0.632455532033676*nuVtSqSum[1]*temp_diff[13]+0.6324555320336759*nuVtSqSum[2]*temp_diff[5]+0.7071067811865475*nuVtSqSum[0]*temp_diff[5]+0.7071067811865475*nuVtSqSum[1]*temp_diff[3]; 
  diff_incr[6] = 0.7071067811865475*nuVtSqSum[2]*temp_diff[17]+0.7071067811865475*nuVtSqSum[1]*temp_diff[10]+0.7071067811865475*nuVtSqSum[0]*temp_diff[6]; 
  diff_incr[7] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[7]+0.7071067811865475*nuVtSqSum[0]*temp_diff[7]+0.7071067811865475*temp_diff[0]*nuVtSqSum[2]+0.6324555320336759*nuVtSqSum[1]*temp_diff[1]; 
  diff_incr[8] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[12]+0.7071067811865475*nuVtSqSum[0]*temp_diff[8]; 
  diff_incr[9] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[15]+0.7071067811865475*nuVtSqSum[0]*temp_diff[9]; 
  diff_incr[10] = 0.6324555320336759*nuVtSqSum[1]*temp_diff[17]+0.6324555320336759*nuVtSqSum[2]*temp_diff[10]+0.7071067811865475*nuVtSqSum[0]*temp_diff[10]+0.7071067811865475*nuVtSqSum[1]*temp_diff[6]; 
  diff_incr[11] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[11]+0.7071067811865475*nuVtSqSum[0]*temp_diff[11]+0.632455532033676*nuVtSqSum[1]*temp_diff[4]+0.7071067811865475*nuVtSqSum[2]*temp_diff[2]; 
  diff_incr[12] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[12]+0.7071067811865475*nuVtSqSum[0]*temp_diff[12]+0.7071067811865475*nuVtSqSum[1]*temp_diff[8]; 
  diff_incr[13] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[13]+0.7071067811865475*nuVtSqSum[0]*temp_diff[13]+0.632455532033676*nuVtSqSum[1]*temp_diff[5]+0.7071067811865475*nuVtSqSum[2]*temp_diff[3]; 
  diff_incr[14] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[18]+0.7071067811865475*nuVtSqSum[0]*temp_diff[14]; 
  diff_incr[15] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[15]+0.7071067811865475*nuVtSqSum[0]*temp_diff[15]+0.7071067811865475*nuVtSqSum[1]*temp_diff[9]; 
  diff_incr[16] = 0.7071067811865475*nuVtSqSum[1]*temp_diff[19]+0.7071067811865475*nuVtSqSum[0]*temp_diff[16]; 
  diff_incr[17] = 0.4517539514526256*nuVtSqSum[2]*temp_diff[17]+0.7071067811865475*nuVtSqSum[0]*temp_diff[17]+0.6324555320336759*nuVtSqSum[1]*temp_diff[10]+0.7071067811865475*nuVtSqSum[2]*temp_diff[6]; 
  diff_incr[18] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[18]+0.7071067811865475*nuVtSqSum[0]*temp_diff[18]+0.7071067811865475*nuVtSqSum[1]*temp_diff[14]; 
  diff_incr[19] = 0.6324555320336759*nuVtSqSum[2]*temp_diff[19]+0.7071067811865475*nuVtSqSum[0]*temp_diff[19]+0.7071067811865475*nuVtSqSum[1]*temp_diff[16]; 

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
} 
